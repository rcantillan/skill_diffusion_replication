# ==============================================================================
# 03_nestedness.R
#
# cs_s = (N_obs - mean(N*_s)) / sd(N*_s)   [Hosseinioun et al. 2025]
#
# N_SIM = 5000 iterations (standard per Hosseinioun et al. 2025)
#
# OPTIMAL APPROACH — sequential with checkpoint:
#   - No mclapply -> no fork -> no OOM
#   - Internal operations are pure C-level (crossprod, matrix indexing)
#   - Checkpoint every 10 skills -> auto-resume if interrupted
#   - Additional RAM: ~200MB regardless of hardware
#
# Estimated time: ~3-4 hours sequential on standard hardware
#
# Input:  data/derived/all_events_adoption_enriched.rds
# Output: data/derived/skill_cs_scores.rds
# Next:   03b_nestedness_merge.R
#         03c_nestedness_merge_ab.R
# ==============================================================================

gc(); gc()
library(data.table)

N_SIM          <- 5000L   # Hosseinioun et al. 2025
SEED           <- 42L
CHECKPOINT_DIR <- "data/derived/_nestedness_ckpt"
CHECKPOINT_N   <- 10L     # save every N skills

out_dir    <- "data/derived"
input_file <- file.path(out_dir, "all_events_adoption_enriched.rds")

stopifnot("Input not found — run 02a_enrich_adoption.R first" =
            file.exists(input_file))

dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
message(sprintf("N_SIM: %d | checkpoint every %d skills", N_SIM, CHECKPOINT_N))

# ==============================================================================
# Step 1 — Binary occ x skill matrix (2015 skillscape from SOURCE occupations)
# ==============================================================================
message("\n>>> Step 1: Building binary matrix...")

dt <- readRDS(input_file)
setDT(dt)
message(sprintf("  Rows: %s | Source occ: %d | Skills: %d",
                format(nrow(dt), big.mark = ","),
                uniqueN(dt$source),
                uniqueN(dt$skill_name)))

occ_skill <- unique(dt[, .(soc = source, skill = skill_name)])
rm(dt); gc()

M_wide      <- dcast(occ_skill, soc ~ skill,
                     fun.aggregate = length, fill = 0L)
occ_ids     <- M_wide$soc
M           <- as.matrix(M_wide[, -1, with = FALSE])
rownames(M) <- occ_ids
rm(M_wide, occ_skill); gc()

# Sort by decreasing marginals (required by NODF)
M <- M[order(rowSums(M), decreasing = TRUE),
       order(colSums(M), decreasing = TRUE)]
storage.mode(M) <- "integer"

skill_names <- colnames(M)
n_col       <- ncol(M)
n_row       <- nrow(M)
message(sprintf("  Matrix: %d occ x %d skills | density: %.3f",
                n_row, n_col, mean(M > 0)))
message(sprintf("  Matrix RAM: %.1f MB",
                object.size(M) / 1e6))

# ==============================================================================
# Step 2 — Vectorized precomputation (done once)
# ==============================================================================
message("\n>>> Step 2: Precomputation...")

cs_col <- colSums(M)
rs_row <- rowSums(M)

# Valid pair indices for NODF
col_idx   <- which(outer(cs_col, cs_col, ">"), arr.ind = TRUE)
row_idx   <- which(outer(rs_row, rs_row, ">"), arr.ind = TRUE)
col_i     <- col_idx[, 1L]; col_j <- col_idx[, 2L]
row_i     <- row_idx[, 1L]; row_j <- row_idx[, 2L]
col_denom <- cs_col[col_j]
row_denom <- rs_row[row_j]
n_cp      <- length(col_i)
n_rp      <- length(row_i)
rm(col_idx, row_idx); gc()

message(sprintf("  Column pairs: %s | row pairs: %s",
                format(n_cp, big.mark = ","),
                format(n_rp, big.mark = ",")))

# Full overlap (crossprod: pure C-level, ~1 second)
message("  Computing crossprod...")
OV_col        <- crossprod(M)
OV_row        <- tcrossprod(M)
col_dots_base <- OV_col[cbind(col_i, col_j)]
row_dots_base <- OV_row[cbind(row_i, row_j)]
rm(OV_col, OV_row); gc()

# Submatrices for rank-1 correction
# NOTE: precomputed here — dimensions are n_rp x n_col
# With n_rp ~200K and n_col=160: 200K x 160 x 4 bytes = ~128MB — manageable
message("  Precomputing Mk_i / Mk_j...")
Mk_i <- M[row_i, , drop = FALSE]
Mk_j <- M[row_j, , drop = FALSE]
message(sprintf("  RAM Mk_i + Mk_j: %.1f MB",
                (object.size(Mk_i) + object.size(Mk_j)) / 1e6))
message("  [OK] Precomputation complete")

gc()

# ==============================================================================
# Step 3 — Observed NODF + verification
# ==============================================================================
message("\n>>> Step 3: Observed NODF...")

N_obs <- {
  nodf_c <- if (n_cp > 0) mean(col_dots_base / col_denom) * 100 else 0
  nodf_r <- if (n_rp > 0) mean(row_dots_base / row_denom) * 100 else 0
  (nodf_c + nodf_r) / 2
}
message(sprintf("  Observed NODF: %.4f", N_obs))
stopifnot("NODF outside expected range [10, 60]" = N_obs > 10 & N_obs < 60)
message("  [OK] NODF within expected range")

if (requireNamespace("vegan", quietly = TRUE)) {
  N_v    <- vegan::nestednodf(M, order = FALSE,
                               weighted = FALSE)$statistic[["NODF"]]
  diff_v <- abs(N_obs - N_v)
  message(sprintf("  Vegan verification: %.4f | diff = %.2e %s",
                  N_v, diff_v,
                  ifelse(diff_v < 0.01, "[OK]", "[CHECK]")))
} else {
  message("  [SKIP] vegan not installed")
}

# ==============================================================================
# Step 4 — Checkpoint: resume from prior progress if available
# ==============================================================================
ckpt_file <- file.path(CHECKPOINT_DIR, "cs_progress.rds")
if (file.exists(ckpt_file)) {
  ckpt      <- readRDS(ckpt_file)
  cs_vec    <- ckpt$cs_vec
  start_k   <- ckpt$last_k + 1L
  message(sprintf("\n  [RESUME] Resuming from skill %d/%d",
                  start_k, n_col))
  rm(ckpt)
} else {
  cs_vec  <- rep(NA_real_, n_col)
  start_k <- 1L
  message("\n  [NEW] Starting from skill 1")
}

# ==============================================================================
# Step 5 — Sequential simulation loop
#
# For each skill k:
#   1. Permute column k N_SIM times
#   2. Update only affected dot products (rank-1 correction)
#   3. cs_k = (N_obs - mean(null)) / sd(null)
#
# No fork -> no OOM -> constant RAM throughout the loop
# ==============================================================================
message(sprintf(
  "\n>>> Step 5: Computing cs for skills %d-%d x %d sims...",
  start_k, n_col, N_SIM))

set.seed(SEED + start_k)   # reproducible even when resuming
t_start  <- proc.time()["elapsed"]
t_report <- t_start

for (k in start_k:n_col) {

  # Report every 10 skills
  if ((k - start_k) %% 10 == 0 || k == start_k) {
    elapsed <- proc.time()["elapsed"] - t_start
    k_done  <- k - start_k
    eta     <- if (k_done > 0)
                 elapsed / k_done * (n_col - k + 1L)
               else NA_real_
    message(sprintf(
      "  [%3d/%d] %-42s | %5.1f min | eta %5.1f min",
      k, n_col,
      substr(skill_names[k], 1L, 42L),
      elapsed / 60,
      ifelse(is.na(eta), NA, eta / 60)))
  }

  focal_col <- M[, k]

  # Column pairs involving k
  inv_k   <- which(col_i == k | col_j == k)
  other_k <- ifelse(col_i[inv_k] == k, col_j[inv_k], col_i[inv_k])

  # Column k slices of Mk_i / Mk_j (vectors, not matrices)
  Mk_i_k <- Mk_i[, k]
  Mk_j_k <- Mk_j[, k]

  # N_SIM permutations — replicate calls C-level sample()
  null_N <- replicate(N_SIM, {
    new_col <- sample(focal_col)

    # --- NODF columns ---
    col_dots_k <- col_dots_base
    if (length(inv_k) > 0)
      col_dots_k[inv_k] <- as.numeric(
        new_col %*% M[, other_k, drop = FALSE]
      )
    nodf_c <- if (n_cp > 0) mean(col_dots_k / col_denom) * 100 else 0

    # --- NODF rows (vectorized rank-1 correction) ---
    row_dots_new <- row_dots_base -
                    Mk_i_k * Mk_j_k +
                    new_col[row_i] * new_col[row_j]
    new_rs <- rs_row + (new_col - focal_col)
    valid  <- new_rs[row_i] > new_rs[row_j]
    nodf_r <- if (any(valid))
                mean(row_dots_new[valid] / new_rs[row_j[valid]]) * 100
              else 0

    (nodf_c + nodf_r) / 2
  })

  mu  <- mean(null_N)
  sig <- sd(null_N)
  cs_vec[k] <- if (sig > 0) (N_obs - mu) / sig else NA_real_

  # Checkpoint every CHECKPOINT_N skills
  if (k %% CHECKPOINT_N == 0 || k == n_col) {
    saveRDS(list(cs_vec = cs_vec, last_k = k),
            ckpt_file)
  }
}

elapsed_total <- (proc.time()["elapsed"] - t_start) / 60
message(sprintf("\n  Simulation complete in %.1f min (%.1f hours)",
                elapsed_total, elapsed_total / 60))

rm(Mk_i, Mk_j, M); gc()

# ==============================================================================
# Step 6 — cs score table and save
# ==============================================================================
message("\n>>> Step 6: Saving results...")

skill_stats <- data.table(skill_name = skill_names, cs = cs_vec)
skill_stats[, nested := fifelse(cs > 0, "nested", "un-nested")]

n_na <- sum(is.na(skill_stats$cs))
message(sprintf("  Skills with cs: %d | NAs (sd=0): %d",
                nrow(skill_stats) - n_na, n_na))
message(sprintf("  Nested    (cs>0): %d (%.1f%%)",
                sum(skill_stats$nested == "nested",  na.rm = TRUE),
                100 * mean(skill_stats$nested == "nested", na.rm = TRUE)))
message(sprintf("  Un-nested (cs<=0): %d (%.1f%%)",
                sum(skill_stats$nested == "un-nested", na.rm = TRUE),
                100 * mean(skill_stats$nested == "un-nested", na.rm = TRUE)))
message(sprintf("  cs: min=%.3f | median=%.3f | max=%.3f",
                min(skill_stats$cs,    na.rm = TRUE),
                median(skill_stats$cs, na.rm = TRUE),
                max(skill_stats$cs,    na.rm = TRUE)))

saveRDS(skill_stats, file.path(out_dir, "skill_cs_scores.rds"))
message("  Saved: data/derived/skill_cs_scores.rds")

# Remove checkpoints
unlink(CHECKPOINT_DIR, recursive = TRUE)
message("  [OK] Checkpoints removed")

message("\n>>> 03_nestedness.R complete.")
message("    Next: 03b_nestedness_merge.R")
message("          03c_nestedness_merge_ab.R")
