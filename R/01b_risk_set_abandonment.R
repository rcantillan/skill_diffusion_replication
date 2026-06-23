# ==============================================================================
# 01b_risk_set_abandonment.R
#
# Builds the dyadic ABANDONMENT risk set.
#
# LOGIC:
#   (i, j, s) | RCA_i(s,t0) > 1 & RCA_j(s,t0) > 1   (BOTH specialized at t0)
#   Outcome: does occupation j fall below specialization in skill s by t1?
#
# Mirrors 01a_risk_set_adoption.R exactly except for Step 4 (risk-set
# construction): adoption draws sources and targets from DISJOINT pools
# (RCA>1 vs RCA<=1); abandonment draws BOTH from the SAME pool (RCA>1),
# with i != j enforced in the chunked loop as before. The outcome direction
# is reversed: adoption asks "did the target cross ABOVE the threshold,"
# abandonment asks "did the target fall BELOW it."
#
# FIXES (identical to 01a — see that file for full rationale):
#   Fix 1: clean_id collapses O*NET variants (.01,.02 → parent 6-digit SOC code)
#   Fix 2: risk set restricted to intersect(occ_2015, occ_2024)
#   Fix 3: cosine distance over binary RCA>1 with symmetric lookup
#   Fix 4: RCA denominator-drift robustness — carries abandonment_fixed_denom
#          (RCA with b_s frozen at 2015) and delta_importance (raw O*NET
#          importance change) as parallel outcomes alongside the standard
#          RCA threshold.
#   Fix 5: skill Element.Name harmonization across O*NET releases (same
#          8-entry map as 01a; CONFIRMED against O*NET's official Content
#          Model documentation for the Administrative/Clerical pair).
#
# Input:  data/raw/onet/db_15_1/
#         data/raw/onet/db_29_2_text/
#         data/crosswalk/2010_to_2019_Crosswalk.csv
#
# Output: data/derived/all_events_abandonment.rds
#         Columns: source, target, skill_name, abandonment,
#                  abandonment_fixed_denom, delta_importance,
#                  domain, structural_distance, zero_overlap
#
# Intermediates: deleted before final save
# Next:          02b_enrich_abandonment.R
# ==============================================================================

gc(); gc()
library(data.table)
library(Matrix)
library(progress)
library(igraph)
library(lsa)

# ==============================================================================
# Paths
# ==============================================================================
if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

path_2015 <- if (exists("PATH_ONET_2015")) PATH_ONET_2015 else
             "data/raw/onet/db_15_1"
path_2024 <- if (exists("PATH_ONET_2024")) PATH_ONET_2024 else
             "data/raw/onet/db_29_2_text"
path_cw   <- "data/crosswalk/2010_to_2019_Crosswalk.csv"

output_data_dir <- "data/derived"
dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot("O*NET 2015 not found" = dir.exists(path_2015))
stopifnot("O*NET 2024 not found" = dir.exists(path_2024))
stopifnot("Crosswalk not found"  = file.exists(path_cw))

message("Paths verified:")
message("  O*NET 2015: ", path_2015)
message("  O*NET 2024: ", path_2024)
message("  Crosswalk:  ", path_cw)

onet_files <- c("Skills.txt", "Abilities.txt", "Knowledge.txt",
                "Work Activities.txt")

# ==============================================================================
# Fix 1 — clean_id
# ==============================================================================
clean_id <- function(x) {
  s <- as.character(x)
  s <- gsub("-", "", s)
  s <- gsub("\\.[0-9]+$", "", s)
  trimws(s)
}
stopifnot(clean_id("15-1212.00") == "151212")
stopifnot(clean_id("15-1212.01") == "151212")
stopifnot(clean_id("11-3071.04") == "113071")
stopifnot(nchar(clean_id("29-1141.01")) == 6)
message("[OK] clean_id verified")

# ==============================================================================
# Fix 5 — skill Element.Name harmonization across O*NET releases
#
# Identical map to 01a_risk_set_adoption.R. If you edit one, edit both (or
# better, factor this out into a shared R/utils_skill_names.R sourced by
# both scripts, to guarantee they never drift apart).
#
# "Administrative" / "Clerical": CONFIRMED via O*NET's official Content
# Model documentation (onetcenter.org/content.html, cross-checked against
# the pre-revision definition archived by LMI For All, 2019) -- both
# describe the identical construct under a relabeled element name.
# ==============================================================================
skill_name_harmonization <- c(
  "Operations Monitoring"                                 = "Operation Monitoring",
  "Monitoring Processes, Materials, or Surroundings"      = "Monitor Processes, Materials, or Surroundings",
  "Inspecting Equipment, Structures, or Materials"        = "Inspecting Equipment, Structures, or Material",
  "Judging the Qualities of Objects, Services, or People" = "Judging the Qualities of Things, Services, or People",
  "Working with Computers"                                = "Interacting With Computers",
  "Communicating with People Outside the Organization"    = "Communicating with Persons Outside Organization",
  "Providing Consultation and Advice to Others"           = "Provide Consultation and Advice to Others",
  "Administrative"                                        = "Clerical"  # confirmed
)

harmonize_skill_name <- function(x) {
  fifelse(x %in% names(skill_name_harmonization),
          skill_name_harmonization[x],
          x)
}
stopifnot(harmonize_skill_name("Operations Monitoring") == "Operation Monitoring")
stopifnot(harmonize_skill_name("Some Other Skill") == "Some Other Skill")
message("[OK] harmonize_skill_name verified")

# ==============================================================================
# Helpers
# ==============================================================================
load_onet <- function(folder, filename) {
  f <- file.path(folder, filename)
  if (!file.exists(f)) {
    f_alt <- list.files(folder, pattern = gsub(".txt", "", filename, fixed = TRUE),
                        full.names = TRUE, ignore.case = TRUE)
    if (length(f_alt) > 0) f <- f_alt[1] else {
      message("  [!!] Not found: ", filename)
      return(NULL)
    }
  }
  d <- fread(f, sep = "\t", quote = "",
             na.strings = c("NA", "n/a", "", "*"),
             showProgress = FALSE)
  setnames(d, names(d), make.names(names(d)))
  if ("O.NET.SOC.Code" %in% names(d))
    d[, soc_code := clean_id(O.NET.SOC.Code)]
  val_col <- grep("Data.?Value", names(d), value = TRUE)[1]
  if (!is.na(val_col)) d[, value := as.numeric(get(val_col))]
  d
}

calc_rca <- function(d) {
  agg     <- d[, .(v = mean(value, na.rm = TRUE)),
               by = .(soc_code, Element.Name)]
  occ_sum <- agg[, .(ot = sum(v)), by = soc_code]
  ski_sum <- agg[, .(st = sum(v)), by = Element.Name]
  gt      <- sum(agg$v)
  agg     <- merge(agg, occ_sum, by = "soc_code")
  agg     <- merge(agg, ski_sum, by = "Element.Name")
  agg[, a_js := v / ot]        # occupation's own within-portfolio share
  agg[, b_s  := st / gt]       # economy-wide share — this is what drifts
  agg[, rca  := a_js / b_s]
  agg[, .(soc_code, skill = Element.Name, v, a_js, b_s, rca)]
}

# ==============================================================================
# Step 1 — O*NET 2015 + crosswalk
# ==============================================================================
message("\n>>> Step 1: O*NET 2015...")

dt_15 <- rbindlist(
  lapply(onet_files, function(x) load_onet(path_2015, x)),
  fill = TRUE
)

scale_col <- grep("^Scale.ID$|^Scale_ID$", names(dt_15), value = TRUE)[1]
if (is.na(scale_col))
  scale_col <- grep("Scale", names(dt_15), value = TRUE)[1]
message("  Scale ID column: ", scale_col)

dt_15 <- dt_15[get(scale_col) == "IM" & !is.na(value)]
message(sprintf("  IM rows: %s", format(nrow(dt_15), big.mark = ",")))
message(sprintf("  Occupations (raw): %d", uniqueN(dt_15$soc_code)))
stopifnot(all(nchar(unique(dt_15$soc_code)) == 6))
message("  [OK] 6-digit SOC codes")

cw <- fread(path_cw)
setnames(cw, make.names(names(cw)))
col_2010 <- grep("2010.*Code|Code.*2010", names(cw), value = TRUE)[1]
col_2019 <- grep("2019.*Code|Code.*2019", names(cw), value = TRUE)[1]
message("  Crosswalk cols: ", col_2010, " -> ", col_2019)

cw_clean <- unique(cw[, .(
  soc10 = clean_id(get(col_2010)),
  soc19 = clean_id(get(col_2019))
)])
message(sprintf("  Unique mappings: %d", nrow(cw_clean)))

dt_15 <- merge(dt_15, cw_clean,
               by.x = "soc_code", by.y = "soc10",
               all.x = FALSE, allow.cartesian = TRUE)
dt_15[, soc_code := soc19]
message(sprintf("  Post-crosswalk: %d occupations",
                uniqueN(dt_15$soc_code)))

rca_15 <- calc_rca(dt_15)
message(sprintf("  RCA 2015: %d occ x %d skills",
                uniqueN(rca_15$soc_code), uniqueN(rca_15$skill)))
rm(dt_15, cw, cw_clean); gc()

# ==============================================================================
# Step 2 — O*NET 2024
# ==============================================================================
message("\n>>> Step 2: O*NET 2024...")

dt_24 <- rbindlist(
  lapply(onet_files, function(x) load_onet(path_2024, x)),
  fill = TRUE
)
scale_col_24 <- grep("Scale", names(dt_24), value = TRUE)[1]
dt_24 <- dt_24[get(scale_col_24) == "IM" & !is.na(value)]

n_distinct_before <- uniqueN(dt_24$Element.Name)
dt_24[, Element.Name := harmonize_skill_name(Element.Name)]
n_distinct_after <- uniqueN(dt_24$Element.Name)
message(sprintf("  Fix 5: Element.Name harmonized (%d -> %d distinct labels in 2024)",
                n_distinct_before, n_distinct_after))

rca_24 <- calc_rca(dt_24)
message(sprintf("  RCA 2024: %d occ x %d skills",
                uniqueN(rca_24$soc_code), uniqueN(rca_24$skill)))
rm(dt_24); gc()

orphan_check <- length(setdiff(unique(rca_24$skill), unique(rca_15$skill))) +
                length(setdiff(unique(rca_15$skill), unique(rca_24$skill)))
if (orphan_check > 0) {
  warning(sprintf(
    "  [!!] %d skill names still don't match between 2015 and 2024 after harmonization.",
    orphan_check))
} else {
  message("  [OK] All skill names match exactly between 2015 and 2024 post-harmonization")
}

# ==============================================================================
# Fix 2 — Restrict to occupations present in both years
# ==============================================================================
message("\n>>> Fix 2: Intersection 2015 n 2024...")

in_both <- intersect(unique(rca_15$soc_code), unique(rca_24$soc_code))
message(sprintf("  2015: %d | 2024: %d | Intersection: %d",
                uniqueN(rca_15$soc_code), uniqueN(rca_24$soc_code),
                length(in_both)))

rca_15 <- rca_15[soc_code %in% in_both]
rca_24 <- rca_24[soc_code %in% in_both]
message("  [OK] RCA restricted")

# ==============================================================================
# Fix 4 — Fixed-denominator RCA placebo
# ==============================================================================
message("\n>>> Fix 4: Fixed-denominator RCA placebo...")

b_s_2015 <- unique(rca_15[, .(skill, b_s_2015 = b_s)])

rca_24_fixed <- merge(
  rca_24[, .(soc_code, skill, a_js_2024 = a_js)],
  b_s_2015,
  by = "skill"
)
rca_24_fixed[, rca_fixed_denom := a_js_2024 / b_s_2015]

n_skills_24    <- uniqueN(rca_24$skill)
n_skills_fixed <- uniqueN(rca_24_fixed$skill)
if (n_skills_fixed < n_skills_24) {
  warning(sprintf(
    "  [!!] %d skills present in 2024 lack a 2015 baseline denominator and were dropped from the placebo",
    n_skills_24 - n_skills_fixed))
} else {
  message(sprintf("  [OK] Frozen-denominator RCA built for all %d skills",
                  n_skills_fixed))
}

rca_compare <- merge(
  rca_24[, .(soc_code, skill, rca_standard = rca)],
  rca_24_fixed[, .(soc_code, skill, rca_fixed = rca_fixed_denom)],
  by = c("soc_code", "skill")
)
message(sprintf("  Correlation, standard vs. fixed-denominator RCA (2024): %.4f",
                cor(rca_compare$rca_standard, rca_compare$rca_fixed)))
message(sprintf("  Threshold agreement (both sides of RCA=1): %.4f",
                mean((rca_compare$rca_standard > 1) == (rca_compare$rca_fixed > 1))))
rm(rca_compare); gc()

# ==============================================================================
# Fix 3 — Cosine distance over binary RCA>1
# ==============================================================================
message("\n>>> Fix 3: Binary cosine distance...")

active       <- rca_15[rca > 1, .(soc_code, skill)]
occ_ids      <- sort(unique(active$soc_code))
skill_ids    <- sort(unique(active$skill))
n_occ        <- length(occ_ids)
n_pairs_dir  <- as.integer(n_occ) * (n_occ - 1L)
n_pairs_uniq <- n_pairs_dir / 2L

message(sprintf("  Occupations with RCA>1: %d", n_occ))
message(sprintf("  Expected directed pairs: %s",
                format(n_pairs_dir, big.mark = ",")))

M_bin <- sparseMatrix(
  i    = match(active$soc_code, occ_ids),
  j    = match(active$skill,    skill_ids),
  x    = 1L,
  dims = c(n_occ, length(skill_ids)),
  dimnames = list(occ_ids, skill_ids)
)
sizes        <- as.integer(Matrix::rowSums(M_bin))
names(sizes) <- occ_ids

cooc    <- tcrossprod(M_bin)
cooc_s  <- summary(cooc)
cooc_dt <- data.table(
  i = cooc_s$i, j = cooc_s$j, inter = as.integer(cooc_s$x)
)[i < j]

cooc_dt[, `:=`(
  size_i      = sizes[occ_ids[i]],
  size_j      = sizes[occ_ids[j]]
)]
cooc_dt[, cosine_dist := 1 - inter / (sqrt(size_i) * sqrt(size_j))]
cooc_dt[!is.finite(cosine_dist), cosine_dist := 1]
cooc_dt[, `:=`(soc_i = occ_ids[i], soc_j = occ_ids[j])]

message(sprintf("  Cosine dist: median=%.3f | max=%.3f",
                median(cooc_dt$cosine_dist),
                max(cooc_dt$cosine_dist)))

dist_lookup <- rbind(
  cooc_dt[, .(source = soc_i, target = soc_j,
              structural_distance = cosine_dist, zero_overlap = 0L)],
  cooc_dt[, .(source = soc_j, target = soc_i,
              structural_distance = cosine_dist, zero_overlap = 0L)]
)

n_zero <- n_pairs_uniq - nrow(cooc_dt)
if (n_zero > 0) {
  all_pairs <- CJ(source = occ_ids, target = occ_ids)[source != target]
  setkey(all_pairs, source, target)
  setkey(dist_lookup, source, target)
  missing   <- all_pairs[!dist_lookup]
  missing[, `:=`(structural_distance = 1.0, zero_overlap = 1L)]
  dist_lookup <- rbind(dist_lookup, missing)
  rm(all_pairs, missing)
  message(sprintf("  Zero-overlap pairs added: %d", n_zero * 2L))
}

setkey(dist_lookup, source, target)
stopifnot(nrow(dist_lookup) == n_pairs_dir)
stopifnot(sum(is.na(dist_lookup$structural_distance)) == 0)
message("  [OK] dist_lookup complete")

rm(M_bin, cooc, cooc_s, cooc_dt, active); gc()

# ==============================================================================
# Step 3 — Domain labels (Louvain over skill similarity network)
# ==============================================================================
message("\n>>> Step 3: Domain labels...")

mat_dt <- rca_15[rca > 1, .(soc_code, skill)]
occ_wide <- dcast(mat_dt, skill ~ soc_code,
                  value.var = "skill", fun.aggregate = length)
mat           <- as.matrix(occ_wide[, -1])
rownames(mat) <- occ_wide$skill

sim_skill <- lsa::cosine(t(mat))
sim_skill[is.nan(sim_skill)] <- 0
diag(sim_skill) <- 0

g_skill <- igraph::graph_from_adjacency_matrix(
  sim_skill, mode = "undirected", weighted = TRUE)
g_skill <- igraph::delete_edges(
  g_skill, igraph::E(g_skill)[weight < 0.1])
lou <- igraph::cluster_louvain(g_skill)

skill_classes <- data.table(
  skill  = lou$names,
  domain = fifelse(lou$membership == 1L, "Cognitive", "Physical")
)
message(sprintf("  Domains: %s",
                paste(skill_classes[, .N, by = domain][,
                      paste(domain, N, sep = "=")], collapse = " | ")))

rm(mat, sim_skill, g_skill, lou, occ_wide, mat_dt); gc()

# ==============================================================================
# Step 4 — Panel and ABANDONMENT dyad construction
#
# THIS is the step that differs from 01a. Both source and target are drawn
# from the SAME pool (RCA>1 at t0) -- "both i and j specialized in s at
# baseline" per the Methods. i != j is enforced later in the chunked loop,
# same as adoption. The outcome direction flips: 1 if the target FALLS
# BELOW the threshold by t1 (loses what it shared with the source).
# ==============================================================================
message("\n>>> Step 4: Abandonment risk set...")

panel <- merge(
  rca_15[, .(soc_code, skill, rca_t0 = rca, v_t0 = v)],
  rca_24[, .(soc_code, skill, rca_t1 = rca, v_t1 = v)],
  by  = c("soc_code", "skill"),
  all = TRUE
)
panel <- merge(
  panel,
  rca_24_fixed[, .(soc_code, skill, rca_t1_fixed = rca_fixed_denom)],
  by  = c("soc_code", "skill"),
  all.x = TRUE
)
panel[is.na(rca_t0),       rca_t0       := 0]
panel[is.na(rca_t1),       rca_t1       := 0]
panel[is.na(v_t0),         v_t0         := 0]
panel[is.na(v_t1),         v_t1         := 0]
panel[is.na(rca_t1_fixed), rca_t1_fixed := 0]

n_full_orphan_rows <- panel[v_t0 == 0 & v_t1 == 0, .N]
message(sprintf("  Rows with zero importance both years (sanity check): %d", n_full_orphan_rows))

# Sources: RCA>1 at t0 (positional referent that still holds the skill)
sources_t0 <- panel[rca_t0 > 1, .(source = soc_code, skill)]

# Targets: ALSO RCA>1 at t0 -- same pool as sources. Outcome: did it fall
# BELOW the threshold by t1 (the mirror image of adoption's "did it cross
# above").
targets_risk <- panel[rca_t0 > 1, .(
  target           = soc_code,
  skill,
  v_t0, v_t1,
  rca_t1_tgt       = rca_t1,
  rca_t1_fixed_tgt = rca_t1_fixed
)]
targets_risk[, abandonment             := fifelse(rca_t1_tgt <= 1, 1L, 0L)]
targets_risk[, abandonment_fixed_denom := fifelse(rca_t1_fixed_tgt <= 1, 1L, 0L)]
targets_risk[, delta_importance        := v_t1 - v_t0]

rm(panel, rca_15, rca_24, rca_24_fixed, b_s_2015); gc()

# NOTE: sources_t0 and the targets_risk candidate pool come from the SAME
# filter (rca_t0 > 1), so n_src * n_tgt below over-counts the eventual risk
# set size by however many i==j pairs exist (one per occupation-skill cell);
# those get dropped by `pairs[source != target]` in the chunked loop. This
# is an upper-bound estimate, not the final N.
n_est <- targets_risk[, .(n_tgt = .N), by = skill]
n_est <- merge(n_est,
               sources_t0[, .(n_src = .N), by = skill],
               by = "skill")
n_est[, n_pairs := as.numeric(n_src) * n_tgt]
message(sprintf("  Estimated risk set (upper bound, pre i!=j filter): ~%s triads",
                format(sum(n_est$n_pairs), big.mark = ",")))
message(sprintf("  Expected abandonment rate (standard RCA):       %.4f",
                targets_risk[, mean(abandonment)]))
message(sprintf("  Expected abandonment rate (fixed denominator):  %.4f",
                targets_risk[, mean(abandonment_fixed_denom)]))
message(sprintf("  Mean raw importance delta among at-risk dyads:  %.4f",
                targets_risk[, mean(delta_importance)]))

# ==============================================================================
# Chunked loop — by skill
# ==============================================================================
FLUSH_EVERY <- 20L
chunk_dir   <- file.path(output_data_dir, "_aband_chunks")
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

skills_vec  <- unique(sources_t0$skill)
buffer      <- vector("list", FLUSH_EVERY)
buf_i       <- 0L
chunk_files <- character(0)
chunk_n     <- 0L

flush_buffer <- function(buf, n) {
  chunk_n <<- chunk_n + 1L
  f <- file.path(chunk_dir, sprintf("chunk_%04d.rds", chunk_n))
  saveRDS(rbindlist(buf[seq_len(n)]), f, compress = "gzip")
  chunk_files <<- c(chunk_files, f)
  if (chunk_n %% 10 == 0)
    message(sprintf("    Chunk %d saved", chunk_n))
}

pb <- progress::progress_bar$new(
  format = "  [:bar] :percent eta: :eta | skill :current/:total",
  total  = length(skills_vec), clear = FALSE
)

for (sk in skills_vec) {
  pb$tick()
  src_sk <- sources_t0[skill == sk]
  tgt_sk <- targets_risk[skill == sk]
  if (nrow(src_sk) == 0 || nrow(tgt_sk) == 0) next

  pairs <- as.data.table(expand.grid(
    source = src_sk$source,
    target = tgt_sk$target,
    stringsAsFactors = FALSE
  ))
  pairs <- pairs[source != target]
  pairs <- merge(
    pairs,
    tgt_sk[, .(target, abandonment, abandonment_fixed_denom, delta_importance)],
    by = "target"
  )
  pairs[, skill_name := sk]

  buf_i <- buf_i + 1L
  buffer[[buf_i]] <- pairs
  if (buf_i == FLUSH_EVERY) {
    flush_buffer(buffer, buf_i)
    buffer <- vector("list", FLUSH_EVERY)
    buf_i  <- 0L
    gc()
  }
}
if (buf_i > 0L) flush_buffer(buffer, buf_i)
rm(buffer, sources_t0, targets_risk); gc()

message(sprintf("\n  Binding %d chunks...", length(chunk_files)))
all_events <- rbindlist(lapply(chunk_files, readRDS))
unlink(chunk_dir, recursive = TRUE)
gc()

message(sprintf("  Total abandonment dyads: %s",
                format(nrow(all_events), big.mark = ",")))
message(sprintf("  Abandonment rate (standard):    %.4f",
                mean(all_events$abandonment)))
message(sprintf("  Abandonment rate (fixed denom.): %.4f",
                mean(all_events$abandonment_fixed_denom)))

# ==============================================================================
# Step 5 — Domain + structural distance
# ==============================================================================
message("\n>>> Step 5: Domain + structural distance...")

setkey(skill_classes, skill)
all_events[skill_classes, on = .(skill_name = skill), domain := i.domain]
all_events <- all_events[!is.na(domain)]
message(sprintf("  No domain match (removed): %d",
                nrow(all_events[is.na(domain)])))

setkey(all_events, source, target)
all_events[dist_lookup, on = .(source, target),
           `:=`(structural_distance = i.structural_distance,
                zero_overlap        = i.zero_overlap)]

stopifnot(sum(is.na(all_events$structural_distance)) == 0)

all_events[structural_distance < 0, structural_distance := 0]
all_events[structural_distance > 1, structural_distance := 1]
message("  [OK] Structural distance complete")

rm(dist_lookup, skill_classes); gc()

# ==============================================================================
# Step 6 — Save and report
# ==============================================================================
message("\n>>> Step 6: Saving...")

message(sprintf("  Final dyads:            %s",
                format(nrow(all_events), big.mark = ",")))
message(sprintf("  Unique source occ:      %d",
                uniqueN(all_events$source)))
message(sprintf("  Unique target occ:      %d",
                uniqueN(all_events$target)))
message(sprintf("  Unique skills:          %d",
                uniqueN(all_events$skill_name)))
message(sprintf("  Abandonment rate (standard, RCA threshold):  %.4f",
                mean(all_events$abandonment)))
message(sprintf("  Abandonment rate (fixed-denominator placebo): %.4f",
                mean(all_events$abandonment_fixed_denom)))
message(sprintf("  Correlation, standard vs. fixed-denom outcome: %.4f",
                cor(all_events$abandonment, all_events$abandonment_fixed_denom)))
message(sprintf("  Mean raw importance delta:                  %.4f",
                mean(all_events$delta_importance)))

message("\n  --- By domain ---")
print(all_events[, .(
  rate_standard   = mean(abandonment),
  rate_fixed      = mean(abandonment_fixed_denom),
  mean_raw_delta  = mean(delta_importance),
  n               = .N
), by = domain])

out_path <- file.path(output_data_dir, "all_events_abandonment.rds")
saveRDS(all_events, out_path)
message(sprintf("\n  Saved: %s", out_path))
message(sprintf("  Size:  %.1f MB", file.size(out_path) / 1e6))

message("\n  NOTE: 'abandonment' (standard RCA) remains the primary outcome for")
message("  the main models. 'abandonment_fixed_denom' and 'delta_importance' are")
message("  carried forward for the SI robustness check (RCA denominator-drift)")
message("  and must be propagated through 02b_enrich_abandonment.R and the gravity")
message("  model script as alternative dependent variables.")

rm(all_events); gc()
message("\n>>> 01b_risk_set_abandonment.R complete.")
message("    Next: 02b_enrich_abandonment.R")