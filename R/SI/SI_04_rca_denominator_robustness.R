# ==============================================================================
# SI_04_rca_denominator_robustness.R  (v3)
#
# SI Robustness Check S3.0: does the signed status-gap gradient survive when
# the adoption/abandonment outcome is redefined to remove RCA's
# denominator-drift channel?
#
# WHY THIS GOES FIRST: every other robustness check (boundary placebo,
# permutation, RCA threshold sensitivity, alternative status measures,
# sub-period stability) assumes the diffusion/abandonment OUTCOME itself is
# validly constructed, and asks whether the MODEL is robust on top of that.
# This check instead interrogates the outcome's construction: RCA is a
# Balassa-style relative index whose weighted mean across occupations equals
# 1 by construction, so a real increase in a skill's prevalence among
# high-status occupations mechanically raises the bar for everyone else,
# which could push low-status occupations below the RCA=1 threshold with
# zero change in their own behavior. If that mechanism explained a
# meaningful share of the result, every downstream robustness check would be
# robust to the wrong thing.
#
# DESIGN: two panels per flow, matching 05a/05b exactly:
#   Panel A: source + skill FE  (focal unit = source occupation)
#   Panel B: target + skill FE  (focal unit = target occupation)
#
# Three parallel outcome definitions per flow x panel:
#   (1) standard         - RCA>1/<=1 threshold crossing (Tables SM3/SM4)
#   (2) fixed_denom      - RCA with economy-wide term frozen at 2015
#   (3) delta_importance - raw O*NET importance change, z-scored within
#                          each flow separately before estimation, so
#                          coefficients are in SD-outcome units and
#                          directionally comparable to the RCA specifications.
#                          NOTE on sign convention: for abandonment, the
#                          sign of raw-importance coefficients is multiplied
#                          by -1 before comparing to the standard, because
#                          abandonment = decrease in importance (negative
#                          raw change), whereas the RCA outcome codes
#                          abandonment as a positive binary event.
#
# Input:
#   data/derived/riskset_adoption.rds
#   data/derived/riskset_abandonment.rds
#   output/tables/main/occ_status_scores.csv   <- from 04_status_pca.R
#
# Output:
#   output/tables/si/tab_S3_0_rca_denom_robustness.csv
#   output/tables/si/tab_S3_0_rca_denom_robustness_full.csv
# ==============================================================================

library(data.table)
library(fixest)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
detect_col <- function(patterns, cols, label) {
  hit <- cols[grepl(paste(patterns, collapse = "|"), cols, ignore.case = TRUE)]
  if (length(hit) == 0) {
    message(sprintf("  [!!] Could not auto-detect %s — check column names.", label))
    return(NA_character_)
  }
  chosen <- hit[1]
  message(sprintf("  [OK] %s -> '%s'", label, chosen))
  chosen
}

# ------------------------------------------------------------------------------
# Status index
# ------------------------------------------------------------------------------
status_file <- "output/tables/main/occ_status_scores.csv"
stopifnot("Canonical status index not found -- run 04_status_pca.R first" =
            file.exists(status_file))
sigma <- fread(status_file, colClasses = c(occ = "character"))
stopifnot(all(c("occ", "status_pc1") %in% names(sigma)))
stopifnot("occ codes are not 6 characters" = all(nchar(sigma$occ) == 6))
message(sprintf(">>> Canonical status index loaded: %d occupations\n", nrow(sigma)))

# ------------------------------------------------------------------------------
# Prepare risk sets
# Attaches status gap variables and z-scores delta_importance within flow.
# Z-scoring is done within each flow separately (adoption and abandonment
# have different distributions of importance change), so coefficients are
# in SD-outcome units and directionally comparable across specifications.
# ------------------------------------------------------------------------------
prep_riskset <- function(path, flow_label, std_outcome_name, fixed_outcome_name) {
  message(sprintf(">>> Loading %s (%s)...", path, flow_label))
  rs <- readRDS(path)
  setDT(rs)
  message(sprintf("  %d rows, %d columns", nrow(rs), ncol(rs)))

  cols <- names(rs)
  src  <- detect_col("^source$",                    cols, "source")
  tgt  <- detect_col("^target$",                    cols, "target")
  skl  <- detect_col(c("^skill_name$", "^skill$"),  cols, "skill")
  dist <- detect_col(c("structural.?dist"),          cols, "distance")
  arch <- detect_col(c("archetype", "skill.?class"), cols, "archetype")

  required_outcomes <- c(std_outcome_name, fixed_outcome_name, "delta_importance")
  missing_outcomes  <- setdiff(required_outcomes, cols)
  if (length(missing_outcomes) > 0) {
    stop(sprintf(
      "%s is missing required outcome columns: %s\n  Re-run 01a/01b and 02a/02b first.",
      path, paste(missing_outcomes, collapse = ", ")))
  }
  message(sprintf("  [OK] Outcomes present: %s", paste(required_outcomes, collapse = ", ")))

  rs[[src]] <- as.character(rs[[src]])
  rs[[tgt]] <- as.character(rs[[tgt]])

  # Status gap
  rs <- merge(rs, sigma[, .(occ, sigma_source = status_pc1)],
              by.x = src, by.y = "occ", all.x = TRUE)
  rs <- merge(rs, sigma[, .(occ, sigma_target = status_pc1)],
              by.x = tgt, by.y = "occ", all.x = TRUE)
  stopifnot("Unmatched occupations against canonical status index" =
              sum(is.na(rs$sigma_source) | is.na(rs$sigma_target)) == 0)

  rs[, G_ij      := sigma_target - sigma_source]
  rs[, pc1_up    := pmax(0,  G_ij)]
  rs[, pc1_down  := pmin(0,  G_ij)]
  rs[, pc1_dummy := as.integer(G_ij > 0)]

  # Z-score delta_importance within this flow
  mu_di <- mean(rs$delta_importance, na.rm = TRUE)
  sd_di <- sd(rs$delta_importance,   na.rm = TRUE)
  rs[, delta_importance := (delta_importance - mu_di) / sd_di]
  message(sprintf("  [OK] delta_importance z-scored (mean=%.4f, sd=%.4f before scaling)",
                  mu_di, sd_di))

  list(data          = rs,
       src           = src,
       tgt           = tgt,
       skl           = skl,
       dist          = dist,
       arch          = arch,
       flow          = flow_label,
       std_outcome   = std_outcome_name,
       fixed_outcome = fixed_outcome_name)
}

adoption    <- prep_riskset(
  "data/derived/riskset_adoption.rds",
  "adoption",    "diffusion",   "diffusion_fixed_denom")

abandonment <- prep_riskset(
  "data/derived/riskset_abandonment.rds",
  "abandonment", "abandonment", "abandonment_fixed_denom")

# ------------------------------------------------------------------------------
# Fit one model: outcome x panel (A or B) x flow
#
# Panel A: source + skill FE  (rs_obj$src + rs_obj$skl)
# Panel B: target + skill FE  (rs_obj$tgt + rs_obj$skl)
#
# Formula mirrors 05a/05b exactly:
#   outcome ~ (pc1_dummy + pc1_up + pc1_down + dist):archetype | fe1 + fe2
# ------------------------------------------------------------------------------
fit_model <- function(rs_obj, outcome, panel, is_continuous) {

  fe1 <- if (panel == "A") rs_obj$src else rs_obj$tgt
  fe2 <- rs_obj$skl

  fml <- as.formula(sprintf(
    "%s ~ (pc1_dummy + pc1_up + pc1_down + %s):%s | %s + %s",
    outcome, rs_obj$dist, rs_obj$arch, fe1, fe2
  ))
  cluster_fml <- as.formula(
    sprintf("~ %s + %s + %s", rs_obj$src, rs_obj$tgt, rs_obj$skl))

  t0 <- Sys.time()
  m <- if (is_continuous) {
    feols(fml, data = rs_obj$data, cluster = cluster_fml)
  } else {
    feglm(fml, data = rs_obj$data,
          family  = binomial(link = "cloglog"),
          cluster = cluster_fml,
          lean    = TRUE)
  }
  elapsed <- as.numeric(Sys.time() - t0, units = "mins")
  message(sprintf("  Fitted Panel %s | %s | %s in %.1f min",
                  panel, outcome, rs_obj$flow, elapsed))

  ct <- coeftable(m)
  dt <- data.table(
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"])

  # Parse e.g. "pc1_up:archetypeSC_Specialized" -> var + archetype
  rx <- sprintf("^(pc1_dummy|pc1_up|pc1_down|%s):%s(.+)$",
                rs_obj$dist, rs_obj$arch)
  dt <- dt[grepl(rx, term)]
  dt[, `:=`(
    var       = sub(rx, "\\1", term),
    archetype = sub(rx, "\\2", term)
  )]
  dt[, term    := NULL]
  dt[, flow    := rs_obj$flow]
  dt[, panel   := panel]
  dt[, outcome := outcome]
  dt[, n       := nobs(m)]
  dt[]
}

# ------------------------------------------------------------------------------
# Run all 12 models: 2 flows x 2 panels x 3 outcomes
# delta_importance is OLS/continuous; the others are cloglog/binary.
# The 8 binary models are identical to v2 — only the 4 OLS models change
# because delta_importance is now z-scored.
# ------------------------------------------------------------------------------
message("\n>>> Fitting 12 models (2 flows x 2 panels x 3 outcomes)...\n")

results <- rbindlist(list(
  # --- Adoption Panel A ---
  fit_model(adoption, adoption$std_outcome,   panel = "A", is_continuous = FALSE),
  fit_model(adoption, adoption$fixed_outcome, panel = "A", is_continuous = FALSE),
  fit_model(adoption, "delta_importance",     panel = "A", is_continuous = TRUE),
  # --- Adoption Panel B ---
  fit_model(adoption, adoption$std_outcome,   panel = "B", is_continuous = FALSE),
  fit_model(adoption, adoption$fixed_outcome, panel = "B", is_continuous = FALSE),
  fit_model(adoption, "delta_importance",     panel = "B", is_continuous = TRUE),
  # --- Abandonment Panel A ---
  fit_model(abandonment, abandonment$std_outcome,   panel = "A", is_continuous = FALSE),
  fit_model(abandonment, abandonment$fixed_outcome, panel = "A", is_continuous = FALSE),
  fit_model(abandonment, "delta_importance",         panel = "A", is_continuous = TRUE),
  # --- Abandonment Panel B ---
  fit_model(abandonment, abandonment$std_outcome,   panel = "B", is_continuous = FALSE),
  fit_model(abandonment, abandonment$fixed_outcome, panel = "B", is_continuous = FALSE),
  fit_model(abandonment, "delta_importance",         panel = "B", is_continuous = TRUE)
))

# ------------------------------------------------------------------------------
# Sign convention correction for delta_importance under abandonment:
# abandonment = decrease in importance -> negative raw change.
# A coefficient predicting more abandonment predicts a more negative raw
# change, so abandonment raw estimates are multiplied by -1 before
# comparing signs to the standard RCA-based specification.
# ------------------------------------------------------------------------------
results[, estimate_adj := fcase(
  flow == "abandonment" & outcome == "delta_importance", -estimate,
  default = estimate
)]

# Label outcome roles uniformly
results[, outcome_role := fcase(
  outcome %in% c("diffusion",   "abandonment"),                       "standard",
  outcome %in% c("diffusion_fixed_denom", "abandonment_fixed_denom"), "fixed_denom",
  outcome == "delta_importance",                                        "delta_importance"
)]

# ------------------------------------------------------------------------------
# Comparison table: slopes only (pc1_up, pc1_down), wide format
# For delta_importance, use the sign-adjusted estimate.
# ------------------------------------------------------------------------------
message("\n>>> Assembling comparison table...")

slopes <- results[var %in% c("pc1_up", "pc1_down")]

slopes[, estimate_for_table := fifelse(
  outcome_role == "delta_importance", estimate_adj, estimate)]

comparison <- dcast(
  slopes,
  flow + panel + archetype + var ~ outcome_role,
  value.var = c("estimate_for_table", "se")
)

setnames(comparison,
  old = c("estimate_for_table_standard",
          "se_standard",
          "estimate_for_table_fixed_denom",
          "se_fixed_denom",
          "estimate_for_table_delta_importance",
          "se_delta_importance"),
  new = c("est_standard", "se_standard",
          "est_fixed",    "se_fixed",
          "est_raw",      "se_raw")
)

setorder(comparison, flow, panel, archetype, var)

# Sign-match diagnostics
comparison[, sign_match_fixed := sign(est_standard) == sign(est_fixed)]
comparison[, sign_match_raw   := sign(est_standard) == sign(est_raw)]

print(comparison)

# Summary
message("\n  SUMMARY (Panel A):")
pa <- comparison[panel == "A"]
message(sprintf("    Fixed denom sign preserved: %d / %d",
                sum(pa$sign_match_fixed, na.rm = TRUE), nrow(pa)))
message(sprintf("    Raw importance sign preserved: %d / %d",
                sum(pa$sign_match_raw, na.rm = TRUE), nrow(pa)))

message("\n  SUMMARY (Panel B):")
pb <- comparison[panel == "B"]
message(sprintf("    Fixed denom sign preserved: %d / %d",
                sum(pb$sign_match_fixed, na.rm = TRUE), nrow(pb)))
message(sprintf("    Raw importance sign preserved: %d / %d",
                sum(pb$sign_match_raw, na.rm = TRUE), nrow(pb)))

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
dir.create("output/tables/si", showWarnings = FALSE, recursive = TRUE)
fwrite(comparison, "output/tables/si/tab_S3_0_rca_denom_robustness.csv")
fwrite(results,    "output/tables/si/tab_S3_0_rca_denom_robustness_full.csv")

message("\n  Saved:")
message("    output/tables/si/tab_S3_0_rca_denom_robustness.csv")
message("    output/tables/si/tab_S3_0_rca_denom_robustness_full.csv")
message("\n  READ AS:")
message("    est_standard = Table SM3/SM4 specification (RCA threshold crossing, cloglog)")
message("    est_fixed    = RCA with denominator frozen at 2015 (cloglog)")
message("    est_raw      = z-scored O*NET importance change, sign-adjusted for abandonment (OLS)")
message("    sign_match_* = TRUE if standard and alternative agree in direction")

fread("output/tables/si/tab_S3_0_rca_denom_robustness.csv") |> print(digits=4)
