# ==============================================================================
# utils.R
#
# Shared utility functions for the analysis pipeline.
# Source this file at the top of any script that needs extract_coefs().
#
# Usage:
#   source("R/utils.R")
# ==============================================================================

# ------------------------------------------------------------------------------
# extract_coefs()
#
# Extract ATC-relevant coefficients from a fitted fixest model and return
# them as a tidy data.table ready for table assembly and plotting.
#
# The function handles both term orderings that fixest may produce:
#   "wage_up:domainCognitive"  or  "domainCognitive:wage_up"
#
# Parameters:
#   m             fixest model object
#   panel_short   character: "Panel A" or "Panel B"
#   threshold_val numeric (optional): RCA threshold value, used in sensitivity
#                 scripts to tag rows (e.g., 0.90, 1.00, 1.10, 1.25)
#
# Returns:
#   data.table with columns:
#     panel_short, coef, estimate, se [, threshold]
#
# Coefficient labels:
#   Theta_up_{domain}  — upward wage gradient (wage_up)
#   Theta_dn_{domain}  — downward wage gradient (wage_down)
#   kappa_{domain}     — directional entry discontinuity (up_dummy)
#   delta_{domain}     — structural distance friction
# ------------------------------------------------------------------------------
extract_coefs <- function(m, panel_short, threshold_val = NULL) {

  ct <- fixest::coeftable(m)
  rn <- rownames(ct)

  # Try each candidate term name in order; return first match
  find_one <- function(pats) {
    for (p in pats) {
      h <- grep(p, rn, value = TRUE)[1L]
      if (!is.na(h)) return(c(ct[h, "Estimate"], ct[h, "Std. Error"]))
    }
    c(NA_real_, NA_real_)
  }

  params <- list(
    Theta_up_Cog = find_one(c("wage_up:domainCognitive",
                               "domainCognitive:wage_up")),
    Theta_dn_Cog = find_one(c("wage_down:domainCognitive",
                               "domainCognitive:wage_down")),
    kappa_Cog    = find_one(c("up_dummy:domainCognitive",
                               "domainCognitive:up_dummy")),
    delta_Cog    = find_one(c("structural_distance:domainCognitive",
                               "domainCognitive:structural_distance")),
    Theta_up_Phy = find_one(c("wage_up:domainPhysical",
                               "domainPhysical:wage_up")),
    Theta_dn_Phy = find_one(c("wage_down:domainPhysical",
                               "domainPhysical:wage_down")),
    kappa_Phy    = find_one(c("up_dummy:domainPhysical",
                               "domainPhysical:up_dummy")),
    delta_Phy    = find_one(c("structural_distance:domainPhysical",
                               "domainPhysical:structural_distance"))
  )

  out <- data.table::data.table(
    panel_short = panel_short,
    coef        = names(params),
    estimate    = vapply(params, `[`, numeric(1), 1),
    se          = vapply(params, `[`, numeric(1), 2)
  )

  if (!is.null(threshold_val)) out[, threshold := threshold_val]

  out
}

# ------------------------------------------------------------------------------
# extract_coefs_archetype()
#
# Variant of extract_coefs() for models with the three-archetype specification:
#   SC_Scaffolding, SC_Specialized (reference), Physical_Terminal
#
# Parameters:
#   m             fixest model object
#   panel_short   character: "Panel A" or "Panel B"
#
# Returns:
#   data.table with columns:
#     panel_short, archetype, coef, estimate, se
# ------------------------------------------------------------------------------
extract_coefs_archetype <- function(m, panel_short) {

  ct <- fixest::coeftable(m)
  rn <- rownames(ct)

  archetypes <- c("SC_Specialized", "SC_Scaffolding", "Physical_Terminal")

  find_one <- function(pats) {
    for (p in pats) {
      h <- grep(p, rn, fixed = FALSE, value = TRUE)[1L]
      if (!is.na(h)) return(c(ct[h, "Estimate"], ct[h, "Std. Error"]))
    }
    c(NA_real_, NA_real_)
  }

  rows <- lapply(archetypes, function(arch) {
    suf <- paste0("atc_archetype", arch)
    params <- list(
      Theta_up = find_one(c(paste0("wage_up:", suf),
                             paste0(suf, ":wage_up"))),
      Theta_dn = find_one(c(paste0("wage_down:", suf),
                             paste0(suf, ":wage_down"))),
      kappa    = find_one(c(paste0("up_dummy:", suf),
                             paste0(suf, ":up_dummy"))),
      delta    = find_one(c(paste0("structural_distance:", suf),
                             paste0(suf, ":structural_distance")))
    )
    data.table::data.table(
      panel_short = panel_short,
      archetype   = arch,
      coef        = names(params),
      estimate    = vapply(params, `[`, numeric(1), 1),
      se          = vapply(params, `[`, numeric(1), 2)
    )
  })

  data.table::rbindlist(rows)
}