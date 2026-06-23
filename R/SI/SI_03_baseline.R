# ==============================================================================
# SI_03_baseline.R  —  Two-domain baseline model (adoption + abandonment)
#
# Estimates the cloglog gravity model with 2 domains (Cognitive / Physical)
# for adoption and abandonment, both panels (Source FE and Target FE).
# This baseline is the reference point for all SI robustness tests.
#
# Output:
#   output/models/si/baseline_adopt_A.rds
#   output/models/si/baseline_adopt_B.rds
#   output/models/si/baseline_aband_A.rds
#   output/models/si/baseline_aband_B.rds
#   output/tables/si/baseline_coefs_2d.csv
# ==============================================================================

source("R/SI/00_setup_SI.R")

COEF_KEY <- c("b_up_Cog", "b_dn_Cog", "b_up_Phy", "b_dn_Phy")

# ==============================================================================
# Function: estimate one panel
# ==============================================================================
run_panel <- function(setup, fe_vars, panel_label) {
  message(sprintf("\n  Estimating %s %s (%s)...",
                  setup$flow, panel_label, paste(fe_vars, collapse="+")))
  t0 <- proc.time()["elapsed"]
  m  <- feglm(
    setup$fml,
    data      = setup$dt,
    family    = binomial("cloglog"),
    fixef     = fe_vars,
    cluster   = c("source", "target", "skill_name"),
    lean      = TRUE,
    mem.clean = TRUE,
    nthreads  = 0
  )
  elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
  message(sprintf("  Completed in %.1f min", elapsed))
  m
}

all_coefs <- list()

# ==============================================================================
# ADOPTION
# ==============================================================================
message("\n", strrep("=", 60))
message(">>> ADOPTION")
message(strrep("=", 60))

setup_adopt <- load_flow("adoption")

# Panel A — Source FE
m_adopt_A <- run_panel(setup_adopt,
                       c("source", "skill_name"), "Panel A")
coefs_adopt_A <- extract_coefs_2d(m_adopt_A, "Panel A", "adoption")
saveRDS(m_adopt_A, file.path(SI_MODELS, "baseline_adopt_A.rds"))
rm(m_adopt_A); gc()

# Panel B — Target FE
m_adopt_B <- run_panel(setup_adopt,
                       c("target", "skill_name"), "Panel B")
coefs_adopt_B <- extract_coefs_2d(m_adopt_B, "Panel B", "adoption")
saveRDS(m_adopt_B, file.path(SI_MODELS, "baseline_adopt_B.rds"))
rm(m_adopt_B); gc()

all_coefs[["adopt_A"]] <- coefs_adopt_A
all_coefs[["adopt_B"]] <- coefs_adopt_B
rm(setup_adopt); gc()

# ==============================================================================
# ABANDONMENT
# ==============================================================================
message("\n", strrep("=", 60))
message(">>> ABANDONMENT")
message(strrep("=", 60))

setup_aband <- load_flow("abandonment")

# Panel A — Source FE
m_aband_A <- run_panel(setup_aband,
                       c("source", "skill_name"), "Panel A")
coefs_aband_A <- extract_coefs_2d(m_aband_A, "Panel A", "abandonment")
saveRDS(m_aband_A, file.path(SI_MODELS, "baseline_aband_A.rds"))
rm(m_aband_A); gc()

# Panel B — Target FE
m_aband_B <- run_panel(setup_aband,
                       c("target", "skill_name"), "Panel B")
coefs_aband_B <- extract_coefs_2d(m_aband_B, "Panel B", "abandonment")
saveRDS(m_aband_B, file.path(SI_MODELS, "baseline_aband_B.rds"))
rm(m_aband_B); gc()

all_coefs[["aband_A"]] <- coefs_aband_A
all_coefs[["aband_B"]] <- coefs_aband_B
rm(setup_aband); gc()

# ==============================================================================
# Save and report
# ==============================================================================
baseline_all <- rbindlist(all_coefs)
fwrite(baseline_all, file.path(SI_TABLES, "baseline_coefs_2d.csv"))

# ==============================================================================
# Re-extract coefficients from saved models
# (robust to any grepl issues — reads directly from coeftable)
# ==============================================================================
message("\n>>> Re-extracting coefficients from saved models...")

baseline_all <- rbind(
  extract_coefs_2d(readRDS(file.path(SI_MODELS, "baseline_adopt_A.rds")),
                   "Panel A", "adoption"),
  extract_coefs_2d(readRDS(file.path(SI_MODELS, "baseline_adopt_B.rds")),
                   "Panel B", "adoption"),
  extract_coefs_2d(readRDS(file.path(SI_MODELS, "baseline_aband_A.rds")),
                   "Panel A", "abandonment"),
  extract_coefs_2d(readRDS(file.path(SI_MODELS, "baseline_aband_B.rds")),
                   "Panel B", "abandonment")
)

fwrite(baseline_all, file.path(SI_TABLES, "baseline_coefs_2d.csv"))

message("\n>>> BASELINE SUMMARY (ATC coefficients):")
print(baseline_all[coef %in% COEF_KEY,
                   .(flow, panel, coef,
                     estimate  = round(estimate,  3),
                     std_error = round(std_error, 3),
                     sig)])

# Sign verification
message("\n>>> SIGN VERIFICATION:")
for (fl in c("adoption", "abandonment")) {
  b_up_cog <- baseline_all[flow==fl & panel=="Panel A" &
                            coef=="b_up_Cog", estimate]
  b_up_phy <- baseline_all[flow==fl & panel=="Panel A" &
                            coef=="b_up_Phy", estimate]
  if (length(b_up_cog) == 0) b_up_cog <- NA_real_
  if (length(b_up_phy) == 0) b_up_phy <- NA_real_
  atc_gap <- b_up_cog - b_up_phy
  exp_cog <- if (fl=="adoption") "> 0" else "< 0"
  exp_phy <- if (fl=="adoption") "< 0" else "> 0"
  ok_cog  <- if (fl=="adoption") isTRUE(b_up_cog > 0) else isTRUE(b_up_cog < 0)
  ok_phy  <- if (fl=="adoption") isTRUE(b_up_phy < 0) else isTRUE(b_up_phy > 0)
  message(sprintf(
    "  %s Panel A: b_up_Cog=%.3f [%s %s] | b_up_Phy=%.3f [%s %s] | ATC gap=%.3f",
    fl,
    b_up_cog, exp_cog, ifelse(ok_cog, "OK", "!!"),
    b_up_phy, exp_phy, ifelse(ok_phy, "OK", "!!"),
    atc_gap))
}

message("\n>>> SI_04_baseline.R complete.")
message("    Outputs → output/models/si/ and output/tables/si/")
message("    Next: SI_06_placebo.R")