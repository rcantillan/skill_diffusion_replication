# ==============================================================================
# 00_run_all.R
#
# Master pipeline script — runs all main analyses in order.
# Execute from the repository root:
#
#   Rscript R/00_run_all.R
#
# Prerequisites:
#   - R/99_paths_local.R must exist (local paths, not tracked)
#   - Raw O*NET and BLS data in data/raw/ (see README for details)
#   - RAM >= 16 GB recommended; configure swap >= 16 GB if limited
#
# Steps 01-02 are computationally heavy (~hours).
# Step 03 (nestedness) takes 2-3 hours; supports checkpoint/resume.
# Steps 05-08 take minutes each once risk sets are built.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("SKILL DIFFUSION — MAIN PIPELINE\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("=", 70), "\n\n")

run_step <- function(script, label) {
  cat(sprintf("\n%s\n>>> %s\n%s\n", strrep("-", 60), label, strrep("-", 60)))
  t0 <- proc.time()["elapsed"]
  tryCatch(
    source(script, local = new.env()),
    error = function(e) {
      cat(sprintf("\n[FAILED] %s\n  Error: %s\n", label, conditionMessage(e)))
      stop(sprintf("Pipeline halted at: %s", label), call. = FALSE)
    }
  )
  elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
  cat(sprintf("\n[DONE] %s  (%.1f min)\n", label, elapsed))
  invisible(elapsed)
}

timings <- list()

# ==============================================================================
# BLOCK 1 — Risk set construction (requires raw O*NET + BLS data)
# ==============================================================================
timings[["01a"]] <- run_step("R/01a_risk_set_adoption.R",    "01a  Risk set — adoption")
timings[["01b"]] <- run_step("R/01b_risk_set_abandonment.R", "01b  Risk set — abandonment")

# ==============================================================================
# BLOCK 2 — Enrichment (nestedness contributions, skill archetypes)
# ==============================================================================
timings[["02a"]] <- run_step("R/02a_enrich_adoption.R",    "02a  Enrich — adoption")
timings[["02b"]] <- run_step("R/02b_enrich_abandonment.R", "02b  Enrich — abandonment")

# ==============================================================================
# BLOCK 3 — Nestedness scores (~2–3 hours; checkpoint/resume supported)
# ==============================================================================
timings[["03"]]  <- run_step("R/03_nestedness.R",        "03   Nestedness NODF scores")
timings[["03b"]] <- run_step("R/03b_nestedness_merge.R", "03b  Merge nestedness — adoption")
timings[["03c"]] <- run_step("R/03c_nestedness_merge_ab.R", "03c  Merge nestedness — abandonment")

# ==============================================================================
# BLOCK 4 — Status index
# ==============================================================================
timings[["04"]] <- run_step("R/04_status_pca.R", "04   Occupational status PCA")

# ==============================================================================
# BLOCK 5 — Gravity models
# ==============================================================================
timings[["05a"]] <- run_step("R/05a_gravity_adoption.R",    "05a  Gravity model — adoption")
timings[["05b"]] <- run_step("R/05b_gravity_abandonment.R", "05b  Gravity model — abandonment")

# ==============================================================================
# BLOCK 6 — Projections and figures
# ==============================================================================
timings[["06"]] <- run_step("R/06_projections.R",      "06   Projections + Fig. 3")
timings[["07"]] <- run_step("R/07_fig1_descriptive.R", "07   Fig. 1 — descriptive gradients")
timings[["08"]] <- run_step("R/08_fig2_regression.R",  "08   Fig. 2 — regression coefficients")

# ==============================================================================
# BLOCK 7 — Export (figures → Overleaf; manuscript tables)
# ==============================================================================
timings[["10"]] <- run_step("R/10_export_manuscript_figures.R", "10   Export figures → Overleaf")
timings[["11"]] <- run_step("R/11_make_manuscript_tables.R",    "11   Manuscript tables")

# ==============================================================================
# Summary
# ==============================================================================
cat(sprintf("\n%s\nPIPELINE COMPLETE\n%s\n", strrep("=", 70), strrep("=", 70)))
total <- round(sum(unlist(timings)) / 60, 1)
cat(sprintf("Total elapsed: %.1f hours\n\n", total))
cat("Step timings (min):\n")
for (nm in names(timings))
  cat(sprintf("  %-6s  %.1f min\n", nm, timings[[nm]]))
cat("\nOutputs:\n")
cat("  output/figures/main/   — Main figures (PDF + PNG)\n")
cat("  output/tables/main/    — Model coefficients, status scores\n")
cat("  output/models/         — Estimated feglm models (.rds)\n")
cat("\nNext: run SI scripts in R/SI/ — see README for order.\n")
