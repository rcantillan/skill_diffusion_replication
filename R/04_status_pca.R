# ==============================================================================
# 04_status_pca.R
#
# Builds the occupational status index via PCA over LEVELS.
#
# APPROACH:
#   Status is an occupation-level attribute, not a dyad-level one.
#   1. Occupation table with (log_wage, log_edu, cog) — 741/741 complete
#      read from occ_covariates_complete.rds
#   2. Standard PCA (center=TRUE, scale.=TRUE) over that table
#   3. PC1 score = occupational status index
#   4. Dyadic gap = status_target - status_source
#
# Input:  data/derived/occ_covariates_complete.rds
#         data/derived/riskset_adoption.rds  (sign check only)
# Output: output/tables/main/pca_status_loadings.csv
#         output/tables/main/pca_status_decision.csv
#         output/tables/main/occ_status_scores.csv
#         output/figures/main/fig_pca_status.pdf
#
# Next: 05a_gravity_adoption.R
#       05b_gravity_abandonment.R
# ==============================================================================

gc()
library(data.table)
library(ggplot2)
library(ggsci)
library(patchwork)

if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

occ_file   <- "data/derived/occ_covariates_complete.rds"
input_file <- if (exists("PATH_ADOPTION")) PATH_ADOPTION else
              "data/derived/riskset_adoption.rds"

stopifnot("occ_covariates_complete.rds not found" = file.exists(occ_file))
stopifnot("riskset_adoption.rds not found"        = file.exists(input_file))

out_tables <- file.path("output", "tables", "main")
out_figs   <- file.path("output", "figures", "main")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 1 — Load complete occupational covariates (741/741)
# ==============================================================================
message(">>> Step 1: Loading occupational covariates...")

occ_levels <- readRDS(occ_file)
setDT(occ_levels)

message(sprintf("  Occupations: %d", nrow(occ_levels)))
message(sprintf("  NAs wage: %d | edu: %d | cog: %d",
                sum(is.na(occ_levels$wage)),
                sum(is.na(occ_levels$edu)),
                sum(is.na(occ_levels$cog))))

stopifnot("NAs in covariates" =
            occ_levels[!is.na(wage) & !is.na(edu) & !is.na(cog), .N] == nrow(occ_levels))
message("  [OK] Complete coverage verified")

message(sprintf("  log_wage: median=%.3f | range [%.3f, %.3f]",
                median(occ_levels$log_wage),
                min(occ_levels$log_wage),
                max(occ_levels$log_wage)))
message(sprintf("  log_edu:  median=%.3f | range [%.3f, %.3f]",
                median(occ_levels$log_edu),
                min(occ_levels$log_edu),
                max(occ_levels$log_edu)))
message(sprintf("  cog:      median=%.3f | range [%.3f, %.3f]",
                median(occ_levels$cog),
                min(occ_levels$cog),
                max(occ_levels$cog)))

# ==============================================================================
# Step 2 — PCA over occupation levels (center=TRUE, scale.=TRUE)
# ==============================================================================
message("\n>>> Step 2: PCA over occupational levels...")
message("  Variables: log_wage, log_edu, cog")
message("  center=TRUE, scale.=TRUE")

pca_mat <- as.matrix(occ_levels[, .(log_wage, log_edu, cog)])
pca_fit <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

message("\n  Loadings:")
print(round(pca_fit$rotation, 4))

pct_var <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 1)
message("\n  Variance explained:")
for (i in seq_along(pct_var))
  message(sprintf("  PC%d: %.1f%%", i, pct_var[i]))

all_pos_loadings <- all(pca_fit$rotation[, 1] > 0)
message(sprintf("\n  PC1 loadings all positive: %s",
                ifelse(all_pos_loadings, "YES [OK]", "NO — flip will be applied")))

# ==============================================================================
# Step 3 — Status score per occupation
# ==============================================================================
message("\n>>> Step 3: Projecting status scores...")

occ_levels[, status_pc1 := predict(pca_fit, pca_mat)[, 1]]

message(sprintf("  Mean:  %.4f [~0 by centering]", mean(occ_levels$status_pc1)))
message(sprintf("  SD:    %.4f", sd(occ_levels$status_pc1)))
message(sprintf("  Range: [%.4f, %.4f]",
                min(occ_levels$status_pc1),
                max(occ_levels$status_pc1)))

message("\n  Top 5 by status:")
print(occ_levels[order(-status_pc1)][1:5,
      .(occ, status_pc1, log_wage, log_edu, cog)])
message("\n  Bottom 5 by status:")
print(occ_levels[order(status_pc1)][1:5,
      .(occ, status_pc1, log_wage, log_edu, cog)])

# ==============================================================================
# Step 4 — Sign verification of the gap
# ==============================================================================
message("\n>>> Step 4: Sign verification...")

mu_w <- mean(occ_levels$log_wage); sd_w <- sd(occ_levels$log_wage)
mu_e <- mean(occ_levels$log_edu);  sd_e <- sd(occ_levels$log_edu)
mu_c <- mean(occ_levels$cog);      sd_c <- sd(occ_levels$cog)

occ_z <- occ_levels[, .(
  occ,
  wz = (log_wage - mu_w) / sd_w,
  ez = (log_edu  - mu_e) / sd_e,
  cz = (cog      - mu_c) / sd_c
)]

dt_chk <- readRDS(input_file)
setDT(dt_chk)
set.seed(42)
if (nrow(dt_chk) > 500000)
  dt_chk <- dt_chk[sample(.N, 500000)]
dt_chk <- dt_chk[, .(source, target)]
gc()

dt_chk <- merge(dt_chk,
                occ_levels[, .(occ, s_pc1 = status_pc1)],
                by.x = "source", by.y = "occ", all.x = TRUE)
dt_chk <- merge(dt_chk,
                occ_levels[, .(occ, t_pc1 = status_pc1)],
                by.x = "target", by.y = "occ", all.x = TRUE)
dt_chk[, gap := t_pc1 - s_pc1]

dt_chk <- merge(dt_chk,
                occ_z[, .(occ, wz_s = wz, ez_s = ez, cz_s = cz)],
                by.x = "source", by.y = "occ", all.x = TRUE)
dt_chk <- merge(dt_chk,
                occ_z[, .(occ, wz_t = wz, ez_t = ez, cz_t = cz)],
                by.x = "target", by.y = "occ", all.x = TRUE)

gap_up   <- dt_chk[wz_t > wz_s & ez_t > ez_s & cz_t > cz_s,
                   mean(gap, na.rm = TRUE)]
gap_down <- dt_chk[wz_t < wz_s & ez_t < ez_s & cz_t < cz_s,
                   mean(gap, na.rm = TRUE)]
gap_peer <- dt_chk[abs(wz_t - wz_s) < 0.1 &
                   abs(ez_t - ez_s) < 0.1 &
                   abs(cz_t - cz_s) < 0.1,
                   mean(gap, na.rm = TRUE)]

message(sprintf("  Gap upward dyads:   %+.4f  [must be POSITIVE]", gap_up))
message(sprintf("  Gap downward dyads: %+.4f  [must be NEGATIVE]", gap_down))
message(sprintf("  Gap peer dyads:     %+.4f  [must be ~0]",       gap_peer))

needs_flip <- (!is.na(gap_up) && gap_up < 0)

if (needs_flip) {
  message("  -> FLIP applied (PC1 was oriented in reverse)")
  occ_levels[, status_pc1 := -status_pc1]
  gap_up   <- -gap_up
  gap_down <- -gap_down
  gap_peer <- -gap_peer
} else {
  message("  -> Sign correct, no flip applied")
}

stopifnot("Upward gap must be positive"  = !is.na(gap_up)   && gap_up   > 0)
stopifnot("Downward gap must be negative" = !is.na(gap_down) && gap_down < 0)
message("  [OK] Sign convention verified")

rm(dt_chk, occ_z); gc()

# ==============================================================================
# Step 5 — Diagnostic figure (English labels)
# ==============================================================================
message("\n>>> Step 5: Diagnostic figure...")

theme_sa <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(size = 11),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 10, colour = "grey10"),
    panel.grid        = element_blank(),
    axis.line         = element_line(linewidth = 0.35, colour = "grey20"),
    axis.ticks        = element_line(linewidth = 0.3),
    axis.ticks.length = unit(2, "pt")
  )

# Panel A: Loadings
load_dt   <- as.data.table(pca_fit$rotation, keep.rownames = "var")
load_long <- melt(load_dt, id.vars = "var",
                  variable.name = "PC", value.name = "loading")
load_long[, var := factor(var,
  levels = c("log_wage", "log_edu", "cog"),
  labels = c("Log wage", "Log education", "Cognitive score"))]
load_long[, pc_label := paste0(as.character(PC), "\n(",
                                pct_var[match(as.character(PC),
                                              paste0("PC", seq_along(pct_var)))],
                                "%)")]

p_load <- ggplot(load_long, aes(x = var, y = loading, fill = var)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_text(aes(label = round(loading, 3),
                vjust = ifelse(loading >= 0, -0.3, 1.2)),
            size = 2.8, family = "Helvetica") +
  facet_wrap(~ pc_label, ncol = 3) +
  scale_fill_aaas() +
  scale_y_continuous(limits = c(-1.05, 1.05),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = NULL, y = "Loading") +
  theme_sa +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

# Panel B: Status score distribution
p_score <- ggplot(occ_levels, aes(x = status_pc1)) +
  geom_histogram(fill = pal_aaas()(1), colour = "white",
                 bins = 35, linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  labs(x = "Status score (PC1)", y = "N occupations") +
  theme_sa

message("\n  Correlations between level variables:")
print(round(cor(occ_levels[, .(log_wage, log_edu, cog)],
                use = "complete.obs"), 3))

fig_out <- (p_load / (p_score | plot_spacer())) +
  plot_annotation(
    subtitle = sprintf(
      "PC1 explains %.1f%% of variance (PC2: %.1f%%). Gap = status_target - status_source. Flip applied: %s.",
      pct_var[1], pct_var[2], ifelse(needs_flip, "YES", "NO")),
    theme = theme(
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

ggsave(
  filename = file.path(out_figs, "fig_pca_status.pdf"),
  plot     = fig_out,
  width    = 12, height = 8, units = "in",
  device   = cairo_pdf, bg = "white"
)
message("  Saved: output/figures/main/fig_pca_status.pdf")

# ==============================================================================
# Step 6 — Save outputs
# ==============================================================================
message("\n>>> Step 6: Saving outputs...")

fwrite(as.data.table(pca_fit$rotation, keep.rownames = "var"),
       file.path(out_tables, "pca_status_loadings.csv"))

fwrite(occ_levels[, .(occ, log_wage, log_edu, cog, status_pc1)],
       file.path(out_tables, "occ_status_scores.csv"))

fwrite(data.table(
  approach         = "occupation_level_pca",
  variables        = "log_wage, log_edu, cog",
  center_scale     = "TRUE / TRUE",
  n_occupations    = nrow(occ_levels),
  pc1_pct_var      = pct_var[1],
  pc2_pct_var      = pct_var[2],
  all_pos_loadings = all_pos_loadings,
  needs_flip       = needs_flip,
  gap_up           = round(gap_up,   4),
  gap_down         = round(gap_down, 4),
  gap_peer         = round(gap_peer, 4)
), file.path(out_tables, "pca_status_decision.csv"))

message("  pca_status_loadings.csv")
message("  occ_status_scores.csv")
message("  pca_status_decision.csv")

message(sprintf("\n  SUMMARY:"))
message(sprintf("  Occupations:       %d", nrow(occ_levels)))
message(sprintf("  PC1 variance:      %.1f%%", pct_var[1]))
message(sprintf("  All loadings +: %s", ifelse(all_pos_loadings, "YES", "NO")))
message(sprintf("  Flip applied:   %s", ifelse(needs_flip, "YES", "NO")))

if (pct_var[1] >= 70 && all_pos_loadings)
  message("  [OK] PCA valid — proceed with 05a_gravity_adoption.R")
if (pct_var[1] < 70)
  message(sprintf("  [!!] PC1 explains only %.1f%% — evaluate adequacy",
                  pct_var[1]))

gc()
message("\n>>> 04_status_pca.R complete.")
message("    Next: 05a_gravity_adoption.R")
message("          05b_gravity_abandonment.R")
