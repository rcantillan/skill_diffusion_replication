# ==============================================================================
# 06_projections.R
#
# Projects estimated models onto the full risk set and compares observed
# vs. projected rates by status quintile.
#
# DESIGN:
#   - status_pc1 read directly from occ_status_scores.csv
#   - Quintiles of the FOCAL occupation (target for both adoption & abandonment)
#   - Panel A = primary specification (Source FE)
#   - Panel B = robustness specification (Target + Skill FE)
#   - Three counterfactuals per panel:
#       DTC model      : full directional model (kappa + b_up + b_dn + delta)
#       Symmetric null : status magnitude only, no directionality (b_avg × |gap| + delta)
#       Distance null  : skill-profile proximity only (delta)
#
# Inputs:
#   data/derived/riskset_adoption.rds
#   data/derived/riskset_abandonment.rds
#   output/tables/main/coefs_pc1_adoption.csv
#   output/tables/abandonment/coefs_pc1_abandonment.csv
#   output/tables/main/occ_status_scores.csv
#
# Outputs:
#   output/figures/main/fig3_unified.pdf   (Panel A + Panel B, shared legend and X-axis)
#   output/figures/main/fig3_panelA.pdf    (backward compatibility with current .tex)
#   output/figures/main/fig3_panelB.pdf    (backward compatibility with current .tex)
#   output/tables/main/proj_quintiles.csv
#   output/tables/main/proj_gradients.csv
#
# Next: 07_fig1_descriptive.R
# ==============================================================================
 
rm(list = ls()); gc(); gc()
library(data.table)
library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)
 
if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")
 
# ==============================================================================
# Paths
# ==============================================================================
adopt_rds   <- "data/derived/riskset_adoption.rds"
aband_rds   <- "data/derived/riskset_abandonment.rds"
coefs_adopt <- "output/tables/main/coefs_pc1_adoption.csv"
coefs_aband <- "output/tables/abandonment/coefs_pc1_abandonment.csv"
scores_path <- "output/tables/main/occ_status_scores.csv"
out_figs    <- "output/figures/main"
out_tables  <- "output/tables/main"
 
for (d in c(out_figs, out_tables))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
 
for (f in c(adopt_rds, aband_rds, coefs_adopt, coefs_aband, scores_path)) {
  if (!file.exists(f)) stop(sprintf(
    "Required file not found: %s\n  Run the pipeline first:\n  Rscript R/01b_risk_set_abandonment.R\n  Rscript R/02b_enrich_abandonment.R\n  Rscript R/03c_nestedness_merge_ab.R", f))
}
message("[OK] All inputs verified")
 
MAIN_ARCHETYPES <- c("SC_General", "SC_Specialized", "Physical_Terminal")
ARCH_LABELS <- c(
  SC_General    = "General socio-cognitive",
  SC_Specialized    = "Specialized socio-cognitive",
  Physical_Terminal = "Sensory-physical"
)
MIN_N <- 30L
 
# ==============================================================================
# Step 1 — Load status scores and build quintiles
# ==============================================================================
message("\n>>> Step 1: Status scores and quintiles...")
 
scores <- fread(scores_path)
setDT(scores)
scores[, occ := as.character(occ)]
 
qs <- quantile(scores$status_pc1, probs = seq(0, 1, 0.2), na.rm = TRUE, type = 7)
qs <- unique(qs)
stopifnot("Not enough unique breaks for 5 quintiles" = length(qs) >= 6)
 
scores[, status_quintile := as.integer(
  cut(status_pc1, breaks = qs, labels = 1:5, include.lowest = TRUE)
)]
message(sprintf("  Occupations: %d | Quintile levels: %d",
                nrow(scores), uniqueN(scores$status_quintile)))
print(scores[, .N, by = status_quintile][order(status_quintile)])
 
# ==============================================================================
# Step 2 — Projection function
# ==============================================================================
calibrate_alpha <- function(H_vec, obs_mean) {
  H_vec <- H_vec[is.finite(H_vec) & H_vec > 0]
  if (!length(H_vec) || !is.finite(obs_mean) ||
      obs_mean <= 0 || obs_mean >= 1) return(1.0)
  f <- function(a) mean(1 - exp(-a * H_vec), na.rm = TRUE) - obs_mean
  if (f(1e-6) * f(100) > 0) return(1.0)
  tryCatch(uniroot(f, c(1e-6, 100), tol = 1e-8)$root, error = function(e) 1.0)
}
 
project_flow <- function(dt, coefs, outcome_col, focal_col,
                         flow_label, panel = "Panel A") {
 
  message(sprintf("\n  Projecting %s (%s)...", flow_label, panel))
 
  .panel <- panel
  ct <- coefs[panel == .panel]
 
  get_coef <- function(arch, var_name) {
    v <- ct[archetype == arch & var == var_name, coef]
    if (!length(v) || all(is.na(v))) return(0)
    v[1]
  }
 
  for (arch in MAIN_ARCHETYPES) {
    b_up  <- get_coef(arch, "b_up")
    b_dn  <- get_coef(arch, "b_dn")
    kappa <- get_coef(arch, "kappa")
    delta <- get_coef(arch, "delta")
    b_avg <- (b_up + b_dn) / 2
 
    dt[as.character(atc_archetype) == arch,
       h_model := exp(kappa * pc1_dummy +
                      b_up  * pc1_up    +
                      b_dn  * pc1_down  +
                      delta * structural_distance)]
 
    dt[as.character(atc_archetype) == arch,
       h_symmetric := exp(b_avg * (pc1_up - pc1_down) +
                          delta * structural_distance)]
 
    dt[as.character(atc_archetype) == arch,
       h_null := exp(delta * structural_distance)]
  }
 
  dt[, focal := get(focal_col)]
  dt[scores, on = .(focal = occ), focal_quintile := i.status_quintile]
 
  agg <- dt[
    !is.na(focal_quintile) & !is.na(h_model) & !is.na(h_null),
    .(h_model     = sum(h_model,     na.rm = TRUE),
      h_symmetric = sum(h_symmetric, na.rm = TRUE),
      h_null      = sum(h_null,      na.rm = TRUE),
      obs_rate    = mean(get(outcome_col), na.rm = TRUE),
      n           = .N),
    by = .(focal, skill_name, focal_quintile, atc_archetype)
  ]
  agg <- agg[n >= MIN_N]
 
  results <- rbindlist(lapply(MAIN_ARCHETYPES, function(arch) {
    sub <- agg[as.character(atc_archetype) == arch]
    if (!nrow(sub)) return(NULL)
    obs_mean    <- mean(sub$obs_rate, na.rm = TRUE)
    a_model     <- calibrate_alpha(sub$h_model,     obs_mean)
    a_symmetric <- calibrate_alpha(sub$h_symmetric, obs_mean)
    a_null      <- calibrate_alpha(sub$h_null,      obs_mean)
    sub[, rate_model     := 1 - exp(-a_model     * h_model)]
    sub[, rate_symmetric := 1 - exp(-a_symmetric * h_symmetric)]
    sub[, rate_null      := 1 - exp(-a_null      * h_null)]
    sub
  }), fill = TRUE)
 
  proj <- results[, .(
    obs       = mean(obs_rate,       na.rm = TRUE),
    model     = mean(rate_model,     na.rm = TRUE),
    symmetric = mean(rate_symmetric, na.rm = TRUE),
    null      = mean(rate_null,      na.rm = TRUE),
    n         = .N
  ), by = .(focal_quintile, atc_archetype)]
 
  proj[, flow  := flow_label]
  proj[, panel := panel]
  setorder(proj, atc_archetype, focal_quintile)
  proj
}
 
# ==============================================================================
# Step 3 — Prepare adoption data
# ==============================================================================
message("\n>>> Step 3: Preparing adoption data...")
 
dt_adopt <- readRDS(adopt_rds)
setDT(dt_adopt)
dt_adopt[, source := as.character(source)]
dt_adopt[, target := as.character(target)]
 
dt_adopt[scores, on = .(source = occ), status_source := i.status_pc1]
dt_adopt[scores, on = .(target = occ), status_target := i.status_pc1]
dt_adopt <- dt_adopt[!is.na(status_source) & !is.na(status_target)]
 
dt_adopt[, status_pc1 := status_target - status_source]
dt_adopt[, pc1_up     := pmax(0, status_pc1)]
dt_adopt[, pc1_down   := pmin(0, status_pc1)]
dt_adopt[, pc1_dummy  := fifelse(status_pc1 > 0, 1L, 0L)]
dt_adopt[, c("status_source", "status_target", "status_pc1") := NULL]
 
if (!"atc_archetype" %in% names(dt_adopt)) {
  cs_med <- dt_adopt[domain == "Cognitive", median(cs, na.rm = TRUE)]
  dt_adopt[, atc_archetype := fcase(
    domain == "Cognitive" & cs >= cs_med, "SC_General",
    domain == "Cognitive" & cs <  cs_med, "SC_Specialized",
    domain == "Physical",                  "Physical_Terminal"
  )]
}
dt_adopt[, atc_archetype := factor(atc_archetype, levels = MAIN_ARCHETYPES)]
dt_adopt <- dt_adopt[atc_archetype %in% MAIN_ARCHETYPES & !is.na(structural_distance)]
 
message(sprintf("  Rows: %s", format(nrow(dt_adopt), big.mark = ",")))
gc()
 
# ==============================================================================
# Step 4 — Project adoption: Panel A and Panel B
# ==============================================================================
message("\n>>> Step 4: Projecting adoption (both panels)...")
 
coefs_a <- fread(coefs_adopt)
 
proj_adopt_A <- project_flow(dt_adopt, coefs_a, "diffusion", "target", "Adoption", "Panel A")
proj_adopt_B <- project_flow(dt_adopt, coefs_a, "diffusion", "target", "Adoption", "Panel B")
 
rm(dt_adopt); gc()
 
# ==============================================================================
# Step 5 — Prepare abandonment data
# ==============================================================================
message("\n>>> Step 5: Preparing abandonment data...")
 
dt_aband <- readRDS(aband_rds)
setDT(dt_aband)
dt_aband[, source := as.character(source)]
dt_aband[, target := as.character(target)]
 
dt_aband[scores, on = .(source = occ), status_source := i.status_pc1]
dt_aband[scores, on = .(target = occ), status_target := i.status_pc1]
dt_aband <- dt_aband[!is.na(status_source) & !is.na(status_target)]
 
dt_aband[, status_pc1 := status_target - status_source]
dt_aband[, pc1_up     := pmax(0, status_pc1)]
dt_aband[, pc1_down   := pmin(0, status_pc1)]
dt_aband[, pc1_dummy  := fifelse(status_pc1 > 0, 1L, 0L)]
dt_aband[, c("status_source", "status_target", "status_pc1") := NULL]
 
if (!"atc_archetype" %in% names(dt_aband)) {
  cs_med <- dt_aband[domain == "Cognitive", median(cs, na.rm = TRUE)]
  dt_aband[, atc_archetype := fcase(
    domain == "Cognitive" & cs >= cs_med, "SC_General",
    domain == "Cognitive" & cs <  cs_med, "SC_Specialized",
    domain == "Physical",                  "Physical_Terminal"
  )]
}
dt_aband[, atc_archetype := factor(atc_archetype, levels = MAIN_ARCHETYPES)]
dt_aband <- dt_aband[atc_archetype %in% MAIN_ARCHETYPES & !is.na(structural_distance)]
 
message(sprintf("  Rows: %s", format(nrow(dt_aband), big.mark = ",")))
gc()
 
# ==============================================================================
# Step 6 — Project abandonment: Panel A and Panel B
# ==============================================================================
message("\n>>> Step 6: Projecting abandonment (both panels)...")
 
coefs_b <- fread(coefs_aband)
 
proj_aband_A <- project_flow(dt_aband, coefs_b, "abandonment", "target", "Abandonment", "Panel A")
proj_aband_B <- project_flow(dt_aband, coefs_b, "abandonment", "target", "Abandonment", "Panel B")
 
rm(dt_aband); gc()
 
# ==============================================================================
# Step 7 — Combine and compute gradients
# ==============================================================================
message("\n>>> Step 7: Combining results and computing Q5-Q1 gradients...")
 
proj_A   <- rbind(proj_adopt_A, proj_aband_A)
proj_B   <- rbind(proj_adopt_B, proj_aband_B)
proj_all <- rbind(proj_A, proj_B)
proj_all[, atc_archetype := factor(atc_archetype, levels = MAIN_ARCHETYPES)]
 
grad <- proj_all[, {
  .(grad_obs       = obs      [focal_quintile == 5] - obs      [focal_quintile == 1],
    grad_model     = model    [focal_quintile == 5] - model    [focal_quintile == 1],
    grad_symmetric = symmetric[focal_quintile == 5] - symmetric[focal_quintile == 1],
    grad_null      = null     [focal_quintile == 5] - null     [focal_quintile == 1])
}, by = .(flow, atc_archetype, panel)]
 
grad[, recovery_model     := round(grad_model     / grad_obs, 3)]
grad[, recovery_symmetric := round(grad_symmetric / grad_obs, 3)]
grad[, advantage_vs_sym   := round(grad_model     / grad_symmetric, 2)]
grad[, advantage_vs_null  := round(grad_model     / grad_null, 2)]
 
message("\n  Q5-Q1 gradients:")
print(grad)
 
fwrite(proj_all, file.path(out_tables, "proj_quintiles.csv"))
fwrite(grad,     file.path(out_tables, "proj_gradients.csv"))
message("  Saved: proj_quintiles.csv | proj_gradients.csv")
 
# ==============================================================================
# Step 8 — Theme y helper de figura
# ==============================================================================
message("\n>>> Step 8: Building unified figure...")
 
Q_LABS  <- c("Q1\n(low)", "Q2", "Q3", "Q4", "Q5\n(high)")
X_LABEL <- "Occupational status quintile (PC1, 75.0% var.)"
 
SERIES_COLORS <- c(
  "Observed"           = "grey20",
  "Model"              = "#00B4D8",
  "Symmetric null"     = "#F4A261",
  "Distance-only null" = "grey60"
)
SERIES_LT <- c(
  "Observed"           = "solid",
  "Model"              = "solid",
  "Symmetric null"     = "solid",
  "Distance-only null" = "dotdash"
)
SERIES_SH <- c(
  "Observed"           = 16L,
  "Model"              = 21L,
  "Symmetric null"     = 22L,
  "Distance-only null" = 4L
)
 
# Base theme: no X-axis title or legend (added globally)
theme_sa_base <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, face = "bold"),
    axis.title.y     = element_text(size = 10),
    axis.title.x     = element_blank(),
    axis.text        = element_text(size = 9, colour = "grey10"),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    axis.line        = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.8, "lines"),
    panel.spacing    = unit(0.7, "lines"),
    plot.margin      = margin(4, 8, 2, 8),
    plot.tag         = element_text(face = "bold", size = 11)
  )
 
# Build subplot (no legend, no X-axis title)
make_subplot <- function(proj_dt, panel_title, fe_subtitle) {
  long <- melt(proj_dt,
               id.vars      = c("focal_quintile", "atc_archetype", "flow", "panel", "n"),
               measure.vars = c("obs", "model", "symmetric", "null"),
               variable.name = "series", value.name = "rate")
  long[, series := fcase(
    series == "obs",       "Observed",
    series == "model",     "Model",
    series == "symmetric", "Symmetric null",
    series == "null",      "Distance-only null"
  )]
  long[, series := factor(series,
    levels = c("Observed", "Model", "Symmetric null", "Distance-only null"))]
  long[, arch_f := factor(as.character(atc_archetype), MAIN_ARCHETYPES, ARCH_LABELS)]
  long[, flow_f := factor(flow, c("Adoption", "Abandonment"))]
  long <- long[!is.na(rate) & is.finite(rate)]

  ggplot(long,
    aes(x = focal_quintile, y = rate,
        colour = series, shape = series, linetype = series, group = series)) +
    geom_line(linewidth = 0.9, lineend = "round") +
    geom_point(size = 2.2, stroke = 0.5, aes(fill = series), colour = "white") +
    facet_grid(rows = vars(flow_f), cols = vars(arch_f), scales = "free_y") +
    scale_colour_manual(values = SERIES_COLORS) +
    scale_fill_manual(  values = SERIES_COLORS) +
    scale_linetype_manual(values = SERIES_LT) +
    scale_shape_manual(   values = SERIES_SH) +
    scale_x_continuous(breaks = 1:5, labels = Q_LABS,
                       expand = expansion(mult = c(.04, .04))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(.02, .08))) +
    labs(x = NULL, y = "Skill rate",
         title    = panel_title,
         subtitle = fe_subtitle) +
    theme_sa_base +
    theme(
      plot.title    = element_text(face = "bold", size = 11,
                                   hjust = 0, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 9.5, colour = "grey30",
                                   hjust = 0, margin = margin(b = 4))
    ) +
    guides(
      colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.0)),
      linetype = guide_legend(nrow = 1),
      shape    = guide_legend(nrow = 1),
      fill     = "none"
    )
}

# ==============================================================================
# Step 9 — Ensamblar figura unificada y guardar
# ==============================================================================
message("\n>>> Step 9: Assembling and saving figures...")

p_A <- make_subplot(proj_A, "Panel A", "Source + Skill FE")
p_B <- make_subplot(proj_B, "Panel B", "Target + Skill FE")

# Extract shared legend from Panel B
p_B_leg    <- p_B + theme(legend.position = "bottom")
leg_shared <- cowplot::get_legend(p_B_leg)

# Shared X-axis label
xlab_grob <- grid::textGrob(
  X_LABEL,
  gp = grid::gpar(fontsize = 10, fontfamily = "Helvetica")
)

# Figura unificada
fig_unified <- (p_A / plot_spacer() / p_B) /
  patchwork::wrap_elements(full = xlab_grob,  clip = FALSE) /
  patchwork::wrap_elements(full = leg_shared, clip = FALSE) +
  patchwork::plot_layout(heights = c(1, 0.04, 1, 0.05, 0.10))

ggsave(file.path(out_figs, "fig3_unified.pdf"),
       fig_unified,
       width = 9.0, height = 11.5, units = "in",
       device = cairo_pdf, bg = "white")
ggsave(file.path(out_figs, "fig3_unified.png"),
       fig_unified,
       width = 9.0, height = 11.5, units = "in",
       dpi = 300, bg = "white")
message("  Saved: fig3_unified.pdf / .png")

# Individual panel figures (backward compatibility)
make_fig_standalone <- function(proj_dt, panel_title, fe_subtitle) {
  make_subplot(proj_dt, panel_title, fe_subtitle) +
    theme(legend.position = "bottom",
          axis.title.x    = element_text(size = 10, margin = margin(t = 10))) +
    labs(x = X_LABEL) +
    guides(
      colour   = guide_legend(nrow = 1, override.aes = list(linewidth = 1.0)),
      linetype = guide_legend(nrow = 1),
      shape    = guide_legend(nrow = 1),
      fill     = "none"
    )
}

fig_A <- make_fig_standalone(proj_A, "Panel A", "Source + Skill FE")
fig_B <- make_fig_standalone(proj_B, "Panel B", "Target + Skill FE")

ggsave(file.path(out_figs, "fig3_panelA.pdf"),
       fig_A, width = 9, height = 6.5, units = "in",
       device = cairo_pdf, bg = "white")
ggsave(file.path(out_figs, "fig3_panelA.png"),
       fig_A, width = 9, height = 6.5, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(out_figs, "fig3_panelB.pdf"),
       fig_B, width = 9, height = 6.5, units = "in",
       device = cairo_pdf, bg = "white")
ggsave(file.path(out_figs, "fig3_panelB.png"),
       fig_B, width = 9, height = 6.5, units = "in",
       dpi = 300, bg = "white")
message("  Saved: fig3_panelA.pdf / .png")
message("  Saved: fig3_panelB.pdf / .png")