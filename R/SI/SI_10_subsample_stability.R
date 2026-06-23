# ==============================================================================
# SI_10_subsample_stability.R  —  Subsample stability (3 skill types)
#
# Re-estimates the baseline with 3 independent seeds (42, 123, 999) to
# verify that directional gravity coefficients do not depend on the specific
# 50% subsample of source occupations drawn for estimation.
#
# Stability criterion: coefficient of variation (CV) < 5% across seeds,
# or absolute SD < 0.02. Both criteria are reported.
#
# Skill types (3 levels):
#   SC_Scaffolding    — specialized socio-cognitive
#   SC_Specialized    — general socio-cognitive
#   Physical_Terminal — physical-sensory
#
# Runs for ADOPTION and ABANDONMENT, Panel A and Panel B.
#
# Output:
#   output/tables/si/table_S3_subsample_stability.csv
#   output/figures/si/fig_SI_subsample_stability.pdf / .png
# ==============================================================================
 
source("R/SI/00_setup_SI.R")
library(ggplot2)
library(patchwork)
 
SKILL_LEVELS <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")
SKILL_LABELS <- c(
  SC_Scaffolding    = "Specialized socio-cognitive",
  SC_Specialized    = "General socio-cognitive",
  Physical_Terminal = "Physical-sensory"
)
 
COEF_KEY <- c(
  "b_up_SC_Scaffolding",    "b_dn_SC_Scaffolding",
  "b_up_SC_Specialized",    "b_dn_SC_Specialized",
  "b_up_Physical_Terminal", "b_dn_Physical_Terminal"
)
 
SEEDS <- c(42L, 123L, 999L)
 
# ==============================================================================
# extract_coefs_3skill()
# ==============================================================================
extract_coefs_3skill <- function(model, panel_label, flow_label) {
  ct_mat <- coeftable(model)
  ct <- data.table(
    term      = rownames(ct_mat),
    estimate  = ct_mat[, 1L],
    std_error = ct_mat[, 2L],
    t_stat    = ct_mat[, 3L],
    p_value   = ct_mat[, 4L]
  )
  term_map <- list(
    b_up_SC_Scaffolding    = "pc1_up:atc_archetypeSC_Scaffolding",
    b_dn_SC_Scaffolding    = "pc1_down:atc_archetypeSC_Scaffolding",
    b_up_SC_Specialized    = "pc1_up:atc_archetypeSC_Specialized",
    b_dn_SC_Specialized    = "pc1_down:atc_archetypeSC_Specialized",
    b_up_Physical_Terminal = "pc1_up:atc_archetypePhysical_Terminal",
    b_dn_Physical_Terminal = "pc1_down:atc_archetypePhysical_Terminal"
  )
  rows <- lapply(names(term_map), function(coef_name) {
    idx <- which(ct[["term"]] == term_map[[coef_name]])
    if (length(idx) == 0L) {
      warning(sprintf("Term not found: %s", term_map[[coef_name]]))
      return(data.table(coef = coef_name, estimate = NA_real_,
                        std_error = NA_real_, panel = panel_label,
                        flow = flow_label))
    }
    data.table(coef = coef_name, estimate = ct[idx, estimate],
               std_error = ct[idx, std_error],
               panel = panel_label, flow = flow_label)
  })
  rbindlist(rows)
}
 
# ==============================================================================
# attach_archetype() — load_flow() returns domain but not atc_archetype
# ==============================================================================
message(">>> Building atc_archetype lookup...")
archetype_lkp <- unique(
  readRDS("data/derived/riskset_adoption.rds")[,
    .(skill_name, atc_archetype = as.character(atc_archetype))]
)[!is.na(atc_archetype)]
stopifnot(uniqueN(archetype_lkp$skill_name) == nrow(archetype_lkp))
 
attach_archetype <- function(setup) {
  dt <- setup$dt
  dt[archetype_lkp, on = "skill_name", atc_archetype := i.atc_archetype]
  dt[, atc_archetype := factor(atc_archetype, levels = SKILL_LEVELS)]
  dt <- dt[!is.na(atc_archetype)]
  setup$dt <- dt
  setup
}
 
# ==============================================================================
# run_seed() — estimate one flow × panel × seed combination
# ==============================================================================
run_seed <- function(flow, seed, fe_vars, panel_label) {
  message(sprintf("\n  %s %s seed=%d...", flow, panel_label, seed))
  setup <- attach_archetype(load_flow(flow, seed = seed))
 
  outcome <- setup$outcome
  fml <- as.formula(sprintf(
    "%s ~ (up_dummy + pc1_up + pc1_down + structural_distance) : atc_archetype",
    outcome))
 
  t0 <- proc.time()["elapsed"]
  m <- feglm(
    fml,
    data      = setup$dt,
    family    = binomial("cloglog"),
    fixef     = fe_vars,
    cluster   = c("source", "target", "skill_name"),
    lean      = TRUE, mem.clean = TRUE, nthreads = 0
  )
  elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
  message(sprintf("  Done in %.1f min", elapsed))
 
  coefs <- extract_coefs_3skill(m, panel_label, flow)
  coefs[, seed := seed]
  rm(m, setup); gc()
  coefs
}
 
# ==============================================================================
# Main estimation loop
# ==============================================================================
all_seeds <- list()
 
for (flow in c("adoption", "abandonment")) {
  message("\n", strrep("=", 60))
  message(sprintf(">>> %s", toupper(flow)))
  message(strrep("=", 60))
 
  for (seed in SEEDS) {
    key_A <- sprintf("stab_%s_A_seed%d", flow, seed)
    key_B <- sprintf("stab_%s_B_seed%d", flow, seed)
    ckpt_A <- file.path(SI_MODELS, sprintf("%s.rds", key_A))
    ckpt_B <- file.path(SI_MODELS, sprintf("%s.rds", key_B))
 
    if (file.exists(ckpt_A)) {
      message(sprintf("  [cache] %s Panel A seed=%d", flow, seed))
      all_seeds[[key_A]] <- readRDS(ckpt_A)
    } else {
      res_A <- run_seed(flow, seed, c("source", "skill_name"), "Panel A")
      saveRDS(res_A, ckpt_A)
      all_seeds[[key_A]] <- res_A
    }
 
    if (file.exists(ckpt_B)) {
      message(sprintf("  [cache] %s Panel B seed=%d", flow, seed))
      all_seeds[[key_B]] <- readRDS(ckpt_B)
    } else {
      res_B <- run_seed(flow, seed, c("target", "skill_name"), "Panel B")
      saveRDS(res_B, ckpt_B)
      all_seeds[[key_B]] <- res_B
    }
  }
}
 
stab_all <- rbindlist(all_seeds)
fwrite(stab_all, file.path(SI_MODELS, "subsample_all_seeds_3skill.csv"))
 
# ==============================================================================
# Stability statistics
# ==============================================================================
stab_stats <- stab_all[coef %in% COEF_KEY, .(
  est_42  = estimate[seed == 42L],
  est_123 = estimate[seed == 123L],
  est_999 = estimate[seed == 999L],
  mean_est = mean(estimate, na.rm = TRUE),
  sd_est   = sd(estimate,   na.rm = TRUE),
  cv       = sd(estimate, na.rm = TRUE) / abs(mean(estimate, na.rm = TRUE))
), by = .(flow, panel, coef)]
 
stab_stats[, stable := cv < 0.05 | sd_est < 0.02]
 
fwrite(stab_stats, file.path(SI_TABLES, "table_S3_subsample_stability.csv"))
 
message("\n>>> SUBSAMPLE STABILITY:")
print(stab_stats[order(flow, panel, coef),
                 .(flow, panel, coef,
                   est_42  = round(est_42,  3),
                   est_123 = round(est_123, 3),
                   est_999 = round(est_999, 3),
                   cv      = round(cv,      3),
                   stable)])
 
n_stable <- sum(stab_stats$stable, na.rm = TRUE)
n_total  <- nrow(stab_stats)
message(sprintf("\nStable: %d/%d coefficients", n_stable, n_total))
 
# ==============================================================================
# Figure — Subsample stability
# Science Advances style: square panels, AAAS palette, theme_classic
# Layout: 2 rows (direction β↑ / β↓) × 3 columns (skill types)
#         stacked: Adoption / Abandonment
# ==============================================================================
message("\n>>> Generating figure...")
 
if (!exists("stab_all")) {
  stab_all <- fread(file.path(SI_MODELS, "subsample_all_seeds_3skill.csv"))
}
 
plot_dt <- stab_all[coef %in% COEF_KEY & !is.na(estimate)]
 
# Direction label
plot_dt[, direction := fifelse(grepl("^b_up", coef), "\u03b2\u2191", "\u03b2\u2193")]
 
# Skill type label
plot_dt[, skill_type := fcase(
  grepl("SC_Scaffolding",    coef), "Specialized socio-cognitive",
  grepl("SC_Specialized",    coef), "General socio-cognitive",
  grepl("Physical_Terminal", coef), "Physical-sensory"
)]
plot_dt[, skill_type := factor(skill_type,
  levels = c("Specialized socio-cognitive",
             "General socio-cognitive",
             "Physical-sensory"))]
plot_dt[, direction  := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
plot_dt[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]
plot_dt[, seed_label := factor(paste0("seed=", seed),
  levels = c("seed=42", "seed=123", "seed=999"))]
 
PANEL_COLS  <- c("Panel A" = "#2E6EA6", "Panel B" = "#C0392B")
PANEL_FILLS <- c("Panel A" = "#6AAED6", "Panel B" = "#E07070")
PANEL_LBLS  <- c("Panel A" = "Panel A (Source + Skill FE)",
                 "Panel B" = "Panel B (Target + Skill FE)")
 
theme_si <- theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12,
                                     margin = margin(t = 2, b = 4)),
    axis.title        = element_text(size = 13),
    axis.title.x      = element_text(size = 12, margin = margin(t = 5)),
    axis.text         = element_text(size = 11, colour = "grey15"),
    axis.text.x       = element_text(size = 10),
    panel.border      = element_rect(colour = "grey30", fill = NA,
                                     linewidth = 0.7),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.7),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 12),
    legend.key.width  = unit(1.6, "cm"),
    panel.spacing.x   = unit(0.7, "lines"),
    panel.spacing.y   = unit(0.5, "lines"),
    plot.title        = element_text(face = "bold", size = 14,
                                     margin = margin(b = 4))
  )
 
make_flow_plot <- function(flow_name, show_legend = FALSE) {
  pd <- plot_dt[flow_label == flow_name]
 
  p <- ggplot(pd, aes(x = seed_label, y = estimate,
                      colour = panel, fill = panel,
                      group  = panel)) +
    geom_hline(yintercept = 0, colour = "grey60",
               linewidth = 0.35, linetype = "dotted") +
    geom_line(linewidth = 0.85, lineend = "round",
              position = position_dodge(width = 0.2)) +
    geom_point(size = 3.2, shape = 21, stroke = 0.9,
               position = position_dodge(width = 0.2)) +
    # 2 rows (direction: β↑ / β↓) × 3 columns (skill type)
    facet_grid(direction ~ skill_type, scales = "free_y") +
    scale_colour_manual(values = PANEL_COLS, labels = PANEL_LBLS, name = NULL) +
    scale_fill_manual(  values = PANEL_FILLS, labels = PANEL_LBLS, name = NULL) +
    scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(title = flow_name,
         x     = "Random draw (seed)",
         y     = "Estimate (cloglog scale)") +
    theme_si
 
  if (!show_legend)
    p <- p + theme(legend.position = "none")
  p
}
 
p_adopt       <- make_flow_plot("Adoption",    show_legend = FALSE)
p_aband       <- make_flow_plot("Abandonment", show_legend = TRUE)
leg           <- cowplot::get_legend(p_aband)
p_aband_noleg <- p_aband + theme(legend.position = "none")
 
fig_stab <- (p_adopt / p_aband_noleg /
  patchwork::wrap_elements(full = leg)) +
  plot_layout(heights = c(1, 1, 0.06))
 
ggsave(file.path(SI_FIGS, "fig_SI_subsample_stability.pdf"),
       fig_stab, width = 13, height = 11,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(SI_FIGS, "fig_SI_subsample_stability.png"),
       fig_stab, width = 13, height = 11,
       units = "in", dpi = 300, bg = "white")
 
message("  Saved: fig_SI_subsample_stability.pdf / .png")
message("\n>>> SI_11_subsample_stability.R complete.")