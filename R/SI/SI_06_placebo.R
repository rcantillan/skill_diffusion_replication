# ==============================================================================
# SI_06_placebo.R  —  Test (i): Threshold placebo
#
# Displaces the status-gap cutoff by c ∈ {-1, -0.5, -0.25, 0, 0.25, 0.5, 1}
# PC1 units and re-estimates the full gravity model at each displacement.
# c = 0 reproduces the baseline (verified against saved baseline models).
#
# Threat addressed: if the directional asymmetry were an artifact of an
# arbitrary discontinuity at G_ij = 0, equally strong effects would emerge
# at displaced cutoff values. The prediction is that the asymmetry peaks near
# c = 0 and decays as |c| → 1, consistent with progressive dyad
# misclassification.
#
# Skill types (3): SC_General (general socio-cognitive),
#                  SC_Specialized (specialized socio-cognitive),
#                  Physical_Terminal (physical-sensory).
# The model includes atc_archetype as the interaction term instead of the
# coarser domain (Cognitive / Physical) used in earlier versions.
#
# Parameters recovered per panel × flow × cutoff:
#   b_up_SC_General,    b_dn_SC_General
#   b_up_SC_Specialized,    b_dn_SC_Specialized
#   b_up_Physical_Terminal, b_dn_Physical_Terminal
#
# Runs for ADOPTION and ABANDONMENT, both panels (A and B).
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
#         output/tables/si/baseline_coefs_2d.csv  (for c=0 verification)
#
# Output: output/tables/si/test_i_placebo.csv
#         output/figures/si/fig_SI_test_i.pdf / .png
# ==============================================================================

source("R/SI/00_setup_SI.R")
library(ggplot2)
library(patchwork)
library(cowplot)

# Three skill types — ordered for figure facets
SKILL_LEVELS <- c("SC_General", "SC_Specialized", "Physical_Terminal")
SKILL_LABELS <- c("General socio-cognitive",
                  "Specialized socio-cognitive",
                  "Sensory-physical")

# Key coefficients: β↑ and β↓ for each of the three skill types
COEF_KEY <- c(
  "b_up_SC_General",    "b_dn_SC_General",
  "b_up_SC_Specialized",    "b_dn_SC_Specialized",
  "b_up_Physical_Terminal", "b_dn_Physical_Terminal"
)

# Baseline for c=0 verification
baseline <- fread(file.path(SI_TABLES, "baseline_coefs_2d.csv"))

# ==============================================================================
# extract_coefs_3skill()
# Extracts β↑ and β↓ for each of the three skill types from a feglm model.
# Uses exact string matching to avoid grepl() mishandling ':' in term names.
#
# FIX: coeftable() returns a matrix with term names as rownames (not a column).
# Correct approach: grab rownames explicitly before converting to data.table.
# ==============================================================================
extract_coefs_3skill <- function(model, panel_label, flow_label) {
  ct_mat <- coeftable(model)           # matrix; rownames = coefficient names
  ct <- data.table(
    term      = rownames(ct_mat),
    estimate  = ct_mat[, 1L],
    std_error = ct_mat[, 2L],
    t_stat    = ct_mat[, 3L],
    p_value   = ct_mat[, 4L]
  )

  # Diagnostic: print available terms to catch name mismatches early
  message(sprintf("    Available terms (%d): %s",
                  nrow(ct),
                  paste(head(ct$term, 20), collapse = " | ")))

  # Expected term names — fixest interaction format:
  # pc1_up:atc_archetypeSC_General  (continuous × factor level)
  term_map <- list(
    b_up_SC_General    = "pc1_up:atc_archetypeSC_General",
    b_dn_SC_General    = "pc1_down:atc_archetypeSC_General",
    b_up_SC_Specialized    = "pc1_up:atc_archetypeSC_Specialized",
    b_dn_SC_Specialized    = "pc1_down:atc_archetypeSC_Specialized",
    b_up_Physical_Terminal = "pc1_up:atc_archetypePhysical_Terminal",
    b_dn_Physical_Terminal = "pc1_down:atc_archetypePhysical_Terminal"
  )

  rows <- lapply(names(term_map), function(coef_name) {
    pattern <- term_map[[coef_name]]
    idx     <- which(ct[["term"]] == pattern)   # exact match — avoids regex on ':'
    if (length(idx) == 0L) {
      warning(sprintf("Term not found: %s", pattern))
      return(data.table(coef = coef_name, estimate = NA_real_,
                        std_error = NA_real_, panel = panel_label,
                        flow = flow_label))
    }
    data.table(
      coef      = coef_name,
      estimate  = ct[idx, estimate],
      std_error = ct[idx, std_error],
      panel     = panel_label,
      flow      = flow_label
    )
  })
  rbindlist(rows)
}

# ==============================================================================
# load_flow_placebo()
# Loads data keeping s_pc1 and t_pc1 as separate columns so the gap
# can be redefined at each cutoff displacement without reloading.
# Uses atc_archetype (3 levels) instead of domain (2 levels).
# ==============================================================================
load_flow_placebo <- function(flow = c("adoption", "abandonment"),
                              seed = SEED,
                              frac = SAMPLE_FRAC) {
  flow <- match.arg(flow)
  message(sprintf("\n>>> load_flow_placebo('%s')...", flow))

  rds_path    <- if (flow == "adoption")
                   "data/derived/riskset_adoption.rds" else
                   "data/derived/riskset_abandonment.rds"
  scores_path <- "output/tables/main/occ_status_scores.csv"
  outcome_col <- if (flow == "adoption") "diffusion" else "abandonment"

  stopifnot(file.exists(rds_path), file.exists(scores_path))

  dt <- readRDS(rds_path); setDT(dt)

  # Use atc_archetype (3 levels) — replaces domain (2 levels)
  keep <- c("source", "target", "skill_name", outcome_col,
            "atc_archetype", "structural_distance")
  dt <- dt[, ..keep]
  dt[, atc_archetype := factor(as.character(atc_archetype),
                                levels = SKILL_LEVELS)]
  dt <- dt[!is.na(atc_archetype) & !is.na(structural_distance)]

  # Merge status scores — keep s_pc1 and t_pc1 separate for cutoff shift
  scores <- fread(scores_path); setDT(scores)
  scores[, occ := as.character(occ)]
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  dt[scores, on = .(source = occ), s_pc1 := i.status_pc1]
  dt[scores, on = .(target = occ), t_pc1 := i.status_pc1]
  dt <- dt[!is.na(s_pc1) & !is.na(t_pc1)]

  # Status quintile of source (for stratum definition, unchanged)
  qs <- quantile(scores$status_pc1, probs = 0:5/5, na.rm = TRUE)
  scores[, status_q := as.integer(cut(status_pc1, breaks = qs,
                                       labels = 1:5,
                                       include.lowest = TRUE))]
  dt[scores, on = .(source = occ), status_q_source := i.status_q]
  dt <- dt[!is.na(status_q_source)]
  rm(scores); gc()

  # 50% subsample — same sources as baseline for comparability
  ckpt <- file.path(SI_MODELS, sprintf("sources_seed%d.rds", seed))
  if (file.exists(ckpt)) {
    sources_sample <- readRDS(ckpt)
    message(sprintf("  Sources from disk: %d", length(sources_sample)))
  } else {
    set.seed(seed)
    sources_sample <- sample(unique(dt$source),
                             size = round(uniqueN(dt$source) * frac))
    saveRDS(sources_sample, ckpt)
    message(sprintf("  Sources generated: %d of %d",
                    length(sources_sample), uniqueN(dt$source)))
  }
  dt <- dt[source %in% sources_sample]
  message(sprintf("  %s: %s triads | %d src | %d tgt",
                  flow, format(nrow(dt), big.mark = ","),
                  uniqueN(dt$source), uniqueN(dt$target)))

  list(dt = dt, outcome = outcome_col, flow = flow)
}

# ==============================================================================
# run_placebo()
# Re-estimates the model at each cutoff displacement c.
# Formula interacts spline terms with atc_archetype (3 levels).
# ==============================================================================
run_placebo <- function(setup, fe_vars, panel_label) {
  dt         <- copy(setup$dt)
  flow_label <- setup$flow
  outcome    <- setup$outcome
  results    <- list()

  # atc_archetype has 3 levels; fixest will produce one interaction per level.
  # Reference level is SC_General (first factor level); all three slopes
  # are recovered via extract_coefs_3skill() using exact term matching.
  fml <- as.formula(sprintf(
    "%s ~ (up_dummy + pc1_up + pc1_down + structural_distance) : atc_archetype",
    outcome))

  for (c_val in CUTOFF_GRID) {
    message(sprintf("  %s %s | c = %+.2f", flow_label, panel_label, c_val))

    # Redefine gap terms with displaced cutoff
    dt[, pc1_up   := pmax(0, (t_pc1 - s_pc1) - c_val)]
    dt[, pc1_down := pmin(0, (t_pc1 - s_pc1) - c_val)]
    dt[, up_dummy := fifelse((t_pc1 - s_pc1) > c_val, 1L, 0L)]

    m <- feglm(
      fml,
      data      = dt,
      family    = binomial("cloglog"),
      fixef     = fe_vars,
      cluster   = c("source", "target", "skill_name"),
      lean      = TRUE, mem.clean = TRUE, nthreads = 0
    )

    res           <- extract_coefs_3skill(m, panel_label, flow_label)
    res[, c_shift := c_val]
    results[[as.character(c_val)]] <- res
    rm(m); gc()
  }

  rbindlist(results)
}

# ==============================================================================
# Main
# ==============================================================================
all_results <- list()

# Adoption
message("\n", strrep("=", 60))
message(">>> ADOPTION — Test (i)")
message(strrep("=", 60))
setup_a <- load_flow_placebo("adoption")
res_a_A <- run_placebo(setup_a, c("source", "skill_name"), "Panel A")
res_a_B <- run_placebo(setup_a, c("target", "skill_name"), "Panel B")
all_results[["adopt"]] <- rbind(res_a_A, res_a_B)
rm(setup_a, res_a_A, res_a_B); gc()

# Abandonment
message("\n", strrep("=", 60))
message(">>> ABANDONMENT — Test (i)")
message(strrep("=", 60))
setup_b <- load_flow_placebo("abandonment")
res_b_A <- run_placebo(setup_b, c("source", "skill_name"), "Panel A")
res_b_B <- run_placebo(setup_b, c("target", "skill_name"), "Panel B")
all_results[["aband"]] <- rbind(res_b_A, res_b_B)
rm(setup_b, res_b_A, res_b_B); gc()

placebo_all <- rbindlist(all_results)
fwrite(placebo_all, file.path(SI_TABLES, "test_i_placebo.csv"))
message("  Saved: test_i_placebo.csv")

# ==============================================================================
# Verification: c=0 should match baseline within 1%
# ==============================================================================
message("\n>>> Verification c=0 vs baseline:")
for (fl in c("adoption", "abandonment")) {
  for (pa in c("Panel A", "Panel B")) {
    c0  <- placebo_all[flow == fl & panel == pa & c_shift == 0 &
                         coef %in% COEF_KEY,
                       .(coef, est_c0 = estimate)]
    bas <- baseline[flow == fl & panel == pa & coef %in% COEF_KEY,
                    .(coef, est_base = estimate)]
    chk <- merge(c0, bas, by = "coef")
    chk <- chk[!is.na(est_c0) & !is.na(est_base) & est_base != 0]
    if (nrow(chk) == 0) {
      message(sprintf("  %s %s: no valid rows", fl, pa)); next
    }
    chk[, dev := abs(est_c0 - est_base) / abs(est_base)]
    max_dev <- chk[, max(dev, na.rm = TRUE)]
    message(sprintf("  %s %s: max_dev=%.2f%% %s",
                    fl, pa, max_dev * 100,
                    ifelse(max_dev < 0.01, "[OK]", "[CHECK]")))
  }
}

# Report key β↑ coefficients (Panel A)
message("\n>>> Key β↑ coefficients (Panel A):")
print(dcast(
  placebo_all[panel == "Panel A" &
                coef %in% c("b_up_SC_General",
                             "b_up_SC_Specialized",
                             "b_up_Physical_Terminal")],
  flow + c_shift ~ coef, value.var = "estimate"
)[order(flow, c_shift)])

# ==============================================================================
# Figure S5 — Threshold placebo: 6 subpanels per flow (3 types × 2 directions)
# Layout: 2 rows (β↑ top, β↓ bottom) × 3 columns (skill types)
#         stacked vertically: Adoption block / Abandonment block
# ==============================================================================
message("\n>>> Generating Fig. S5...")

# Autonomous reload block
if (!exists("placebo_all")) {
  message("Reloading placebo_all from CSV...")
  placebo_all <- fread(file.path(SI_TABLES, "test_i_placebo.csv"))
}
if (!exists("baseline")) {
  baseline <- fread(file.path(SI_TABLES, "baseline_coefs_2d.csv"))
}
if (!exists("CUTOFF_GRID")) {
  CUTOFF_GRID <- c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1)
}

plot_data <- placebo_all[coef %in% COEF_KEY & !is.na(estimate)]
plot_data[, ci_lo := estimate - 1.96 * std_error]
plot_data[, ci_hi := estimate + 1.96 * std_error]

# Direction label (β↑ / β↓)
plot_data[, direction := fifelse(grepl("^b_up", coef), "\u03b2\u2191", "\u03b2\u2193")]

# Skill type label — extract from coef name
plot_data[, skill_type := fcase(
  grepl("SC_General",    coef), "General socio-cognitive",
  grepl("SC_Specialized",    coef), "Specialized socio-cognitive",
  grepl("Physical_Terminal", coef), "Sensory-physical"
)]
plot_data[, skill_type := factor(skill_type,
  levels = c("General socio-cognitive",
             "Specialized socio-cognitive",
             "Sensory-physical"))]
plot_data[, direction := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
plot_data[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]

# Baseline reference diamonds (Panel A, c=0) — mapped to new skill labels
base_ref <- baseline[panel == "Panel A" & coef %in% COEF_KEY]
base_ref[, direction := fifelse(grepl("^b_up", coef), "\u03b2\u2191", "\u03b2\u2193")]
base_ref[, skill_type := fcase(
  grepl("SC_General",    coef), "General socio-cognitive",
  grepl("SC_Specialized",    coef), "Specialized socio-cognitive",
  grepl("Physical_Terminal", coef), "Sensory-physical"
)]
base_ref[, skill_type := factor(skill_type,
  levels = c("General socio-cognitive",
             "Specialized socio-cognitive",
             "Sensory-physical"))]
base_ref[, direction := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
base_ref[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]

PANEL_COLS  <- c("Panel A" = "#1B4F8A", "Panel B" = "#B03030")
PANEL_FILLS <- c("Panel A" = "#4A90D9", "Panel B" = "#E07070")
PANEL_LBLS  <- c("Panel A" = "Panel A (Source + Skill FE)",
                 "Panel B" = "Panel B (Target + Skill FE)")

theme_si <- theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12,
                                      margin = margin(t = 2, b = 4)),
    axis.title         = element_text(size = 13),
    axis.title.x       = element_text(size = 12, margin = margin(t = 6)),
    axis.text          = element_text(size = 11, colour = "grey15"),
    axis.text.x        = element_text(size = 9),
    panel.border       = element_rect(colour = "grey30", fill = NA,
                                      linewidth = 0.7),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.7),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 12),
    legend.key.width   = unit(1.8, "cm"),
    legend.key.height  = unit(0.4, "cm"),
    panel.spacing.x    = unit(0.7, "lines"),
    panel.spacing.y    = unit(0.5, "lines")
  )

# make_flow_plot(): 2 rows (direction) × 3 columns (skill type)
make_flow_plot <- function(flow_name, show_legend = FALSE) {
  pd <- plot_data[flow_label == flow_name]
  br <- base_ref[flow_label == flow_name]

  p <- ggplot(pd, aes(x = c_shift, y = estimate,
                      colour = panel, fill = panel,
                      linetype = panel)) +
    geom_hline(yintercept = 0, colour = "grey55",
               linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "grey40",
               linewidth = 0.45, linetype = "dotted") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9, lineend = "round",
              position = position_dodge(width = 0.06)) +
    geom_point(size = 3.0, shape = 21, stroke = 0.8,
               position = position_dodge(width = 0.06)) +
    # Diamond: baseline Panel A at c=0
    geom_point(data = br,
               aes(x = 0, y = estimate),
               colour = PANEL_COLS["Panel A"],
               fill   = PANEL_COLS["Panel A"],
               shape = 23, size = 4.0, stroke = 0.8,
               inherit.aes = FALSE) +
    # 3 columns (skill type) × 2 rows (direction), free y-scales
    facet_grid(direction ~ skill_type, scales = "free_y") +
    scale_colour_manual(values = PANEL_COLS, labels = PANEL_LBLS,
                        name = NULL) +
    scale_fill_manual(  values = PANEL_FILLS, labels = PANEL_LBLS,
                        name = NULL) +
    scale_linetype_manual(
      values = c("Panel A" = "solid", "Panel B" = "dashed"),
      labels = PANEL_LBLS, name = NULL) +
    scale_x_continuous(
      breaks = CUTOFF_GRID,
      labels = c("-1.00", "-0.50", "-0.25", "0",
                 "+0.25", "+0.50", "+1.00"),
      expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    labs(title = flow_name,
         x = "Cutoff displacement c (PC1 units)",
         y = "Estimate (cloglog scale)") +
    theme_si +
    theme(plot.title = element_text(face = "bold", size = 14,
                                     margin = margin(b = 4)))

  if (!show_legend)
    p <- p + theme(legend.position = "none")
  p
}

p_adopt       <- make_flow_plot("Adoption",    show_legend = FALSE)
p_aband       <- make_flow_plot("Abandonment", show_legend = TRUE)
leg           <- cowplot::get_legend(p_aband)
p_aband_noleg <- p_aband + theme(legend.position = "none")

fig_S5 <- (p_adopt / p_aband_noleg /
  patchwork::wrap_elements(full = leg)) +
  plot_layout(heights = c(1, 1, 0.06))

ggsave(file.path(SI_FIGS, "fig_SI_test_i.pdf"),
       fig_S5, width = 12, height = 14,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(SI_FIGS, "fig_SI_test_i.png"),
       fig_S5, width = 12, height = 14,
       units = "in", dpi = 300, bg = "white")

message("  Saved: fig_SI_test_i.pdf / .png")
message("\n>>> SI_05_S7_test_i_placebo.R complete.")
message("    Next: SI_07_stratum_permutation.R")

# ------------------------------------------------------------------------------
# Quick check: domain permutation null (unchanged, kept for pipeline continuity)
# ------------------------------------------------------------------------------
if (file.exists("output/tables/si/test_ii_domain_perm.csv")) {
  perm <- fread("output/tables/si/test_ii_domain_perm.csv")
  print(perm[coef %in% c("b_up_SC_General",
                          "b_up_SC_Specialized",
                          "b_up_Physical_Terminal"),
             .(null_mean = round(mean(estimate), 3),
               null_sd   = round(sd(estimate),   3)),
             by = .(flow, panel, coef)][order(flow, panel, coef)])
}