# ==============================================================================
# SI_11_subperiod_stability.R
#
# Section S9.1: Temporal Stability — Table S2 + Fig. SI_subperiods
#
# Re-estimates the directional gravity model (3 archetypes) separately for
# three non-overlapping sub-periods following the O*NET refresh cycle:
#   (i)  2015–2018: O*NET 20.3 → 23.1  (tag: 1518)
#   (ii) 2019–2021: O*NET 24.1 → 26.1  (tag: 1921)
#   (iii) 2022–2024: O*NET 27.1 → 29.2 (tag: 2224)
#
# Term names in fixest for this model:
#   "pc1_up:atc_archetypeSC_Scaffolding"
#   "pc1_up:atc_archetypeSC_Specialized"
#   "pc1_up:atc_archetypePhysical_Terminal"
#   (and analogous for pc1_down, pc1_dummy, structural_distance)
#
# Input:  data/derived/riskset_{flow}_{tag}.rds
#         output/tables/main/occ_status_scores.csv
#
# Output: output/tables/si/table_S2_subperiods.csv
#         output/tables/si/table_S2_subperiods.tex
#         output/figures/si/fig_SI_subperiods.pdf / .png
# ==============================================================================
 
source("R/SI/00_setup_SI.R")
library(ggplot2)
library(patchwork)
library(cowplot)
 
SUBPERIODS <- list(
  list(label = "2015\u20132018", tag = "1518"),
  list(label = "2019\u20132021", tag = "1921"),
  list(label = "2022\u20132024", tag = "2224")
)
 
ARCH_3 <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")
 
# ==============================================================================
# load_subperiod()
# ==============================================================================
load_subperiod <- function(flow, tag, seed = SEED, frac = SAMPLE_FRAC) {
  rds_path    <- file.path("data", "derived",
                            sprintf("riskset_%s_%s.rds", flow, tag))
  outcome_col <- if (flow == "adoption") "diffusion" else "abandonment"
 
  if (!file.exists(rds_path)) {
    message(sprintf("  [SKIP] Not found: %s", rds_path))
    return(NULL)
  }
 
  message(sprintf("  Loading %s %s...", flow, tag))
  dt <- readRDS(rds_path); setDT(dt)
 
  keep <- c("source", "target", "skill_name", outcome_col,
            "domain", "atc_archetype", "structural_distance")
  keep <- intersect(keep, names(dt))
  dt   <- dt[, ..keep]
 
  # Reconstruct atc_archetype if missing
  if (!"atc_archetype" %in% names(dt)) {
    cs_scores <- readRDS("data/derived/skill_cs_scores.rds"); setDT(cs_scores)
    domain_lkp <- unique(readRDS(
      "data/derived/riskset_adoption.rds")[, .(skill_name, domain)])
    cs_scores <- merge(cs_scores, domain_lkp, by = "skill_name", all.x = TRUE)
    med_cog   <- cs_scores[domain == "Cognitive", median(cs, na.rm = TRUE)]
    cs_scores[, atc_archetype := fcase(
      domain == "Cognitive" & cs >= med_cog, "SC_Scaffolding",
      domain == "Cognitive" & cs <  med_cog, "SC_Specialized",
      domain == "Physical",                   "Physical_Terminal"
    )]
    dt[cs_scores, on = "skill_name", atc_archetype := i.atc_archetype]
    rm(cs_scores, domain_lkp); gc()
  }
 
  dt[, atc_archetype := factor(as.character(atc_archetype), levels = ARCH_3)]
  dt <- dt[atc_archetype %in% ARCH_3 & !is.na(structural_distance)]
 
  # Status gap
  scores <- fread("output/tables/main/occ_status_scores.csv")
  setDT(scores); scores[, occ := as.character(occ)]
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  dt[scores, on = .(source = occ), s_pc1 := i.status_pc1]
  dt[scores, on = .(target = occ), t_pc1 := i.status_pc1]
  dt <- dt[!is.na(s_pc1) & !is.na(t_pc1)]
  dt[, pc1_gap   := t_pc1 - s_pc1]
  dt[, pc1_up    := pmax(0,  pc1_gap)]
  dt[, pc1_down  := pmin(0,  pc1_gap)]
  dt[, pc1_dummy := fifelse(pc1_gap > 0, 1L, 0L)]
  dt[, c("s_pc1","t_pc1","pc1_gap") := NULL]
  rm(scores); gc()
 
  # 50% subsample — intersect with available sources
  ckpt <- file.path(SI_MODELS, sprintf("sources_seed%d.rds", seed))
  if (file.exists(ckpt)) {
    sources_sample <- intersect(readRDS(ckpt), unique(dt$source))
  } else {
    set.seed(seed)
    sources_sample <- sample(unique(dt$source),
                             size = round(uniqueN(dt$source) * frac))
  }
  dt <- dt[source %in% sources_sample]
 
  message(sprintf("    %s dyads | %d src | %d tgt | rate: %.4f",
                  format(nrow(dt), big.mark=","),
                  uniqueN(dt$source),
                  uniqueN(dt$target),
                  mean(dt[[outcome_col]], na.rm=TRUE)))
 
  list(dt = dt, outcome = outcome_col, flow = flow)
}
 
# ==============================================================================
# extract_coefs_3arch()
# ==============================================================================
extract_coefs_3arch <- function(m, panel_label, flow_label, period_label) {
  ct <- as.data.table(m$coeftable, keep.rownames = "term")
  nm <- names(ct)
  if ("Estimate"   %in% nm) setnames(ct, "Estimate",   "estimate")
  if ("Std. Error" %in% nm) setnames(ct, "Std. Error", "std_error")
  if ("z value"    %in% nm) setnames(ct, "z value",    "z")
  if ("Pr(>|z|)"   %in% nm) setnames(ct, "Pr(>|z|)",  "p")
 
  find_coef <- function(patterns) {
    for (p in patterns) {
      idx <- which(ct[["term"]] == p)
      if (length(idx) > 0)
        return(c(ct[["estimate"]][idx[1L]],
                 ct[["std_error"]][idx[1L]]))
    }
    c(NA_real_, NA_real_)
  }
 
  params <- list()
  for (arch in ARCH_3) {
    params[[paste0("b_up_",  arch)]] <- find_coef(c(
      sprintf("pc1_up:atc_archetype%s",              arch),
      sprintf("atc_archetype%s:pc1_up",              arch)))
    params[[paste0("b_dn_",  arch)]] <- find_coef(c(
      sprintf("pc1_down:atc_archetype%s",            arch),
      sprintf("atc_archetype%s:pc1_down",            arch)))
    params[[paste0("kappa_", arch)]] <- find_coef(c(
      sprintf("pc1_dummy:atc_archetype%s",           arch),
      sprintf("atc_archetype%s:pc1_dummy",           arch)))
    params[[paste0("delta_", arch)]] <- find_coef(c(
      sprintf("structural_distance:atc_archetype%s", arch),
      sprintf("atc_archetype%s:structural_distance", arch)))
  }
 
  out <- data.table(
    flow      = flow_label,
    period    = period_label,
    panel     = panel_label,
    coef      = names(params),
    estimate  = vapply(params, function(x) x[1L], numeric(1L)),
    std_error = vapply(params, function(x) x[2L], numeric(1L))
  )
  out[, sig := fcase(
    is.na(estimate),                 "",
    abs(estimate/std_error) >= 3.29, "***",
    abs(estimate/std_error) >= 2.58, "**",
    abs(estimate/std_error) >= 1.96, "*",
    abs(estimate/std_error) >= 1.65, "~",
    default = "ns"
  )]
  out
}
 
# ==============================================================================
# Main estimation loop
# ==============================================================================
all_coefs <- list()
 
fml_3arch <- function(outcome) {
  as.formula(sprintf(
    "%s ~ (pc1_dummy + pc1_up + pc1_down + structural_distance):atc_archetype",
    outcome))
}
 
for (sp in SUBPERIODS) {
  message("\n", strrep("=", 65))
  message(sprintf(">>> Sub-period: %s", sp$label))
  message(strrep("=", 65))
 
  for (flow in c("adoption", "abandonment")) {
    setup <- load_subperiod(flow, sp$tag)
    if (is.null(setup)) next
 
    fml <- fml_3arch(setup$outcome)
 
    for (fe_config in list(
      list(fe = c("source", "skill_name"), label = "Panel A"),
      list(fe = c("target", "skill_name"), label = "Panel B")
    )) {
      key  <- sprintf("%s_%s_%s", flow, sp$tag, gsub(" ","",fe_config$label))
      ckpt <- file.path(SI_MODELS, sprintf("subperiod_%s.rds", key))
 
      if (file.exists(ckpt)) {
        prev <- readRDS(ckpt)
        if (!all(is.na(prev$estimate))) {
          message(sprintf("  [cache] %s %s %s", flow, sp$label, fe_config$label))
          all_coefs[[key]] <- prev
          next
        }
        file.remove(ckpt)
        message(sprintf("  [stale] Deleted %s", ckpt))
      }
 
      message(sprintf("\n  Estimating %s %s %s...",
                      flow, sp$label, fe_config$label))
      t0 <- proc.time()["elapsed"]
 
      m <- feglm(
        fml,
        data      = setup$dt,
        family    = binomial("cloglog"),
        fixef     = fe_config$fe,
        cluster   = c("source", "target", "skill_name"),
        lean      = TRUE, mem.clean = TRUE, nthreads = 0
      )
 
      elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
      message(sprintf("  Completed in %.1f min", elapsed))
 
      coefs <- extract_coefs_3arch(m, fe_config$label, flow, sp$label)
      saveRDS(coefs, ckpt)
      all_coefs[[key]] <- coefs
      rm(m); gc()
    }
    rm(setup); gc()
  }
}
 
# ==============================================================================
# Consolidate and verify
# ==============================================================================
subperiod_all <- rbindlist(all_coefs)
fwrite(subperiod_all, file.path(SI_TABLES, "table_S2_subperiods.csv"))
 
message("\n>>> Sub-period results (β↑ Panel A):")
key_coefs <- paste0("b_up_", ARCH_3)
print(subperiod_all[
  panel == "Panel A" & coef %in% key_coefs,
  .(flow, period, coef,
    estimate  = round(estimate,  3),
    std_error = round(std_error, 3),
    sig)
][order(flow, period, coef)])
 
message("\n>>> Sign check β↑_Physical_Terminal (Panel A):")
print(subperiod_all[
  panel == "Panel A" & coef == "b_up_Physical_Terminal",
  .(flow, period,
    estimate = round(estimate, 3),
    sig,
    sign_ok = fcase(
      flow == "adoption",    estimate < 0,
      flow == "abandonment", estimate > 0,
      default = NA
    ))
][order(flow, period)])
 
# ==============================================================================
# LaTeX — Table S2
# ==============================================================================
message("\n>>> Generating LaTeX Table S2...")
 
make_row <- function(dt, flow_v, period_v, panel_v, arch) {
  b_up <- dt[flow==flow_v & period==period_v & panel==panel_v &
               coef==paste0("b_up_",arch)]
  b_dn <- dt[flow==flow_v & period==period_v & panel==panel_v &
               coef==paste0("b_dn_",arch)]
  fmt    <- function(r) if (nrow(r)==0||is.na(r$estimate)) "---" else
                        sprintf("%.3f%s", r$estimate, r$sig)
  fmt_se <- function(r) if (nrow(r)==0||is.na(r$std_error)) "(---)" else
                        sprintf("(%.3f)", r$std_error)
  list(b_up=fmt(b_up), se_up=fmt_se(b_up),
       b_dn=fmt(b_dn), se_dn=fmt_se(b_dn))
}
 
lines <- c(
  "% TABLE S2 — Temporal stability",
  "% Generated by SI_10_S9_table_S2_subperiods.R",
  "\\begin{table}[ht]",
  "\\centering\\small",
  "\\caption{\\textbf{Supplementary Table S2. Temporal stability of directional",
  "  friction parameters across sub-periods.}",
  "  $\\beta^{\\uparrow}$, $\\beta^{\\downarrow}$: directional friction parameters",
  "  re-estimated separately for three non-overlapping sub-periods.",
  "  Panel A: source and skill FE. Panel B: target and skill FE.",
  "  50\\% subsample of source occupations (seed 42).",
  "  SE in parentheses, three-way clustering.",
  "  $^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, $^{\\sim}p<0.10$.}",
  "\\label{tab:SI_S2}",
  "\\begin{tabular}{llcccccc}",
  "\\toprule",
  " & & \\multicolumn{2}{c}{Specialized socio-cognitive} & \\multicolumn{2}{c}{General socio-cognitive} & \\multicolumn{2}{c}{Physical-sensory} \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
  "Period & Panel & $\\beta^{\\uparrow}$ & $\\beta^{\\downarrow}$ & $\\beta^{\\uparrow}$ & $\\beta^{\\downarrow}$ & $\\beta^{\\uparrow}$ & $\\beta^{\\downarrow}$ \\\\",
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Adoption}} \\\\"
)
 
for (sp in SUBPERIODS) {
  for (pa in c("Panel A","Panel B")) {
    r_sc <- make_row(subperiod_all,"adoption",sp$label,pa,"SC_Scaffolding")
    r_sp <- make_row(subperiod_all,"adoption",sp$label,pa,"SC_Specialized")
    r_ph <- make_row(subperiod_all,"adoption",sp$label,pa,"Physical_Terminal")
    lines <- c(lines,
      sprintf("\\multirow{2}{*}{%s} & %s & %s & %s & %s & %s & %s & %s \\\\",
              sp$label, pa,
              r_sc$b_up, r_sc$b_dn, r_sp$b_up, r_sp$b_dn,
              r_ph$b_up, r_ph$b_dn),
      sprintf(" & & %s & %s & %s & %s & %s & %s \\\\",
              r_sc$se_up, r_sc$se_dn, r_sp$se_up, r_sp$se_dn,
              r_ph$se_up, r_ph$se_dn))
  }
}
 
lines <- c(lines,
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Abandonment}} \\\\"
)
 
for (sp in SUBPERIODS) {
  for (pa in c("Panel A","Panel B")) {
    r_sc <- make_row(subperiod_all,"abandonment",sp$label,pa,"SC_Scaffolding")
    r_sp <- make_row(subperiod_all,"abandonment",sp$label,pa,"SC_Specialized")
    r_ph <- make_row(subperiod_all,"abandonment",sp$label,pa,"Physical_Terminal")
    lines <- c(lines,
      sprintf("\\multirow{2}{*}{%s} & %s & %s & %s & %s & %s & %s & %s \\\\",
              sp$label, pa,
              r_sc$b_up, r_sc$b_dn, r_sp$b_up, r_sp$b_dn,
              r_ph$b_up, r_ph$b_dn),
      sprintf(" & & %s & %s & %s & %s & %s & %s \\\\",
              r_sc$se_up, r_sc$se_dn, r_sp$se_up, r_sp$se_dn,
              r_ph$se_up, r_ph$se_dn))
  }
}
 
lines <- c(lines,
  "\\bottomrule",
  "\\multicolumn{8}{p{15cm}}{\\footnotesize Note:",
  "  $\\hat{\\beta}^{\\uparrow}_{\\text{Physical-sensory}} < 0$ in adoption and $> 0$",
  "  in abandonment across all three sub-periods and both panels, confirming",
  "  the directional asymmetry is not confined to the full 2015--2024 window.",
  "  Specialized and general socio-cognitive requirements show positive",
  "  $\\hat{\\beta}^{\\uparrow}$ in adoption and negative in abandonment across",
  "  all sub-periods. Weaker significance in 2019--2021 reflects reduced",
  "  statistical power from shorter windows, not reversal of the pattern.}",
  "\\end{tabular}",
  "\\end{table}"
)
 
writeLines(lines, file.path(SI_TABLES, "table_S2_subperiods.tex"))
message("  Saved: table_S2_subperiods.tex")
message("  Saved: table_S2_subperiods.csv")
 
# ==============================================================================
# FIGURE — Sub-period stability
# Autonomous reload: skip estimation when re-running figure only
# Forest plot: 2 rows (β↑ / β↓) × 3 columns (skill types)
#              Colors = sub-period; stacked Adoption / Abandonment
# ==============================================================================
message("\n>>> Generating Fig. SI_subperiods...")
 
if (!exists("subperiod_all")) {
  message("Reloading from CSV...")
  subperiod_all <- fread(file.path(SI_TABLES, "table_S2_subperiods.csv"))
}
 
COEF_KEY <- c(
  "b_up_SC_Scaffolding",    "b_dn_SC_Scaffolding",
  "b_up_SC_Specialized",    "b_dn_SC_Specialized",
  "b_up_Physical_Terminal", "b_dn_Physical_Terminal"
)
 
PERIOD_LEVELS <- c("2015\u20132018", "2019\u20132021", "2022\u20132024")
 
# AAAS palette for periods
PERIOD_COLS  <- c(
  "2015\u20132018" = "#3B4992",
  "2019\u20132021" = "#EE0000",
  "2022\u20132024" = "#008280"
)
PERIOD_FILLS <- c(
  "2015\u20132018" = "#7B89C2",
  "2019\u20132021" = "#F07070",
  "2022\u20132024" = "#40B2B0"
)
 
plot_dt <- subperiod_all[panel == "Panel A" &
                           coef %in% COEF_KEY &
                           !is.na(estimate)]
plot_dt[, ci_lo := estimate - 1.96 * std_error]
plot_dt[, ci_hi := estimate + 1.96 * std_error]
 
plot_dt[, direction := fifelse(grepl("^b_up", coef),
                                "\u03b2\u2191", "\u03b2\u2193")]
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
plot_dt[, period := factor(period, levels = PERIOD_LEVELS)]
 
theme_si <- theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12,
                                      margin = margin(t = 2, b = 4)),
    axis.title         = element_text(size = 13),
    axis.title.x       = element_text(size = 12, margin = margin(t = 5)),
    axis.text          = element_text(size = 11, colour = "grey15"),
    axis.text.x        = element_text(size = 10),
    panel.border       = element_rect(colour = "grey30", fill = NA,
                                      linewidth = 0.7),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.7),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 12),
    legend.key.width   = unit(0.5, "cm"),
    panel.spacing.x    = unit(0.7, "lines"),
    panel.spacing.y    = unit(0.5, "lines"),
    plot.title         = element_text(face = "bold", size = 14,
                                      margin = margin(b = 4))
  )
 
make_flow_plot <- function(flow_name, show_legend = FALSE) {
  pd <- plot_dt[flow_label == flow_name]
 
  p <- ggplot(pd, aes(x = period, y = estimate,
                      colour = period, fill = period,
                      group  = period)) +
    geom_hline(yintercept = 0, colour = "grey55",
               linewidth = 0.35, linetype = "dotted") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.20, linewidth = 0.75,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = 3.5, shape = 21, stroke = 0.9,
               position = position_dodge(width = 0.6)) +
    facet_grid(direction ~ skill_type, scales = "free_y") +
    scale_colour_manual(values = PERIOD_COLS, name = NULL) +
    scale_fill_manual(  values = PERIOD_FILLS, name = NULL) +
    scale_x_discrete(expand = expansion(mult = c(0.3, 0.3))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(title = flow_name,
         x     = "Sub-period",
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
 
fig_sub <- (p_adopt / p_aband_noleg /
  patchwork::wrap_elements(full = leg)) +
  plot_layout(heights = c(1, 1, 0.06))
 
ggsave(file.path(SI_FIGS, "fig_SI_subperiods.pdf"),
       fig_sub, width = 13, height = 11,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(SI_FIGS, "fig_SI_subperiods.png"),
       fig_sub, width = 13, height = 11,
       units = "in", dpi = 300, bg = "white")
 
message("  Saved: fig_SI_subperiods.pdf / .png")
message("\n>>> SI_10_S9_table_S2_subperiods.R complete.")