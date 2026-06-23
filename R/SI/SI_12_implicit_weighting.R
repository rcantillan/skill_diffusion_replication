# ==============================================================================
# SI_12_implicit_weighting.R  (v3 — frequency weights)
#
# PROBLEM. The standard riskset has one row per (target j, source i, skill s).
# A realized event with n_js eligible sources contributes n_js rows to the
# likelihood — implicitly weighting it n_js times. The skill FE absorbs the
# skill-level mean of n_js but not within-skill between-target variation in
# n_js, which is precisely where the gap slopes β↑/β↓ are identified from.
# If source-rich targets sit systematically at one end of the status hierarchy,
# β↑/β↓ could reflect row multiplicity rather than occupational behaviour.
#
# SOLUTION (option b — frequency weights). Re-estimate the standard dyadic
# model with weights = 1/n_js in feglm, so each realized event (j,s) counts
# equally regardless of how many sources it has. The total contribution of
# event (j,s) to the likelihood becomes Σ_i 1/n_js = 1.
#
# This is a one-argument change to the main model:
#   feglm(..., weights = ~w_ij, ...)
# Everything else — formula, FE structure, clustering, data — is identical
# to 05a/05b. Both Panel A (source + skill FE) and Panel B (target + skill FE)
# are estimated, preserving the full dyadic identification structure.
#
# The comparison is: weighted model vs. unweighted baseline (SM3/SM4).
# Any difference in coefficients is attributable solely to removing the
# implicit source-multiplicity weighting.
#
# Input:
#   data/derived/riskset_adoption.rds
#   data/derived/riskset_abandonment.rds
#   output/tables/main/occ_status_scores.csv
#   output/tables/main/coefs_pc1_adoption.csv    <- Panel A/B from 05a
#   output/tables/main/coefs_pc1_abandonment.csv <- Panel A/B from 05b
#
# Output:
#   output/tables/si/tab_SI20_weighted_coefs.csv
#   output/figures/si/fig_SI20_weighted.pdf
#   output/figures/si/fig_SI20_weighted.png
# ==============================================================================

library(data.table)
library(fixest)
library(ggplot2)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
DIR_DERIVED <- "data/derived"
DIR_SCORES  <- "output/tables/main"
DIR_TABLES  <- "output/tables/si"
DIR_FIGS    <- "output/figures/si"
dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGS,   showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Status index
# ------------------------------------------------------------------------------
scores_path <- file.path(DIR_SCORES, "occ_status_scores.csv")
stopifnot("Status index not found — run 04_status_pca.R first" =
            file.exists(scores_path))
sigma <- fread(scores_path, colClasses = c(occ = "character"))
stopifnot(all(c("occ", "status_pc1") %in% names(sigma)))
message(sprintf(">>> Status index: %d occupations", nrow(sigma)))

# ------------------------------------------------------------------------------
# prep_riskset()
#
# Loads riskset, attaches status gap variables (identical to 05a/05b),
# and computes w_ij = 1/n_js where n_js = number of sources for (target, skill).
# ------------------------------------------------------------------------------
prep_riskset <- function(flow = c("adoption", "abandonment")) {
  flow        <- match.arg(flow)
  rds_path    <- file.path(DIR_DERIVED, sprintf("riskset_%s.rds", flow))
  outcome_col <- if (flow == "adoption") "diffusion" else "abandonment"
  stopifnot(file.exists(rds_path))

  message(sprintf("\n>>> prep_riskset('%s')...", flow))
  dt <- readRDS(rds_path)
  setDT(dt)
  message(sprintf("  Raw riskset: %s rows", format(nrow(dt), big.mark = ",")))

  # Keep only needed columns
  keep <- c("source", "target", "skill_name", outcome_col,
            "atc_archetype", "structural_distance")
  missing_cols <- setdiff(keep, names(dt))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing columns in riskset: %s",
                 paste(missing_cols, collapse = ", ")))
  dt <- dt[, ..keep]
  dt <- dt[!is.na(structural_distance)]

  # Attach status scores
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  dt[sigma, on = .(source = occ), s_pc1 := i.status_pc1]
  dt[sigma, on = .(target = occ), t_pc1 := i.status_pc1]
  n_before <- nrow(dt)
  dt <- dt[!is.na(s_pc1) & !is.na(t_pc1)]
  if (n_before > nrow(dt))
    message(sprintf("  Dropped %d rows with unmatched status",
                    n_before - nrow(dt)))

  # Status gap — identical to 05a/05b
  dt[, G_ij      := t_pc1 - s_pc1]
  dt[, pc1_up    := pmax(0,  G_ij)]
  dt[, pc1_down  := pmin(0,  G_ij)]
  dt[, pc1_dummy := as.integer(G_ij > 0)]
  dt[, c("s_pc1", "t_pc1", "G_ij") := NULL]

  # Frequency weight: w_ij = 1/n_js
  # n_js = number of source rows per (target, skill_name)
  dt[, n_js := .N, by = .(target, skill_name)]
  dt[, w_ij := 1 / n_js]

  message(sprintf("  n_js: median=%.0f  mean=%.1f  max=%.0f",
                  median(dt$n_js), mean(dt$n_js), max(dt$n_js)))
  message(sprintf("  w_ij: median=%.4f  mean=%.4f  min=%.4f",
                  median(dt$w_ij), mean(dt$w_ij), min(dt$w_ij)))
  message(sprintf("  Final riskset: %s rows", format(nrow(dt), big.mark = ",")))

  list(data = dt, outcome = outcome_col, flow = flow)
}

# ------------------------------------------------------------------------------
# fit_weighted()
#
# Estimates cloglog model with w_ij = 1/n_js frequency weights.
# Formula and FE structure identical to 05a/05b.
# Both Panel A (source + skill FE) and Panel B (target + skill FE).
# ------------------------------------------------------------------------------
fit_weighted <- function(rs_obj, panel) {
  dt      <- rs_obj$data
  outcome <- rs_obj$outcome
  flow    <- rs_obj$flow

  fe1 <- if (panel == "A") "source" else "target"
  fe2 <- "skill_name"

  fml <- as.formula(sprintf(
    "%s ~ (pc1_dummy + pc1_up + pc1_down + structural_distance):atc_archetype",
    outcome
  ))
  cluster_fml <- ~ source + target + skill_name

  message(sprintf("  Fitting Panel %s | %s | weighted...", panel, flow))
  t0 <- proc.time()["elapsed"]
  m <- feglm(
    fml,
    data      = dt,
    family    = binomial(link = "cloglog"),
    fixef     = c(fe1, fe2),
    cluster   = cluster_fml,
    weights   = ~w_ij,          # <- frequency weight: 1/n_js
    lean      = TRUE,
    mem.clean = TRUE
  )
  elapsed <- round((proc.time()["elapsed"] - t0) / 60, 1)
  message(sprintf("  Done in %.1f min | N = %s",
                  elapsed, format(nobs(m), big.mark = ",")))
  m
}

# ------------------------------------------------------------------------------
# extract_coefs()
#
# Extracts β↑ (b_up) and β↓ (b_dn) per archetype from a fitted model.
# ------------------------------------------------------------------------------
extract_coefs <- function(m, flow_label, panel_label, model_label) {
  ct <- coeftable(m)
  dt <- data.table(term     = rownames(ct),
                   estimate = ct[, "Estimate"],
                   se       = ct[, "Std. Error"])

  archetypes <- c("SC_Specialized", "SC_Scaffolding", "Physical_Terminal")
  rows <- rbindlist(lapply(archetypes, function(arch) {
    find_val <- function(var_name) {
      # fixest can order interaction terms either way
      patterns <- c(
        sprintf("%s:atc_archetype%s", var_name, arch),
        sprintf("atc_archetype%s:%s", arch, var_name)
      )
      for (p in patterns) {
        idx <- which(dt$term == p)
        if (length(idx) > 0)
          return(c(dt$estimate[idx[1]], dt$se[idx[1]]))
      }
      c(NA_real_, NA_real_)
    }
    up   <- find_val("pc1_up")
    down <- find_val("pc1_down")
    data.table(
      flow      = flow_label,
      panel     = panel_label,
      model     = model_label,
      archetype = arch,
      term      = c("b_up", "b_dn"),
      estimate  = c(up[1],   down[1]),
      se        = c(up[2],   down[2])
    )
  }))

  rows[, sig := fcase(
    is.na(estimate),                    "",
    abs(estimate / se) >= 3.29,         "***",
    abs(estimate / se) >= 2.58,         "**",
    abs(estimate / se) >= 1.96,         "*",
    abs(estimate / se) >= 1.65,         "~",
    default = "ns"
  )]
  rows
}

# ------------------------------------------------------------------------------
# Main loop — adoption + abandonment, Panel A + Panel B
# ------------------------------------------------------------------------------
all_coefs <- list()

for (fl in c("adoption", "abandonment")) {
  message("\n", strrep("=", 60))
  message(sprintf("  %s", toupper(fl)))
  message(strrep("=", 60))

  rs <- prep_riskset(fl)

  for (pn in c("A", "B")) {
    m     <- fit_weighted(rs, panel = pn)
    label <- sprintf("Panel %s", pn)
    coefs <- extract_coefs(m, fl, label, "Weighted (1/n_js)")
    all_coefs[[paste(fl, pn)]] <- coefs

    message(sprintf("\n  Coefficients — %s Panel %s (weighted):", fl, pn))
    print(coefs[, .(archetype, term,
                    estimate = round(estimate, 3),
                    se       = round(se, 3), sig)])
    rm(m); gc()
  }
  rm(rs); gc()
}

coef_tbl <- rbindlist(all_coefs)

# ------------------------------------------------------------------------------
# Load baseline (unweighted) coefficients from 05a/05b for comparison
# Expected format: columns flow, panel (Panel A / Panel B), archetype,
#                  var (b_up / b_dn), coef, se
# ------------------------------------------------------------------------------
base_adopt <- file.path(DIR_SCORES, "coefs_pc1_adoption.csv")
base_aband <- file.path(DIR_SCORES, "coefs_pc1_abandonment.csv")

has_baseline <- file.exists(base_adopt) && file.exists(base_aband)

if (has_baseline) {
  bl_a <- fread(base_adopt)
  bl_b <- fread(base_aband)
  bl_a[, flow := "adoption"]
  bl_b[, flow := "abandonment"]
  baseline <- rbind(bl_a, bl_b)

  # Harmonise column names from 05a/05b extract_coefs format:
  # panel, var (b_up/b_dn), archetype, coef, se
  setnames(baseline,
           old = c("var",  "coef"),
           new = c("term", "estimate"),
           skip_absent = TRUE)
  baseline[, model := "Baseline (unweighted)"]
  baseline <- baseline[term %in% c("b_up", "b_dn")]

  compare <- rbind(
    baseline[, .(flow, panel, model, archetype, term, estimate, se)],
    coef_tbl[, .(flow, panel, model, archetype, term, estimate, se)]
  )
} else {
  warning("Baseline coefficient files not found — figure will show weighted only.")
  compare <- coef_tbl[, .(flow, panel, model, archetype, term, estimate, se)]
}

# ------------------------------------------------------------------------------
# Save coefficient table
# ------------------------------------------------------------------------------
out_path <- file.path(DIR_TABLES, "tab_SI20_weighted_coefs.csv")
fwrite(coef_tbl, out_path)
message(sprintf("\n  Saved: %s", out_path))

# ------------------------------------------------------------------------------
# Sign-match summary vs baseline
# ------------------------------------------------------------------------------
if (has_baseline) {
  message("\n>>> Sign-match summary (weighted vs. unweighted):")
  for (fl in c("adoption", "abandonment")) {
    for (pn in c("A", "B")) {
      label <- sprintf("Panel %s", pn)
      w  <- coef_tbl[flow == fl & panel == label,
                     .(archetype, term, est_w = estimate)]
      b  <- baseline[flow == fl & panel == label,
                     .(archetype, term, est_b = estimate)]
      cm <- merge(w, b, by = c("archetype", "term"))
      cm[, match := sign(est_w) == sign(est_b)]
      message(sprintf("  %s Panel %s: %d/%d sign preserved",
                      fl, pn, sum(cm$match, na.rm = TRUE), nrow(cm)))
    }
  }
}

# ------------------------------------------------------------------------------
# Comparison figure — forest plot
# ------------------------------------------------------------------------------
message("\n>>> Generating comparison figure...")

compare[, ci_lo := estimate - 1.96 * se]
compare[, ci_hi := estimate + 1.96 * se]

compare[, flow_lab := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  levels = c("Adoption", "Abandonment"))]

compare[, panel_lab := factor(
  fifelse(panel == "Panel A",
          "Panel A (source + skill FE)",
          "Panel B (target + skill FE)"),
  levels = c("Panel A (source + skill FE)",
             "Panel B (target + skill FE)"))]

compare[, coef_lab := fcase(
  term == "b_up" & archetype == "SC_Specialized",   "β↑ Spec. SC",
  term == "b_dn" & archetype == "SC_Specialized",   "β↓ Spec. SC",
  term == "b_up" & archetype == "SC_Scaffolding",   "β↑ Gen. SC",
  term == "b_dn" & archetype == "SC_Scaffolding",   "β↓ Gen. SC",
  term == "b_up" & archetype == "Physical_Terminal", "β↑ Physical",
  term == "b_dn" & archetype == "Physical_Terminal", "β↓ Physical"
)]
compare[, coef_lab := factor(coef_lab, levels = c(
  "β↓ Physical", "β↑ Physical",
  "β↓ Gen. SC",  "β↑ Gen. SC",
  "β↓ Spec. SC", "β↑ Spec. SC"
))]
compare[, model := factor(model, levels = c(
  "Baseline (unweighted)",
  "Weighted (1/n_js)"
))]

MODEL_COLS <- c(
  "Baseline (unweighted)" = "#3B4992",
  "Weighted (1/n_js)"     = "#EE0000"
)
MODEL_SHAPES <- c(
  "Baseline (unweighted)" = 16L,
  "Weighted (1/n_js)"     = 17L
)

theme_si <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10, colour = "grey15"),
    panel.border     = element_rect(colour = "grey30", fill = NA,
                                    linewidth = 0.6),
    axis.line        = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10),
    panel.spacing    = unit(1.0, "lines"),
    plot.title       = element_text(face = "bold", size = 12,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, colour = "grey30",
                                    margin = margin(b = 8))
  )

fig <- ggplot(compare,
              aes(x = estimate, y = coef_lab,
                  colour = model, shape = model)) +
  geom_vline(xintercept = 0, colour = "grey50",
             linewidth = 0.6, linetype = "dotted") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.22, linewidth = 0.8,
                 position = position_dodge(0.6)) +
  geom_point(size = 3.2, position = position_dodge(0.6)) +
  facet_grid(panel_lab ~ flow_lab, scales = "free_x") +
  scale_colour_manual(values = MODEL_COLS) +
  scale_shape_manual( values = MODEL_SHAPES) +
  labs(
    title    = "Robustness: implicit source-multiplicity weighting",
    subtitle = paste0(
      "Weighted model re-estimates the dyadic model with w_ij = 1/n_js, ",
      "so each (target, skill) event contributes equally to the likelihood.\n",
      "Formula and clustering identical to Tables SM3/SM4."),
    x = "Coefficient (cloglog scale)",
    y = NULL
  ) +
  theme_si

pdf_path <- file.path(DIR_FIGS, "fig_SI20_weighted.pdf")
png_path <- file.path(DIR_FIGS, "fig_SI20_weighted.png")

ggsave(pdf_path, fig, width = 10, height = 7,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(png_path, fig, width = 10, height = 7,
       units = "in", dpi = 300, bg = "white")
message(sprintf("  Saved: %s", pdf_path))
message(sprintf("  Saved: %s", png_path))

message("\n>>> SI_12_implicit_weighting.R (v3) complete.")
message(sprintf("    Weighted coefs: %s", out_path))