# ==============================================================================
# SI_08_archetype_misclassification.R
# Robustness to skill-archetype misclassification (10% and 20%)
#
# DESIGN:
#   For each contamination level p ∈ {0.10, 0.20}, B = 200 replications.
#
#   Each replication:
#     1. Selects floor(p × 160) skills without replacement.
#     2. For each selected skill, draws a new label from the empirical marginal
#        distribution — Multinomial(1, [49/160, 48/160, 63/160]).
#        A skill CAN receive its original label (no guaranteed misassignment).
#        Marginals are preserved in expectation across replications.
#     3. Propagates the new labels to all triads via join on skill_name.
#     4. Re-estimates the full model under two FE specifications:
#          Panel A: source + skill FE
#          Panel B: target + skill FE
#        Both models are estimated on the same perturbed dt_b.
#     5. Extracts β↑ and β↓ for each archetype in both panels and saves.
#
#   Summary: mean, SD, 5th/95th percentiles, sign preservation rate.
#   Figure: violin distributions vs. baseline reference diamonds — Panels A and B.
#
# OPTIMIZATIONS:
#   - No clustering (vcov = "iid"): sign preservation uses point estimates only.
#     Clustering affects SEs, not β. Saves ~40% per replication.
#   - only.coef = TRUE: returns named vector, not a full fixest object.
#   - Both panels per replication on the same dt_b: single copy of the data.
#   - Checkpoint every CKPT_EVERY replications — safe to interrupt and resume.
#
# NOTE: if checkpoints from a prior Panel-A-only run exist, remove them first:
#         rm -f output/models/si/missclass_ckpt_*.rds
#
# Input:
#   data/derived/riskset_adoption.rds
#   data/derived/riskset_abandonment.rds
#   output/tables/main/occ_status_scores.csv
#   output/tables/si/baseline_coefs_2d.csv   (optional)
#
# Output:
#   output/tables/si/test_archetype_misclassification.csv
#   output/tables/si/tab_archetype_misclassification_summary.csv
#   output/figures/si/fig_archetype_misclassification.pdf / .png
#
# LAUNCH:
#   cd /path/to/skill_diffusion/
#   rm -f output/models/si/missclass_ckpt_*.rds
#   mkdir -p logs
#   systemd-run --user --scope \
#     -p MemoryMax=12G -p MemoryHigh=11G -p CPUWeight=50 \
#     nice -n 10 ionice -c3 \
#     Rscript --vanilla R/SI/SI_08_archetype_misclassification.R \
#     > "logs/si08_$(date +%Y%m%d_%H%M).log" 2>&1 &
#   disown
# ==============================================================================

t_global_start <- proc.time()
source("R/SI/00_setup_SI.R")
library(ggplot2)
library(patchwork)

# ==============================================================================
# Parameters
# ==============================================================================
CONTAM_RATES <- c(0.10, 0.20)    # 10% and 20% contamination rates
B_MISSCLASS  <- 200L
CKPT_EVERY   <- 25L

SKILL_LEVELS <- c("SC_General", "SC_Specialized", "Physical_Terminal")
SKILL_LABELS <- c("General socio-cognitive",
                  "Specialized socio-cognitive",
                  "Sensory-physical")

COEF_KEY <- c(
  "b_up_SC_General",    "b_dn_SC_General",
  "b_up_SC_Specialized",    "b_dn_SC_Specialized",
  "b_up_Physical_Terminal", "b_dn_Physical_Terminal"
)

TERM_MAP <- list(
  b_up_SC_General    = "pc1_up:atc_archetypeSC_General",
  b_dn_SC_General    = "pc1_down:atc_archetypeSC_General",
  b_up_SC_Specialized    = "pc1_up:atc_archetypeSC_Specialized",
  b_dn_SC_Specialized    = "pc1_down:atc_archetypeSC_Specialized",
  b_up_Physical_Terminal = "pc1_up:atc_archetypePhysical_Terminal",
  b_dn_Physical_Terminal = "pc1_down:atc_archetypePhysical_Terminal"
)

baseline_path <- file.path(SI_TABLES, "baseline_coefs_3skill.csv")
baseline <- if (file.exists(baseline_path)) {
  fread(baseline_path)
} else {
  message("  [note] baseline_coefs_3skill.csv not found — reference diamonds omitted")
  NULL
}

# ==============================================================================
# elapsed_str()
# ==============================================================================
elapsed_str <- function(t0) {
  s <- as.numeric(proc.time() - t0)[3]
  sprintf("%02dh %02dm %02ds",
          floor(s / 3600), floor((s %% 3600) / 60), round(s %% 60))
}

# ==============================================================================
# load_flow_missclass()
# ==============================================================================
load_flow_missclass <- function(flow = c("adoption", "abandonment"),
                                seed = SEED,
                                frac = SAMPLE_FRAC) {
  flow <- match.arg(flow)
  message(sprintf("\n>>> load_flow_missclass('%s') | %s",
                  flow, format(Sys.time(), "%H:%M:%S")))

  rds_path    <- if (flow == "adoption")
                   "data/derived/riskset_adoption.rds" else
                   "data/derived/riskset_abandonment.rds"
  scores_path <- "output/tables/main/occ_status_scores.csv"
  outcome_col <- if (flow == "adoption") "diffusion" else "abandonment"
  stopifnot(file.exists(rds_path), file.exists(scores_path))

  dt <- readRDS(rds_path); setDT(dt)
  keep <- c("source", "target", "skill_name", outcome_col,
            "atc_archetype", "structural_distance")
  dt <- dt[, ..keep]
  dt[, atc_archetype := factor(as.character(atc_archetype),
                                levels = SKILL_LEVELS)]
  dt <- dt[!is.na(atc_archetype) & !is.na(structural_distance)]

  scores <- fread(scores_path); setDT(scores)
  scores[, occ := as.character(occ)]
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  dt[scores, on = .(source = occ), s_pc1 := i.status_pc1]
  dt[scores, on = .(target = occ), t_pc1 := i.status_pc1]
  dt <- dt[!is.na(s_pc1) & !is.na(t_pc1)]
  dt[, pc1_gap  := t_pc1 - s_pc1]
  dt[, pc1_up   := pmax(0,  pc1_gap)]
  dt[, pc1_down := pmin(0,  pc1_gap)]
  dt[, up_dummy := fifelse(pc1_gap > 0, 1L, 0L)]
  dt[, c("s_pc1", "t_pc1", "pc1_gap") := NULL]
  rm(scores); gc()

  # Canonical skill → archetype map
  skill_map <- unique(dt[, .(skill_name,
                              archetype_orig = as.character(atc_archetype))])
  stopifnot(skill_map[, uniqueN(archetype_orig),
                      by = skill_name][, all(V1 == 1)])

  # Empirical marginal proportions for Multinomial sampling
  marginal_props <- skill_map[, .N, by = archetype_orig]
  marginal_props[, prop := N / sum(N)]
  setorder(marginal_props, archetype_orig)

  message(sprintf("  Skills: %d | Marginal proportions:",
                  nrow(skill_map)))
  for (i in seq_len(nrow(marginal_props)))
    message(sprintf("    %s: %d skills (%.1f%%)",
                    marginal_props$archetype_orig[i],
                    marginal_props$N[i],
                    marginal_props$prop[i] * 100))

  # 50% subsample — mismo seed para comparabilidad con baseline
  ckpt_src <- file.path(SI_MODELS, sprintf("sources_seed%d.rds", seed))
  if (file.exists(ckpt_src)) {
    sources_sample <- readRDS(ckpt_src)
    message(sprintf("  Sources from disk: %d", length(sources_sample)))
  } else {
    set.seed(seed)
    sources_sample <- sample(unique(dt$source),
                             size = round(uniqueN(dt$source) * frac))
    saveRDS(sources_sample, ckpt_src)
    message(sprintf("  Sources generated: %d of %d",
                    length(sources_sample), uniqueN(dt$source)))
  }
  dt <- dt[source %in% sources_sample]
  message(sprintf("  %s sample: %s triadas | %d src | %d tgt",
                  flow, format(nrow(dt), big.mark = ","),
                  uniqueN(dt$source), uniqueN(dt$target)))

  list(dt             = dt,
       outcome        = outcome_col,
       flow           = flow,
       skill_map      = skill_map,
       marginal_props = marginal_props)
}

# ==============================================================================
# extract_coefs()
# Extrae coeficientes desde el named vector de only.coef = TRUE.
# ==============================================================================
extract_coefs <- function(coef_vec, panel_label, flow_label) {
  rows <- lapply(names(TERM_MAP), function(cname) {
    term <- TERM_MAP[[cname]]
    val  <- if (term %in% names(coef_vec)) coef_vec[[term]] else NA_real_
    data.table(coef = cname, estimate = val,
               panel = panel_label, flow = flow_label)
  })
  rbindlist(rows)
}

# ==============================================================================
# perturb_labels()
# Selects floor(p × N_skills) skills and assigns each a new label
# drawn from the empirical marginal distribution (Multinomial).
# A skill CAN receive its original label.
# ==============================================================================
perturb_labels <- function(skill_map, marginal_props, p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  map      <- copy(skill_map)
  n_contam <- floor(p * nrow(map))
  idx      <- sample(nrow(map), n_contam, replace = FALSE)

  levels_ord <- marginal_props[order(archetype_orig), archetype_orig]
  props_ord  <- marginal_props[order(archetype_orig), prop]
  new_labels <- sample(levels_ord, size = n_contam,
                       replace = TRUE, prob = props_ord)

  map[, archetype_new := archetype_orig]
  map[idx, archetype_new := new_labels]
  map[]
}

# ==============================================================================
# run_missclass()
# B replications for one flow × contamination level p.
# Each replication estimates Panel A (source+skill FE) AND Panel B (target+skill FE)
# on the same perturbed dt_b — a single copy of the data per replication.
# ==============================================================================
run_missclass <- function(setup, p, B = B_MISSCLASS) {
  dt             <- setup$dt
  outcome        <- setup$outcome
  flow           <- setup$flow
  skill_map      <- setup$skill_map
  marginal_props <- setup$marginal_props

  n_contam <- floor(p * nrow(skill_map))
  tag  <- sprintf("%s_p%02d", flow, as.integer(p * 100))
  ckpt <- file.path(SI_MODELS, sprintf("missclass_ckpt_%s.rds", tag))

  done    <- if (file.exists(ckpt)) readRDS(ckpt) else list()
  start_b <- length(done) + 1L

  if (start_b > B) {
    message(sprintf("  [skip] %s already complete (%d reps)", tag, B))
    return(rbindlist(done))
  }

  message(sprintf(
    "\n  %s | p=%.0f%% (%d/%d skills) | Panel A + B | reps %d-%d | %s",
    flow, p * 100, n_contam, nrow(skill_map),
    start_b, B, format(Sys.time(), "%H:%M:%S")))

  fml <- as.formula(sprintf(
    "%s ~ (up_dummy + pc1_up + pc1_down + structural_distance) : atc_archetype",
    outcome))

  t_run <- proc.time()

  for (b in start_b:B) {

    if (b == start_b || b %% 25 == 0L || b == B)
      message(sprintf("    rep %d/%d | %s | %s",
                      b, B, elapsed_str(t_run),
                      format(Sys.time(), "%H:%M:%S")))

    b_seed <- b * 1000L + as.integer(p * 100) +
              ifelse(flow == "adoption", 0L, 50000L)

    # 1. Perturba etiquetas a nivel skill
    map_b <- perturb_labels(skill_map, marginal_props, p, seed = b_seed)

    # 2. Propaga a todas las triadas via join
    dt_b <- copy(dt)
    dt_b[map_b, on = "skill_name",
         atc_archetype := factor(i.archetype_new, levels = SKILL_LEVELS)]

    # 3a. Panel A — source + skill FE
    coef_A <- tryCatch(
      feglm(fml, data = dt_b,
            family    = binomial("cloglog"),
            fixef     = c("source", "skill_name"),
            vcov      = "iid",
            lean      = TRUE, mem.clean = TRUE,
            nthreads  = getFixest_nthreads(),
            only.coef = TRUE),
      error = function(e) {
        message(sprintf("    [!] rep %d Panel A failed: %s",
                        b, conditionMessage(e)))
        NULL
      }
    )

    # 3b. Panel B — target + skill FE (mismo dt_b)
    coef_B <- tryCatch(
      feglm(fml, data = dt_b,
            family    = binomial("cloglog"),
            fixef     = c("target", "skill_name"),
            vcov      = "iid",
            lean      = TRUE, mem.clean = TRUE,
            nthreads  = getFixest_nthreads(),
            only.coef = TRUE),
      error = function(e) {
        message(sprintf("    [!] rep %d Panel B failed: %s",
                        b, conditionMessage(e)))
        NULL
      }
    )

    # 4. Extrae y combina ambos paneles
    res_list <- list()
    if (!is.null(coef_A)) {
      res_A <- extract_coefs(coef_A, "Panel A", flow)
      res_A[, `:=`(rep = b, p_contam = p)]
      res_list[["A"]] <- res_A
    }
    if (!is.null(coef_B)) {
      res_B <- extract_coefs(coef_B, "Panel B", flow)
      res_B[, `:=`(rep = b, p_contam = p)]
      res_list[["B"]] <- res_B
    }

    if (length(res_list) > 0L)
      done[[b]] <- rbindlist(res_list)

    rm(coef_A, coef_B, dt_b, map_b, res_list); gc()

    if (b %% CKPT_EVERY == 0L || b == B) {
      saveRDS(done, ckpt)
      message(sprintf("    checkpoint: %d/%d reps | %s",
                      b, B, elapsed_str(t_run)))
    }
  }

  message(sprintf("  [DONE] %s p=%.0f%% — %s",
                  flow, p * 100, elapsed_str(t_run)))
  rbindlist(done)
}

# ==============================================================================
# Main
# ==============================================================================
reps_path <- file.path(SI_TABLES, "test_archetype_misclassification.csv")

if (file.exists(reps_path)) {
  message("  [cache] Loading saved replications from: ", reps_path)
  reps <- fread(reps_path)
  message(sprintf("  Loaded %s rows.", format(nrow(reps), big.mark = ",")))
} else {
  setup_a <- load_flow_missclass("adoption")
  setup_b <- load_flow_missclass("abandonment")

  all_reps <- list()

  for (p in CONTAM_RATES) {
    message(sprintf("\n%s\n>>> Contamination: %.0f%% | %s\n%s",
                    strrep("=", 60), p * 100,
                    format(Sys.time(), "%H:%M:%S"),
                    strrep("=", 60)))
    all_reps[[sprintf("adopt_p%02d", as.integer(p * 100))]] <-
      run_missclass(setup_a, p)
    all_reps[[sprintf("aband_p%02d", as.integer(p * 100))]] <-
      run_missclass(setup_b, p)
  }

  reps <- rbindlist(all_reps)
  fwrite(reps, reps_path)
  message("  Saved: test_archetype_misclassification.csv")
}

# ==============================================================================
# Summary — sign preservation and relative attenuation vs. baseline
# ==============================================================================
message("\n>>> Summary table...")

summ <- reps[coef %in% COEF_KEY & !is.na(estimate),
  .(mean_est = mean(estimate),
    sd_est   = sd(estimate),
    q05      = quantile(estimate, 0.05),
    q95      = quantile(estimate, 0.95),
    n_reps   = .N),
  by = .(flow, panel, p_contam, coef)]

reps[coef %in% COEF_KEY & !is.na(estimate),
     sign_expected := sign(mean(estimate)),
     by = .(flow, panel, p_contam, coef)]

summ_sign <- reps[coef %in% COEF_KEY & !is.na(estimate),
  .(sign_pres = mean(sign(estimate) == sign_expected)),
  by = .(flow, panel, p_contam, coef)]

summ <- merge(summ, summ_sign, by = c("flow", "panel", "p_contam", "coef"))

if (!is.null(baseline)) {
  base_ref <- baseline[coef %in% COEF_KEY,
                       .(flow, panel, coef, base_est = estimate)]
  summ     <- merge(summ, base_ref,
                    by = c("flow", "panel", "coef"), all.x = TRUE)
  summ[!is.na(base_est) & base_est != 0,
       pct_attenuation := round(100 * (1 - abs(mean_est) / abs(base_est)), 1)]
}

setorder(summ, flow, panel, p_contam, coef)
fwrite(summ, file.path(SI_TABLES,
                       "tab_archetype_misclassification_summary.csv"))
message("  Saved: tab_archetype_misclassification_summary.csv")

# ==============================================================================
# Figura — violin + diamonds baseline
# Filas: direction (β↑ / β↓)   Columnas: 10% / 20% mislabeled
# Secciones: flow × panel (Adoption A / Adoption B / Abandonment A / Abandonment B)
# ==============================================================================
message("\n>>> Generando figura...")

plot_dt <- reps[coef %in% COEF_KEY & !is.na(estimate)]
plot_dt[, direction := fifelse(grepl("^b_up", coef),
                               "\u03b2\u2191", "\u03b2\u2193")]
plot_dt[, skill_type := fcase(
  grepl("SC_General",    coef), "Gen.\nsocio-cog.",
  grepl("SC_Specialized",    coef), "Spec.\nsocio-cog.",
  grepl("Physical_Terminal", coef), "Sensory-\nphysical")]
plot_dt[, skill_type := factor(skill_type,
  levels = c("Gen.\nsocio-cog.", "Spec.\nsocio-cog.", "Sensory-\nphysical"))]
plot_dt[, direction  := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
plot_dt[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]
plot_dt[, contam_label := factor(
  sprintf("%.0f%% mislabeled", p_contam * 100),
  levels = sprintf("%.0f%% mislabeled", sort(CONTAM_RATES) * 100))]
plot_dt[, section := sprintf("%s — %s",
                             fifelse(flow == "adoption",
                                     "Adoption", "Abandonment"),
                             panel)]

theme_si <- theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    panel.border       = element_rect(colour = "grey30", fill = NA,
                                      linewidth = 0.5),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 11))

FLOW_COLS  <- c("Adoption" = "#1B4F8A", "Abandonment" = "#B03030")
FLOW_FILLS <- c("Adoption" = "#4A90D9", "Abandonment" = "#E07070")

make_section <- function(flow_name, panel_name) {
  pd <- plot_dt[flow_label == flow_name & panel == panel_name]
  if (nrow(pd) == 0L) return(NULL)

  p <- ggplot(pd, aes(x = skill_type, y = estimate,
                      fill = flow_label, colour = flow_label)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_violin(alpha = 0.35, linewidth = 0.4,
                draw_quantiles = c(0.05, 0.50, 0.95)) +
    geom_boxplot(width = 0.10, outlier.shape = NA, alpha = 0.8,
                 colour = FLOW_COLS[flow_name], fill = "white",
                 linewidth = 0.4) +
    facet_grid(direction ~ contam_label, scales = "free_y") +
    scale_fill_manual(  values = FLOW_FILLS) +
    scale_colour_manual(values = FLOW_COLS) +
    labs(title = sprintf("%s — %s", flow_name, panel_name),
         x = NULL, y = "Estimate (cloglog)") +
    theme_si

  # Diamonds del baseline
  if (!is.null(baseline)) {
    base_ref <- baseline[panel == panel_name &
                           flow  == tolower(flow_name) &
                           coef  %in% COEF_KEY]
    if (nrow(base_ref) > 0L) {
      base_ref[, direction := fifelse(grepl("^b_up", coef),
                                      "\u03b2\u2191", "\u03b2\u2193")]
      base_ref[, skill_type := fcase(
        grepl("SC_General",    coef), "Gen.\nsocio-cog.",
        grepl("SC_Specialized",    coef), "Spec.\nsocio-cog.",
        grepl("Physical_Terminal", coef), "Sensory-\nphysical")]
      base_ref[, skill_type := factor(skill_type,
        levels = c("Gen.\nsocio-cog.", "Spec.\nsocio-cog.",
                   "Sensory-\nphysical"))]
      base_ref[, direction := factor(direction,
        levels = c("\u03b2\u2191", "\u03b2\u2193"))]

      base_rep <- rbindlist(lapply(
        sprintf("%.0f%% mislabeled", sort(CONTAM_RATES) * 100),
        function(cl) {
          tmp <- copy(base_ref)
          tmp[, contam_label := factor(
            cl, levels = sprintf("%.0f%% mislabeled",
                                 sort(CONTAM_RATES) * 100))]
          tmp
        }
      ))
      p <- p +
        geom_point(data = base_rep,
                   aes(x = skill_type, y = estimate),
                   shape = 23, size = 3, stroke = 0.8,
                   colour = FLOW_COLS[flow_name], fill = "white",
                   inherit.aes = FALSE)
    }
  }
  p
}

# 4 secciones: Adoption A / Adoption B / Abandonment A / Abandonment B
s_AA <- make_section("Adoption",    "Panel A")
s_AB <- make_section("Adoption",    "Panel B")
s_BA <- make_section("Abandonment", "Panel A")
s_BB <- make_section("Abandonment", "Panel B")

fig <- s_AA / s_AB / s_BA / s_BB +
  plot_layout(heights = c(1, 1, 1, 1))

ggsave(file.path(SI_FIGS, "fig_archetype_misclassification.pdf"),
       fig, width = 12, height = 20, units = "in",
       device = cairo_pdf, bg = "white")
ggsave(file.path(SI_FIGS, "fig_archetype_misclassification.png"),
       fig, width = 12, height = 20, units = "in",
       dpi = 300, bg = "white")
message("  Saved: fig_archetype_misclassification.pdf / .png")

# ==============================================================================
# Console summary
# ==============================================================================
message("\n", strrep("=", 70))
message("  ARCHETYPE MISCLASSIFICATION ROBUSTNESS — SUMMARY")
message(strrep("=", 70))
message("  sign_pres = fraction of 200 replications where β retains its sign.")
message("  Values >= 0.90 across both panels and contamination levels confirm robustness.")
message(strrep("-", 70))

for (fl in c("adoption", "abandonment")) {
  for (pnl in c("Panel A", "Panel B")) {
    for (p in CONTAM_RATES) {
      message(sprintf("\n  %s | %s | contamination = %.0f%%",
                      toupper(fl), pnl, p * 100))
      tbl <- summ[flow == fl & panel == pnl & p_contam == p &
                    coef %in% c("b_up_SC_General",
                                "b_up_SC_Specialized",
                                "b_up_Physical_Terminal"),
                  .(coef,
                    mean      = round(mean_est,  4),
                    q05       = round(q05,       4),
                    q95       = round(q95,       4),
                    sign_pres = round(sign_pres, 3))]
      if ("pct_attenuation" %in% names(summ)) {
        atten <- summ[flow == fl & panel == pnl & p_contam == p &
                        coef %in% c("b_up_SC_General",
                                    "b_up_SC_Specialized",
                                    "b_up_Physical_Terminal"),
                      .(coef, pct_atten = round(pct_attenuation, 1))]
        tbl[atten, on = "coef", pct_atten := i.pct_atten]
      }
      print(tbl, row.names = FALSE)
    }
  }
}

message(sprintf(
  "\n>>> SI_08_archetype_misclassification.R complete | %s | %s",
  elapsed_str(t_global_start),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
message("    Next: SI_09_alternative_status.R")