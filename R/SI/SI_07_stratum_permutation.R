# ==============================================================================
# SI_07_stratum_permutation.R  —  Test (iii): Within-stratum permutation
#
# Permutes the outcome (diffusion / abandonment) within cells defined by
# atc_archetype × status_q_source. Preserves:
#   - Outcome rate by skill type (3 levels)
#   - Outcome rate by source status quintile
#   - Skill composition within each cell
# Only destroys the within-stratum association between gap direction and outcome.
#
# Skill types (3 levels):
#   SC_Scaffolding    — specialized socio-cognitive
#   SC_Specialized    — general socio-cognitive
#   Physical_Terminal — physical-sensory
#
# If the directional friction coefficients collapse under permutation, the
# directional signal is genuine and not a compositional artifact within strata.
#
# This test additionally addresses the correlated-technological-shocks
# alternative: if a common exogenous force pushed occupations toward the same
# skills independently of network position, source identity would be
# uninformative and permuted coefficients would match observed ones.
#
# CONSERVATIVE design: permutation acts at the observation (triad) level
# within strata, not at the dyad level — it does not perfectly preserve
# dependence between triads of the same dyad. This makes the test stricter.
#
# B = 1000 replications. Panel A only (source + skill FE) — primary
# identification for both flows (TARGET is the focal unit in both risk sets).
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
#         output/tables/si/baseline_coefs_2d.csv
#
# Output: output/tables/si/test_iii_stratum_perm.csv
#         output/tables/si/test_iii_null_stats.csv
#         output/figures/si/fig_SI_test_iii.pdf / .png
# ==============================================================================

source("R/SI/00_setup_SI.R")
library(ggplot2)
library(patchwork)

# Three skill types — ordered for factor levels and figure facets
SKILL_LEVELS <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")
SKILL_LABELS <- c("Specialized socio-cognitive",
                  "General socio-cognitive",
                  "Physical-sensory")

COEF_KEY <- c(
  "b_up_SC_Scaffolding",    "b_dn_SC_Scaffolding",
  "b_up_SC_Specialized",    "b_dn_SC_Specialized",
  "b_up_Physical_Terminal", "b_dn_Physical_Terminal"
)

baseline <- fread(file.path(SI_TABLES, "baseline_coefs_2d.csv"))

# ==============================================================================
# extract_coefs_3skill()
# Extracts β↑ and β↓ for each of the three skill types from a feglm model.
# Uses exact string matching — avoids grepl() mishandling ':' in term names.
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

  term_map <- list(
    b_up_SC_Scaffolding    = "pc1_up:atc_archetypeSC_Scaffolding",
    b_dn_SC_Scaffolding    = "pc1_down:atc_archetypeSC_Scaffolding",
    b_up_SC_Specialized    = "pc1_up:atc_archetypeSC_Specialized",
    b_dn_SC_Specialized    = "pc1_down:atc_archetypeSC_Specialized",
    b_up_Physical_Terminal = "pc1_up:atc_archetypePhysical_Terminal",
    b_dn_Physical_Terminal = "pc1_down:atc_archetypePhysical_Terminal"
  )

  rows <- lapply(names(term_map), function(coef_name) {
    pattern <- term_map[[coef_name]]
    idx     <- which(ct[["term"]] == pattern)   # exact match
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
# attach_archetype()
# load_flow() returns 'domain' (Cognitive/Physical) but not 'atc_archetype'.
# Build lookup from riskset_adoption.rds and merge by skill_name.
# ==============================================================================
message(">>> Building atc_archetype lookup from riskset_adoption.rds...")
archetype_lkp <- unique(
  readRDS("data/derived/riskset_adoption.rds")[,
    .(skill_name,
      atc_archetype = as.character(atc_archetype))]
)[!is.na(atc_archetype)]
stopifnot(uniqueN(archetype_lkp$skill_name) == nrow(archetype_lkp))
message(sprintf("  Lookup: %d skills | SC_Scaffolding=%d | SC_Specialized=%d | Physical_Terminal=%d",
                nrow(archetype_lkp),
                archetype_lkp[atc_archetype == "SC_Scaffolding",    .N],
                archetype_lkp[atc_archetype == "SC_Specialized",    .N],
                archetype_lkp[atc_archetype == "Physical_Terminal", .N]))

attach_archetype <- function(setup) {
  dt <- setup$dt
  dt[archetype_lkp, on = "skill_name", atc_archetype := i.atc_archetype]
  dt[, atc_archetype := factor(atc_archetype, levels = SKILL_LEVELS)]
  dt <- dt[!is.na(atc_archetype)]
  setup$dt <- dt
  setup
}

# ==============================================================================
# PART 1 — ESTIMATION (skip if test_iii_stratum_perm.csv already exists)
# ==============================================================================
perm_csv <- file.path(SI_TABLES, "test_iii_stratum_perm.csv")

if (file.exists(perm_csv)) {

  message(">>> [SKIP ESTIMATION] Found existing ", perm_csv)
  message(">>> Loading saved permutation results...")
  perm_all <- fread(perm_csv)
  message(sprintf("    %d rows loaded | %d unique replications",
                  nrow(perm_all), uniqueN(perm_all$b_rep)))

} else {

  # --------------------------------------------------------------------------
  # run_stratum_perm()
  # Permutes outcome within (atc_archetype × status_q_source) cells.
  # Formula interacts spline terms with atc_archetype (3 levels).
  # --------------------------------------------------------------------------
  run_stratum_perm <- function(setup) {

    dt         <- setup$dt
    flow_label <- setup$flow
    outcome    <- setup$outcome
    ckpt_file  <- file.path(SI_MODELS,
                             sprintf("test_iii_%s_ckpt.rds", flow_label))

    # Cell distribution diagnostic
    cell_dist <- dt[, .N, by = .(atc_archetype, status_q_source)]
    message("  Cell distribution (atc_archetype x status_q_source):")
    print(cell_dist[order(atc_archetype, status_q_source)])
    message(sprintf("  Min obs per cell: %d", min(cell_dist$N)))

    # Checkpoint resume
    if (file.exists(ckpt_file)) {
      done    <- readRDS(ckpt_file)
      done    <- done[!is.na(estimate)]   # discard stale NAs
      if (nrow(done) == 0L) {
        file.remove(ckpt_file)
        b_start <- 1L
        results <- list()
      } else {
        b_start <- max(done$b_rep) + 1L
        results <- list(prev = done)
        message(sprintf("  [RESUME] %s: from b=%d", flow_label, b_start))
      }
    } else {
      results <- list()
      b_start <- 1L
    }

    # Formula — interactions with atc_archetype (3 levels)
    fml <- as.formula(sprintf(
      "%s ~ (up_dummy + pc1_up + pc1_down + structural_distance) : atc_archetype",
      outcome))

    for (b in seq(b_start, B_PERM)) {
      if (b %% 100 == 0 || b == b_start)
        message(sprintf("  %s Panel A | b=%d/%d | %s",
                        flow_label, b, B_PERM,
                        format(Sys.time(), "%H:%M:%S")))

      set.seed(SEED * 1000L + b)

      # Permute outcome within (atc_archetype × status_q_source)
      dt[, outcome_orig := get(outcome)]
      dt[, (outcome) := sample(get(outcome)),
         by = .(atc_archetype, status_q_source)]

      m <- feglm(
        fml,
        data      = dt,
        family    = binomial("cloglog"),
        fixef     = c("source", "skill_name"),
        cluster   = c("source", "target", "skill_name"),
        lean      = TRUE, mem.clean = TRUE, nthreads = 0
      )

      res        <- extract_coefs_3skill(m, "Panel A", flow_label)
      res[, b_rep := b]
      results[[as.character(b)]] <- res

      # Restore original outcome
      dt[, (outcome) := outcome_orig]
      dt[, outcome_orig := NULL]
      rm(m); gc()

      # Checkpoint every 50 replicas
      if (b %% 50 == 0)
        saveRDS(rbindlist(results), ckpt_file)
    }

    out <- rbindlist(results)
    if (file.exists(ckpt_file)) file.remove(ckpt_file)
    out
  }

  all_perm <- list()

  message("\n", strrep("=", 60))
  message(">>> ADOPTION — Test (iii)")
  message(strrep("=", 60))
  setup_a <- attach_archetype(load_flow("adoption"))
  all_perm[["adopt"]] <- run_stratum_perm(setup_a)
  rm(setup_a); gc()

  message("\n", strrep("=", 60))
  message(">>> ABANDONMENT — Test (iii)")
  message(strrep("=", 60))
  setup_b <- attach_archetype(load_flow("abandonment"))
  all_perm[["aband"]] <- run_stratum_perm(setup_b)
  rm(setup_b); gc()

  perm_all <- rbindlist(all_perm)
  fwrite(perm_all, perm_csv)
  message("  Saved: test_iii_stratum_perm.csv")
}

# ==============================================================================
# PART 2 — NULL DISTRIBUTION STATISTICS
# ==============================================================================
null_stats <- perm_all[coef %in% COEF_KEY & !is.na(estimate), .(
  null_mean = mean(estimate, na.rm = TRUE),
  null_sd   = sd(estimate,   na.rm = TRUE),
  ci_lo_95  = quantile(estimate, 0.025, na.rm = TRUE),
  ci_hi_95  = quantile(estimate, 0.975, na.rm = TRUE),
  ci_lo_90  = quantile(estimate, 0.05,  na.rm = TRUE),
  ci_hi_90  = quantile(estimate, 0.95,  na.rm = TRUE),
  n_reps    = .N
), by = .(flow, coef)]

# Merge baseline (Panel A)
null_stats <- merge(
  null_stats,
  baseline[panel == "Panel A" & coef %in% COEF_KEY,
           .(flow, coef, base_est = estimate)],
  by = c("flow", "coef"),
  all.x = TRUE
)
null_stats[, sep_sds := (base_est - null_mean) / null_sd]
null_stats[, p_val   := pmin(1, 2 * pnorm(-abs(sep_sds)))]

fwrite(null_stats, file.path(SI_TABLES, "test_iii_null_stats.csv"))

message("\n>>> Baseline separation from null (in null SDs):")
print(null_stats[order(flow, coef),
                 .(flow, coef,
                   base_est  = round(base_est,  3),
                   null_mean = round(null_mean, 4),
                   null_sd   = round(null_sd,   4),
                   sep_sds   = round(sep_sds,   0),
                   p_lt_001  = p_val < 0.001)])

# ==============================================================================
# PART 3 — FIGURE S7
# Histogram of permuted null distributions with baseline vertical lines.
# Layout: 2 rows (direction β↑ / β↓) × 3 columns (skill types), per flow.
# ==============================================================================
message("\n>>> Generating Fig. S7...")

# Autonomous reload block
if (!exists("perm_all")) {
  message("Reloading from CSV...")
  perm_all   <- fread(file.path(SI_TABLES, "test_iii_stratum_perm.csv"))
  baseline   <- fread(file.path(SI_TABLES, "baseline_coefs_2d.csv"))
  null_stats <- fread(file.path(SI_TABLES, "test_iii_null_stats.csv"))
}

plot_dt <- perm_all[coef %in% COEF_KEY & !is.na(estimate)]

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

# Baseline reference lines (Panel A only — Test iii runs Panel A)
base_ref <- baseline[panel == "Panel A" & coef %in% COEF_KEY]
base_ref[, direction := fifelse(grepl("^b_up", coef), "\u03b2\u2191", "\u03b2\u2193")]
base_ref[, skill_type := fcase(
  grepl("SC_Scaffolding",    coef), "Specialized socio-cognitive",
  grepl("SC_Specialized",    coef), "General socio-cognitive",
  grepl("Physical_Terminal", coef), "Physical-sensory"
)]
base_ref[, skill_type := factor(skill_type,
  levels = c("Specialized socio-cognitive",
             "General socio-cognitive",
             "Physical-sensory"))]
base_ref[, direction  := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
base_ref[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]

# Null stats for SD annotation
null_ann <- null_stats[coef %in% COEF_KEY & !is.na(base_est)]
null_ann[, direction := fifelse(grepl("^b_up", coef), "\u03b2\u2191", "\u03b2\u2193")]
null_ann[, skill_type := fcase(
  grepl("SC_Scaffolding",    coef), "Specialized socio-cognitive",
  grepl("SC_Specialized",    coef), "General socio-cognitive",
  grepl("Physical_Terminal", coef), "Physical-sensory"
)]
null_ann[, skill_type := factor(skill_type,
  levels = c("Specialized socio-cognitive",
             "General socio-cognitive",
             "Physical-sensory"))]
null_ann[, direction  := factor(direction,
  levels = c("\u03b2\u2191", "\u03b2\u2193"))]
null_ann[, flow_label := factor(
  fifelse(flow == "adoption", "Adoption", "Abandonment"),
  c("Adoption", "Abandonment"))]
null_ann[, ann_label := sprintf("%+.3f\n(%.0f SDs)", base_est, abs(sep_sds))]
null_ann[, hjust_val := fifelse(base_est >= 0, -0.08, 1.08)]

theme_si <- theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12,
                                    margin = margin(t = 2, b = 4)),
    axis.title       = element_text(size = 13),
    axis.title.x     = element_text(size = 12, margin = margin(t = 5)),
    axis.text        = element_text(size = 11, colour = "grey15"),
    panel.border     = element_rect(colour = "grey30", fill = NA,
                                    linewidth = 0.7),
    axis.line        = element_blank(),
    axis.ticks       = element_line(linewidth = 0.7),
    legend.position  = "none",
    panel.spacing.x  = unit(0.7, "lines"),
    panel.spacing.y  = unit(0.5, "lines"),
    plot.title       = element_text(face = "bold", size = 14,
                                    margin = margin(b = 4))
  )

# make_flow_plot(): 2 rows (direction) × 3 columns (skill type)
make_flow_plot <- function(flow_name) {
  pd  <- plot_dt[flow_label == flow_name]
  br  <- base_ref[flow_label == flow_name]
  ann <- null_ann[flow_label == flow_name]

  ggplot(pd, aes(x = estimate)) +
    geom_histogram(bins = 50, fill = "#4A90D9",
                   colour = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0,
               colour = "grey55", linewidth = 0.35,
               linetype = "dotted") +
    geom_vline(data = br,
               aes(xintercept = estimate),
               colour = "black", linewidth = 0.9) +
    geom_text(data = ann,
              aes(x = base_est, y = Inf,
                  label = ann_label,
                  hjust = hjust_val),
              vjust      = 1.2,
              size       = 3.8,
              colour     = "black",
              lineheight = 0.9,
              inherit.aes = FALSE) +
    # 2 rows (direction: β↑ / β↓) × 3 columns (skill type)
    facet_grid(direction ~ skill_type, scales = "free") +
    labs(title = flow_name,
         x     = "Permuted estimate (cloglog scale)",
         y     = "Count") +
    theme_si
}

p_adopt <- make_flow_plot("Adoption")
p_aband <- make_flow_plot("Abandonment")

fig_S7 <- (p_adopt / p_aband)

ggsave(file.path(SI_FIGS, "fig_SI_test_iii.pdf"),
       fig_S7, width = 13, height = 11,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(SI_FIGS, "fig_SI_test_iii.png"),
       fig_S7, width = 13, height = 11,
       units = "in", dpi = 300, bg = "white")

message("  Saved: fig_SI_test_iii.pdf / .png")
message("\n>>> SI_07_S7_test_iii_stratum_perm.R complete.")
message("    Next: SI_08_archetype_misclassification.R")