# ==============================================================================
# SI_09_alternative_status.R
#
# Section S3.3.1: Alternative Status Measures
#
# Standardizes each component of the status index at the occupation level
# (z-score, SD=1) before computing dyadic gaps, so that coefficients are
# comparable across measures and with the main model
# (PC1, which is standardized by construction).
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
# Output: output/figures/SI/fig_SI_alt_status.pdf / .png
#         output/tables/SI/tab_SI_alt_status.csv
# ==============================================================================

library(data.table)
library(fixest)
library(ggplot2)
library(ggsci)
library(patchwork)

source("R/SI/00_setup_SI.R")

dir.create("output/figures/SI", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables/SI",  recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load and standardize occupation scores
# ==============================================================================
message(">>> Step 1: Loading and standardizing occupation scores...")

occ <- fread("output/tables/main/occ_status_scores.csv")
occ[, occ := as.character(occ)]

# Standardize each component at the occupation level (z-score)
occ[, log_wage_z := as.numeric(scale(log_wage))]
occ[, log_edu_z  := as.numeric(scale(log_edu))]
occ[, cog_z      := as.numeric(scale(cog))]

message(sprintf("  Occupations: %d", nrow(occ)))
message(sprintf("  log_wage_z: mean=%.4f sd=%.4f", mean(occ$log_wage_z), sd(occ$log_wage_z)))
message(sprintf("  log_edu_z:  mean=%.4f sd=%.4f", mean(occ$log_edu_z),  sd(occ$log_edu_z)))
message(sprintf("  cog_z:      mean=%.4f sd=%.4f", mean(occ$cog_z),      sd(occ$cog_z)))

# ==============================================================================
# 2. Function to enrich a risk set with standardized gaps
# ==============================================================================
enrich_standardized <- function(dt, occ_scores) {
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]

  # Merge standardized scores for source and target
  for (var in c("log_wage_z", "log_edu_z", "cog_z")) {
    s_var <- paste0("s_", var)
    t_var <- paste0("t_", var)
    dt[occ_scores, on = .(source = occ), (s_var) := get(paste0("i.", var))]
    dt[occ_scores, on = .(target = occ), (t_var) := get(paste0("i.", var))]
  }

  # Compute standardized dyadic gaps
  for (var in c("log_wage_z", "log_edu_z", "cog_z")) {
    gap_var  <- paste0(var, "_gap")
    up_var   <- paste0(var, "_up")
    down_var <- paste0(var, "_down")
    dummy_var <- paste0(var, "_dummy")
    s_var    <- paste0("s_", var)
    t_var    <- paste0("t_", var)
    dt[, (gap_var)   := get(t_var) - get(s_var)]
    dt[, (up_var)    := pmax(0, get(gap_var))]
    dt[, (down_var)  := pmin(0, get(gap_var))]
    dt[, (dummy_var) := fifelse(!is.na(get(gap_var)) & get(gap_var) > 0, 1L, 0L)]
  }
  dt
}

# ==============================================================================
# 3. Setup modelos
# ==============================================================================
flow_files <- c(
  Adoption    = "data/derived/riskset_adoption.rds",
  Abandonment = "data/derived/riskset_abandonment.rds"
)
dep_vars <- c(Adoption = "diffusion", Abandonment = "abandonment")

measure_cols <- list(
  `Log wage`        = c(up = "log_wage_z_up",  down = "log_wage_z_down",  dummy = "log_wage_z_dummy"),
  `Log education`   = c(up = "log_edu_z_up",   down = "log_edu_z_down",   dummy = "log_edu_z_dummy"),
  `Cognitive score` = c(up = "cog_z_up",       down = "cog_z_down",       dummy = "cog_z_dummy")
)

form_rhs_template <- paste0(
  "~ 0 + i(atc_archetype, %s) + i(atc_archetype, %s) + ",
  "i(atc_archetype, %s) + i(atc_archetype, structural_distance) | source + skill_name"
)

KEEP_BASE <- c("source", "target", "skill_name", "atc_archetype", "structural_distance")

# ==============================================================================
# 4. Estimation
# ==============================================================================
models <- list()

for (flow_name in names(flow_files)) {
  cat(sprintf("\n>>> Loading %s risk set...\n", flow_name))
  y_var        <- dep_vars[[flow_name]]
  dt_flow_full <- readRDS(flow_files[[flow_name]])
  setDT(dt_flow_full)
  dt_flow_full[, atc_archetype := as.factor(atc_archetype)]

  # Enriquecer con gaps estandarizados
  cat("  Enriching with standardized gaps...\n")
  dt_flow_full <- enrich_standardized(dt_flow_full, occ)
  gc()

  for (measure in names(measure_cols)) {
    cat(sprintf("  Estimating %s — %s...\n", flow_name, measure))
    cols   <- measure_cols[[measure]]
    keep   <- c(KEEP_BASE, y_var, unname(cols))
    dt_sub <- dt_flow_full[, ..keep]

    form_rhs <- sprintf(form_rhs_template,
                        cols[["up"]], cols[["down"]], cols[["dummy"]])

    mod <- feglm(
      as.formula(paste0(y_var, form_rhs)),
      data      = dt_sub,
      family    = binomial("cloglog"),
      cluster   = ~ source + target + skill_name,
      mem.clean = TRUE,
      nthreads  = 0,
      lean      = TRUE
    )

    models[[paste(flow_name, measure, sep = "__")]] <- mod
    rm(dt_sub); gc()
  }
  rm(dt_flow_full); gc()
}

# ==============================================================================
# 5. Extraer coeficientes
# ==============================================================================
extract_coefs <- function(mod, flow, measure) {
  ct <- as.data.table(coeftable(mod), keep.rownames = "term")
  setnames(ct, c("term", "coef", "se", "t", "p"))
  ct <- ct[grepl("_up|_down", term) & !grepl("dummy|structural", term)]
  ct[, direction := fcase(
    grepl("_up",   term), "Upward gap slope",
    grepl("_down", term), "Downward gap slope"
  )]
  ct[, archetype := fcase(
    grepl("SC_Scaffolding",    term), "General SC",
    grepl("SC_Specialized",    term), "Specialized SC",
    grepl("Physical_Terminal", term), "Physical-sensory"
  )]
  ct[, flow    := flow]
  ct[, measure := measure]
  ct[, ci_lo   := coef - 1.96 * se]
  ct[, ci_hi   := coef + 1.96 * se]
  ct[, .(flow, measure, archetype, direction, coef, se, ci_lo, ci_hi, p)]
}

coef_dt <- rbindlist(mapply(
  function(mod, nm) {
    parts <- strsplit(nm, "__")[[1]]
    extract_coefs(mod, parts[1], parts[2])
  },
  models, names(models), SIMPLIFY = FALSE
))

# Factor ordering
coef_dt[, archetype := factor(archetype,
  levels = c("Specialized SC", "General SC", "Physical-sensory"))]
coef_dt[, measure := factor(measure,
  levels = c("Log wage", "Log education", "Cognitive score"))]
coef_dt[, direction := factor(direction,
  levels = c("Upward gap slope", "Downward gap slope"))]
coef_dt[, flow := factor(flow, levels = c("Adoption", "Abandonment"))]

fwrite(coef_dt, "output/tables/SI/tab_SI_alt_status.csv")
message("  Saved: output/tables/SI/tab_SI_alt_status.csv")

# ==============================================================================
# 6. Figura
# ==============================================================================
pal3 <- pal_aaas("default")(3)
names(pal3) <- c("Log wage", "Log education", "Cognitive score")

theme_sa <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, face = "bold"),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9, colour = "grey10"),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(linewidth = 0.35, colour = "grey20"),
    axis.ticks         = element_line(linewidth = 0.3),
    legend.position    = "bottom",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 9)
  )

pd <- position_dodge(width = 0.55)

p <- ggplot(coef_dt,
            aes(x = archetype, y = coef,
                colour = measure, shape = measure,
                ymin = ci_lo, ymax = ci_hi)) +
  geom_hline(yintercept = 0, linewidth = 0.35,
             colour = "grey40", linetype = "dashed") +
  geom_errorbar(position = pd, width = 0.25, linewidth = 0.5) +
  geom_point(position = pd, size = 2.2) +
  facet_grid(direction ~ flow, scales = "free_y") +
  scale_colour_manual(name = "Status measure", values = pal3) +
  scale_shape_manual(name  = "Status measure", values = c(16, 17, 15)) +
  scale_x_discrete(labels = c(
    "Specialized SC"   = "Spec. SC",
    "General SC"       = "Gen. SC",
    "Physical-sensory" = "Physical"
  )) +
  labs(
    x = NULL,
    y = "Estimate (cloglog scale, per SD of status gap)",
    #subtitle = paste0(
    #  "All status measures standardized to SD = 1 before computing directional gaps. ",
    #  "Source and skill fixed effects (Panel A). ",
    #  "Bars are 95% CIs; standard errors clustered three-way."
    #)
  ) +
  theme_sa

ggsave("output/figures/SI/fig_SI_alt_status.pdf",
       plot = p, width = 9, height = 7,
       units = "in", device = cairo_pdf, bg = "white")
ggsave("output/figures/SI/fig_SI_alt_status.png",
       plot = p, width = 9, height = 7,
       units = "in", dpi = 300, bg = "white")
message("  Saved: output/figures/SI/fig_SI_alt_status.pdf / .png")
message("\n>>> SI_11_alternative_status.R complete.")
