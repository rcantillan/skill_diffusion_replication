# ==============================================================================
# SI_09_alternative_status.R
#
# Section S3.3.1: Alternative Status Measures
#
# Estima modelos cloglog con tres medidas alternativas de estatus (Log Wage,
# Log Educ, Cognitive), estandarizadas en z-score al nivel de ocupación, para
# adoption y abandonment (Panel A: Source+Skill FE, Panel B: Target+Skill FE).
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
# Output: output/figures/si/fig_SI_alt_status_Adoption.pdf/.png
#         output/figures/si/fig_SI_alt_status_Abandonment.pdf/.png
#         output/tables/si/tab_SI_alt_status.csv
# ==============================================================================

library(data.table)
library(fixest)
library(ggplot2)
library(patchwork)

source("R/SI/00_setup_SI.R")

dir.create("output/figures/si", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables/si",  recursive = TRUE, showWarnings = FALSE)

coef_csv <- "output/tables/si/tab_SI_alt_status.csv"
if (file.exists(coef_csv)) {
  file.remove(coef_csv)
  message(">>> [cache] removed: ", coef_csv)
}

# ==============================================================================
# 1. Load & standardise occupation scores
# ==============================================================================
message(">>> Step 1: occupation scores...")
occ <- fread("output/tables/main/occ_status_scores.csv")
occ[, occ       := as.character(occ)]
occ[, log_wage_z := as.numeric(scale(log_wage))]
occ[, log_edu_z  := as.numeric(scale(log_edu))]
occ[, cog_z      := as.numeric(scale(cog))]
message(sprintf("  n_occ=%d | log_wage_z sd=%.4f | log_edu_z sd=%.4f | cog_z sd=%.4f",
                nrow(occ), sd(occ$log_wage_z), sd(occ$log_edu_z), sd(occ$cog_z)))

# ==============================================================================
# 2. Enrich risk set with standardised dyadic gaps
# ==============================================================================
enrich_standardized <- function(dt, occ_scores) {
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  for (var in c("log_wage_z", "log_edu_z", "cog_z")) {
    dt[occ_scores, on = .(source = occ), (paste0("s_", var)) := get(paste0("i.", var))]
    dt[occ_scores, on = .(target = occ), (paste0("t_", var)) := get(paste0("i.", var))]
    dt[, (paste0(var, "_gap"))   := get(paste0("t_", var)) - get(paste0("s_", var))]
    dt[, (paste0(var, "_up"))    := pmax(0, get(paste0(var, "_gap")))]
    dt[, (paste0(var, "_down"))  := pmin(0, get(paste0(var, "_gap")))]
    dt[, (paste0(var, "_dummy")) := fifelse(
      !is.na(get(paste0(var, "_gap"))) & get(paste0(var, "_gap")) > 0, 1L, 0L)]
  }
  dt
}

# ==============================================================================
# 3. Setup
# ==============================================================================
flow_files <- c(Adoption    = "data/derived/riskset_adoption.rds",
                Abandonment = "data/derived/riskset_abandonment.rds")
dep_vars   <- c(Adoption = "diffusion", Abandonment = "abandonment")

measure_cols <- list(
  `Log Wage`  = c(up = "log_wage_z_up",  down = "log_wage_z_down",  dummy = "log_wage_z_dummy"),
  `Log Educ`  = c(up = "log_edu_z_up",   down = "log_edu_z_down",   dummy = "log_edu_z_dummy"),
  `Cognitive` = c(up = "cog_z_up",       down = "cog_z_down",       dummy = "cog_z_dummy")
)

fe_specs <- list(
  "Panel A \u2014 Source FE" = "| source + skill_name",
  "Panel B \u2014 Target FE" = "| target + skill_name"
)

KEEP_BASE <- c("source", "target", "skill_name", "atc_archetype", "structural_distance")

form_rhs_template <- paste0(
  "~ 0 + i(atc_archetype, %s) + i(atc_archetype, %s) + ",
  "i(atc_archetype, %s) + i(atc_archetype, structural_distance)"
)

# ==============================================================================
# 4. Extract coefficients
# ==============================================================================
extract_coefs <- function(mod, flow, measure, panel) {
  ct <- as.data.table(coeftable(mod), keep.rownames = "term")
  setnames(ct, c("term", "coef", "se", "t", "p"))


  ct <- ct[grepl("_up|_down", term) & !grepl("dummy|structural_distance", term)]

  ct[, direction := fcase(
    grepl("_up",   term), "up",
    grepl("_down", term), "down",
    default = NA_character_
  )]
  ct[, archetype := fcase(
    grepl("SC_Specialized",    term, fixed = TRUE), "Specialized socio-cognitive",
    grepl("SC_General",        term, fixed = TRUE), "General socio-cognitive",
    grepl("Physical_Terminal", term, fixed = TRUE), "Sensory-physical",
    default = NA_character_
  )]

  if (any(is.na(ct$archetype))) {
    warning(sprintf("[%s|%s|%s] unmatched archetype: %s",
                    flow, measure, panel,
                    paste(ct[is.na(archetype), term], collapse = "; ")))
  }

  ct <- ct[!is.na(archetype) & !is.na(direction)]
  ct[, flow    := flow]
  ct[, measure := measure]
  ct[, panel   := panel]
  ct[, ci_lo   := coef - 1.96 * se]
  ct[, ci_hi   := coef + 1.96 * se]
  ct[, .(flow, panel, measure, archetype, direction, coef, se, ci_lo, ci_hi, p)]
}

# ==============================================================================
# 5. Estimation loop
# ==============================================================================
all_coefs <- list()

for (flow_name in names(flow_files)) {
  message(sprintf("\n>>> %s risk set", flow_name))
  y_var        <- dep_vars[[flow_name]]
  dt_flow_full <- readRDS(flow_files[[flow_name]])
  setDT(dt_flow_full)
  dt_flow_full[, atc_archetype := as.factor(atc_archetype)]


  message("  enriching...")
  dt_flow_full <- enrich_standardized(dt_flow_full, occ)

  new_cols <- c("log_wage_z_up", "log_wage_z_down",
                "log_edu_z_up",  "log_edu_z_down",
                "cog_z_up",      "cog_z_down")
  missing <- setdiff(new_cols, names(dt_flow_full))
  if (length(missing)) stop("Missing columns after enrich: ", paste(missing, collapse = ", "))
  message("  [enriched columns OK]")
  gc()

  for (measure in names(measure_cols)) {
    cols   <- measure_cols[[measure]]
    keep   <- c(KEEP_BASE, y_var, unname(cols))
    dt_sub <- dt_flow_full[, ..keep]
    message(sprintf("  [%s] nrow=%d", measure, nrow(dt_sub)))

    for (panel_name in names(fe_specs)) {
      message(sprintf("  Estimating: %s | %s | %s", flow_name, measure, panel_name))

      form_rhs <- paste0(
        sprintf(form_rhs_template, cols[["up"]], cols[["down"]], cols[["dummy"]]),
        fe_specs[[panel_name]]
      )

      mod <- feglm(
        as.formula(paste0(y_var, form_rhs)),
        data      = dt_sub,
        family    = binomial("cloglog"),
        cluster   = ~ source + target + skill_name,
        mem.clean = TRUE,
        nthreads  = 0,
        lean      = TRUE
      )

      key <- paste(flow_name, measure, panel_name, sep = "__")
      res <- extract_coefs(mod, flow_name, measure, panel_name)
      all_coefs[[key]] <- res
      rm(mod); gc()
    }
    rm(dt_sub); gc()
  }
  rm(dt_flow_full); gc()
}

message(sprintf("\n>>> Combining %d result sets...", length(all_coefs)))
coef_dt <- rbindlist(all_coefs)
message(sprintf("  Total rows: %d (expected 72)", nrow(coef_dt)))
fwrite(coef_dt, coef_csv)
message("  Saved: ", coef_csv)

# ==============================================================================
# 6. Visualisation
# ==============================================================================

# ── Niveles canónicos ─────────────────────────────────────────────────────────
arch_levels    <- c("General socio-cognitive",
                    "Specialized socio-cognitive",
                    "Sensory-physical")
measure_levels <- c("Log Wage", "Log Educ", "Cognitive")
panel_levels   <- c("Panel A \u2014 Source FE", "Panel B \u2014 Target FE")

coef_dt[, archetype := factor(archetype, levels = arch_levels)]
coef_dt[, measure   := factor(measure,   levels = measure_levels)]
coef_dt[, panel     := factor(panel,     levels = panel_levels)]
coef_dt[, flow      := factor(flow,      levels = c("Adoption", "Abandonment"))]

pal_arch <- c(
  "General socio-cognitive"     = "#3B4992",
  "Specialized socio-cognitive" = "#008280",
  "Sensory-physical"            = "#EE0000"
)

# ── Wide ──────────────────────────────────────────────────────────────────────
coef_wide <- dcast(coef_dt,
                   flow + panel + measure + archetype ~ direction,
                   value.var = c("coef", "se"))

# ── Slope data: V-shape broken at 0 ──────────────────────────────────────────
x_vals   <- c(-5, -2.5, 0, 2.5, 5)
slope_dt <- coef_wide[rep(seq_len(.N), each = length(x_vals))]
slope_dt[, x := rep(x_vals, times = nrow(coef_wide))]
slope_dt[, y := fifelse(x <= 0, coef_down * x, coef_up * x)]
slope_dt[, archetype := factor(archetype, levels = arch_levels)]
slope_dt[, measure   := factor(measure,   levels = measure_levels)]
slope_dt[, panel     := factor(panel,     levels = panel_levels)]

# ── Anotaciones ───────────────────────────────────────────────────────────────
# Stack: Physical (rojo) arriba, General (azul) medio, Specialized (teal) abajo
# Posición Y absoluta fija, proporcional al rango del eje de cada medida:
#   Log Wage / Log Educ → eje ±2.5, sep = 0.32
#   Cognitive           → eje ±4.0, sep = 0.51 (= 0.32 × 8/5)
# β↓: hjust=1, x=-0.20  |  β↑: hjust=0, x=+0.20

stack_order <- c("Sensory-physical",
                 "General socio-cognitive",
                 "Specialized socio-cognitive")

y_top_by_measure <- c("Log Wage" = -1.70, "Log Educ" = -1.70, "Cognitive" = -2.90)
y_sep_by_measure <- c("Log Wage" =  0.32, "Log Educ" =  0.32, "Cognitive" =  0.51)

make_annot <- function(flow_sel) {
  cw <- coef_wide[flow == flow_sel]
  rows <- list()

  for (pn in panel_levels) {
    for (ms in measure_levels) {
      for (rank_i in seq_along(stack_order)) {
        arch <- stack_order[[rank_i]]
        row  <- cw[panel == pn & measure == ms & archetype == arch]
        if (nrow(row) == 0) next

        y_pos <- y_top_by_measure[[ms]] - (rank_i - 1) * y_sep_by_measure[[ms]]

        for (dir in c("down", "up")) {
          sym    <- if (dir == "down") "\u03b2\u2193" else "\u03b2\u2191"
          coef_  <- if (dir == "down") row$coef_down  else row$coef_up
          se_    <- if (dir == "down") row$se_down     else row$se_up
          x_pos  <- if (dir == "down") -0.20 else  0.20
          hjust_ <- if (dir == "down")  1    else  0

          rows[[length(rows) + 1]] <- data.table(
            flow      = flow_sel,
            panel     = pn,
            measure   = ms,
            archetype = arch,
            direction = dir,
            label     = sprintf("%s = %+.3f (%.3f)", sym, coef_, se_),
            y_pos     = y_pos,
            x_pos     = x_pos,
            hjust_val = hjust_
          )
        }
      }
    }
  }
  ann <- rbindlist(rows)
  ann[, archetype := factor(archetype, levels = arch_levels)]
  ann[, measure   := factor(measure,   levels = measure_levels)]
  ann[, panel     := factor(panel,     levels = panel_levels)]
  ann
}

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_sa <- theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    panel.border      = element_rect(colour = "grey10", linewidth = 1.4, fill = NA),
    panel.background  = element_rect(fill = "white"),
    panel.grid.major  = element_line(linewidth = 0.2, colour = "grey88"),
    panel.grid.minor  = element_blank(),
    strip.background  = element_blank(),
    strip.text.x      = element_text(size = 10, face = "bold"),
    strip.text.y      = element_text(size = 10, face = "bold", angle = 90),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9, colour = "grey10"),
    axis.ticks        = element_line(linewidth = 0.3),
    plot.title        = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 9),
    legend.key.width  = unit(2.2, "lines"),
    legend.spacing.x  = unit(0.5, "cm"),
    legend.background = element_blank()
  )

# ── Build: 3 filas por patchwork (eje Y propio por medida) ───────────────────
build_row <- function(flow_sel, measure_sel,
                      show_xaxis  = FALSE,
                      show_legend = FALSE) {

  ylim_val <- if (measure_sel == "Cognitive") c(-4, 4) else c(-2.5, 2.5)
  ybreaks  <- if (measure_sel == "Cognitive") c(-4, -2, 0, 2, 4) else
                                               c(-2.5, -1.25, 0, 1.25, 2.5)

  sd  <- slope_dt[flow == flow_sel & measure == measure_sel]
  ann <- make_annot(flow_sel)[measure == measure_sel]

  p <- ggplot(sd, aes(x = x, y = y, colour = archetype, group = archetype)) +
    geom_vline(xintercept = 0, linewidth = 0.4,
               colour = "grey50", linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.3,
               colour = "grey50", linetype = "dotted") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.6, shape = 21, fill = "white", stroke = 1.3) +
    geom_text(
      show.legend = FALSE,
      data        = ann,
      aes(x = x_pos, y = y_pos, label = label,
          colour = archetype, hjust = hjust_val),
      inherit.aes = FALSE,
      size        = 2.7,
      family      = "Helvetica"
    ) +
    facet_grid(measure ~ panel, switch = "y") +
    scale_colour_manual(values = pal_arch, breaks = arch_levels) +
    scale_x_continuous(
      breaks = c(-5, -2.5, 0, 2.5, 5),
      labels = c("-5", "-2.5", "0", "+2.5", "+5"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(breaks = ybreaks) +
    coord_cartesian(ylim = ylim_val) +
    labs(x = if (show_xaxis) "Status gap (Z-score)" else NULL,
         y = "Log-hazard") +
    theme_sa

  if (!show_xaxis) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  }
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + guides(colour = guide_legend(
      nrow         = 1,
      override.aes = list(shape = 21, size = 3.5,
                          fill = "white", stroke = 1.4)
    ))
  }
  p
}

build_fig <- function(flow_sel) {
  p1 <- build_row(flow_sel, "Log Wage",  show_xaxis = FALSE, show_legend = FALSE)
  p2 <- build_row(flow_sel, "Log Educ",  show_xaxis = FALSE, show_legend = FALSE)
  p3 <- build_row(flow_sel, "Cognitive", show_xaxis = TRUE,  show_legend = TRUE)
  p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
}

# ==============================================================================
# 7. Save
# ==============================================================================
for (flow_sel in c("Adoption", "Abandonment")) {
  p     <- build_fig(flow_sel)
  fname <- paste0("output/figures/si/fig_SI_alt_status_", flow_sel)
  ggsave(paste0(fname, ".pdf"), plot = p, width = 10, height = 11,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(paste0(fname, ".png"), plot = p, width = 10, height = 11,
         units = "in", dpi = 300, bg = "white")
  message("Saved: ", fname, ".pdf / .png")
}

message("\n>>> SI_09_alternative_status.R complete.")