# ==============================================================================
# 08_fig2_regression.R
#
# Fig. 2 — Fitted log-hazard of skill adoption and abandonment as a function
# of signed status gap, under two complementary fixed-effect strategies.
#
# Layout:
#   Row 1 (Adoption):    Panel A Source FE | Panel B Target FE
#   Row 2 (Abandonment): Panel C Source FE | Panel D Target FE
#
# Input:  output/tables/main/coefs_pc1_adoption.csv
#         output/tables/abandonment/coefs_pc1_abandonment.csv
#         output/tables/main/pca_status_decision.csv
#
# Output: output/figures/main/fig2_regression.pdf / .png
# ==============================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
})

if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

out_figs   <- "output/figures/main"
out_tables <- "output/tables/main"
dir.create(out_figs,   showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Verify inputs
# ==============================================================================
inputs <- c(
  "output/tables/main/coefs_pc1_adoption.csv",
  "output/tables/abandonment/coefs_pc1_abandonment.csv",
  file.path(out_tables, "pca_status_decision.csv")
)
for (f in inputs) stopifnot(file.exists(f))
message("[OK] All inputs verified")

pca_dec     <- fread(file.path(out_tables, "pca_status_decision.csv"))
pc1_pct_var <- round(pca_dec$pc1_pct_var, 1)

# ==============================================================================
# Constants
# ==============================================================================
MAIN_ARCHETYPES <- c("SC_General", "SC_Specialized", "Physical_Terminal")
ARCH_LABELS <- c(
  SC_General    = "General socio-cognitive",
  SC_Specialized    = "Specialized socio-cognitive",
  Physical_Terminal = "Sensory-physical"
)
ARCH_COLOURS <- c(
  SC_General    = "#3B4992",
  SC_Specialized    = "#008280",
  Physical_Terminal = "#EE0000"
)

FONT       <- "Helvetica"
X_MAX      <-  5.0;  X_MIN  <- -5.0
N_PTS      <- 300
TICK_X     <- c(-5, -2.5, 0, 2.5, 5)
GLOBAL_Y_LO <- -2.0
GLOBAL_Y_HI <-  2.0
BORDER_LW  <- 1.2

# Rug geometry — adjusted for Y=[-2,+2]
Y_RUG_BASE <- -2.65
BAND_H     <-  0.11
TICK_H     <-  BAND_H * 0.75
y_band <- c(
  Physical_Terminal = Y_RUG_BASE + BAND_H * 0.5,
  SC_General    = Y_RUG_BASE + BAND_H * 1.5,
  SC_Specialized    = Y_RUG_BASE + BAND_H * 2.5
)

# Annotation geometry — adjusted for Y=[-2,+2]
X_RIGHT        <-  1.90
X_LEFT         <- -1.90
y_floor_labels <- -1.80
Y_SEP          <-  0.16

XLAB_SHARED <- sprintf(
  "Status gap between occupations (target \u2212 source)\nPC1 (%.1f%% variance explained)",
  pc1_pct_var)

# ==============================================================================
# Helpers
# ==============================================================================
extract_coefs <- function(path) {
  m  <- readRDS(path)
  ct <- as.data.table(m$coeftable, keep.rownames = "term")
  rm(m); gc()
  nms <- names(ct)
  if ("Estimate"   %in% nms) setnames(ct, "Estimate",   "estimate")
  if ("Std. Error" %in% nms) setnames(ct, "Std. Error", "std_error")
  if ("z value"    %in% nms) setnames(ct, "z value",    "z")
  if ("Pr(>|z|)"   %in% nms) setnames(ct, "Pr(>|z|)",  "p")
  ct[, archetype := fcase(
    grepl("SC_General",    term), "SC_General",
    grepl("SC_Specialized",    term), "SC_Specialized",
    grepl("Physical_Terminal", term), "Physical_Terminal",
    default = NA_character_
  )]
  ct[, var := fcase(
    grepl("pc1_dummy",       term) & !grepl("pc1_up|pc1_down", term), "kappa",
    grepl("pc1_up",          term),                                    "b_up",
    grepl("pc1_down",        term),                                    "b_dn",
    grepl("structural_dist", term),                                    "delta",
    default = "other"
  )]
  ct[!is.na(archetype) & var != "other",
     .(var, archetype,
       coef = round(estimate,  4),
       se   = round(std_error, 4),
       p    = round(p, 4))]
}

build_curves <- function(ct) {
  rbindlist(lapply(MAIN_ARCHETYPES, function(arch) {
    b_up  <- ct[var == "b_up"  & archetype == arch, coef[1]]
    b_dn  <- ct[var == "b_dn"  & archetype == arch, coef[1]]
    kappa <- ct[var == "kappa" & archetype == arch, coef[1]]
    if (!length(b_up) || is.na(b_up)) return(NULL)
    if (!length(kappa) || is.na(kappa)) kappa <- 0
    rbind(
      data.table(archetype = arch, segment = "up",
                 x  = seq(0,    X_MAX, length.out = N_PTS),
                 lp = kappa + b_up * seq(0, X_MAX, length.out = N_PTS)),
      data.table(archetype = arch, segment = "dn",
                 x  = seq(X_MIN, 0,   length.out = N_PTS),
                 lp = b_dn * seq(X_MIN, 0, length.out = N_PTS))
    )
  }))
}

build_ticks <- function(ct) {
  rbindlist(lapply(MAIN_ARCHETYPES, function(arch) {
    b_up  <- ct[var == "b_up"  & archetype == arch, coef[1]]
    b_dn  <- ct[var == "b_dn"  & archetype == arch, coef[1]]
    kappa <- ct[var == "kappa" & archetype == arch, coef[1]]
    if (!length(b_up) || is.na(b_up)) return(NULL)
    if (!length(kappa) || is.na(kappa)) kappa <- 0
    data.table(archetype = arch, x = TICK_X,
               lp = ifelse(TICK_X >= 0,
                           kappa + b_up * TICK_X,
                           b_dn  * TICK_X))
  }))
}

fmt_coef_se <- function(coef, se) {
  paste0(formatC(coef, format = "f", digits = 3, flag = "+"),
         " (", sprintf("%.3f", se), ")")
}

# ==============================================================================
# Theme
# ==============================================================================
theme_sa <- function() {
  theme_classic(base_size = 14, base_family = FONT) +
    theme(
      panel.border      = element_rect(colour = "grey20", fill = NA,
                                       linewidth = BORDER_LW),
      axis.line         = element_blank(),
      axis.title        = element_text(size = 13),
      axis.text         = element_text(size = 12, colour = "grey10"),
      legend.title      = element_blank(),
      legend.text       = element_text(size = 12),
      legend.position   = "bottom",
      legend.key.width  = unit(1.6, "cm"),
      panel.grid        = element_blank(),
      axis.ticks        = element_line(linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      plot.margin       = margin(8, 10, 38, 8, "pt"),
      aspect.ratio      = 1
    )
}

# ==============================================================================
# Panel builder
# ==============================================================================
build_panel <- function(ct, rug_segs, title, y_label,
                        show_legend = FALSE, panel_tag = NULL) {

  curves <- build_curves(ct)
  ticks  <- build_ticks(ct)
  curves[, archetype := factor(archetype, MAIN_ARCHETYPES, ARCH_LABELS)]
  ticks[,  archetype := factor(archetype, MAIN_ARCHETYPES, ARCH_LABELS)]

  # Annotations
  ann_up <- rbindlist(lapply(MAIN_ARCHETYPES, function(arch) {
    b  <- ct[var == "b_up" & archetype == arch, coef[1]]
    se <- ct[var == "b_up" & archetype == arch, se[1]]
    if (!length(b) || is.na(b)) return(NULL)
    data.table(archetype = arch, coef_val = b,
               label = paste0("\u03b2\u2191 = ", fmt_coef_se(b, se)))
  }))
  ann_dn <- rbindlist(lapply(MAIN_ARCHETYPES, function(arch) {
    b  <- ct[var == "b_dn" & archetype == arch, coef[1]]
    se <- ct[var == "b_dn" & archetype == arch, se[1]]
    if (!length(b) || is.na(b)) return(NULL)
    data.table(archetype = arch, coef_val = b,
               label = paste0("\u03b2\u2193 = ", fmt_coef_se(b, se)))
  }))

  ann_up <- ann_up[order(-coef_val)]
  ann_up[, y_pos := y_floor_labels + (.I - 1) * Y_SEP]
  ann_up[, x_pos := X_RIGHT]

  ann_dn <- ann_dn[order(-coef_val)]
  ann_dn[, y_pos := y_floor_labels + (.I - 1) * Y_SEP]
  ann_dn[, x_pos := X_LEFT]

  ann_all <- rbind(ann_up, ann_dn)
  ann_all[, archetype := factor(archetype, MAIN_ARCHETYPES, ARCH_LABELS)]

  # Rug — new format has columns: x, archetype
  rug_dt <- NULL
  if (!is.null(rug_segs)) {
    rug_dt <- copy(rug_segs)
    setDT(rug_dt)
    # Adapt format: may have column 'x' or 'status_pc1'
    if (!"x" %in% names(rug_dt) && "status_pc1" %in% names(rug_dt))
      setnames(rug_dt, "status_pc1", "x")
    rug_dt[, archetype_chr := as.character(archetype)]
    rug_dt[, y_ctr := y_band[archetype_chr]]
    rug_dt[, y     := y_ctr - TICK_H / 2]
    rug_dt[, yend  := y_ctr + TICK_H / 2]
    rug_dt[, xend  := x]
    rug_dt[, c("y_ctr", "archetype_chr") := NULL]
    rug_dt[, archetype := factor(archetype, MAIN_ARCHETYPES, ARCH_LABELS)]
    rug_dt <- rug_dt[!is.na(x) & is.finite(x) &
                       x >= X_MIN & x <= X_MAX]
  }

  y_breaks <- seq(GLOBAL_Y_LO, GLOBAL_Y_HI, by = 0.5)

  p <- ggplot() +
    annotate("rect",
             xmin = 0, xmax = X_MAX,
             ymin = GLOBAL_Y_LO, ymax = GLOBAL_Y_HI,
             fill = "grey97", alpha = 1) +
    geom_vline(xintercept = 0, colour = "grey70",
               linewidth = 0.35, linetype = "dotted") +
    geom_hline(yintercept = 0, colour = "grey70",
               linewidth = 0.35, linetype = "dotted")

  if (!is.null(rug_dt))
    p <- p + geom_segment(
      data = rug_dt,
      aes(x = x, xend = xend, y = y, yend = yend, colour = archetype),
      alpha = 0.18, linewidth = 0.25
    )

  p <- p +
    geom_line(data = curves[segment == "up"],
              aes(x = x, y = lp, colour = archetype),
              linewidth = 1.6, lineend = "round") +
    geom_line(data = curves[segment == "dn"],
              aes(x = x, y = lp, colour = archetype),
              linewidth = 1.6, lineend = "round") +
    geom_point(data = ticks,
               aes(x = x, y = lp, colour = archetype),
               size = 2.8, shape = 21, fill = "white", stroke = 1.0) +
    geom_text(data = ann_all,
              aes(x = x_pos, y = y_pos, label = label, colour = archetype),
              hjust = 0.5, size = 3.6, family = FONT,
              show.legend = FALSE) +
    scale_colour_manual(
      values = setNames(ARCH_COLOURS, ARCH_LABELS),
      labels = ARCH_LABELS, name = NULL) +
    scale_x_continuous(
      limits = c(X_MIN, X_MAX),
      breaks = c(-5, -2.5, 0, 2.5, 5),
      labels = c("-5", "-2.5", "0", "+2.5", "+5"),
      expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))) +
    coord_cartesian(
      xlim = c(X_MIN, X_MAX),
      ylim = c(GLOBAL_Y_LO, GLOBAL_Y_HI),
      clip = "off") +
    labs(title = title, x = NULL, y = y_label, tag = panel_tag) +
    theme_sa() +
    theme(
      plot.title = element_text(size = 13, face = "plain",
                                colour = "grey20", hjust = 0.5,
                                margin = margin(b = 4)),
      plot.tag   = element_text(face = "bold", size = 14)
    )

  if (!show_legend)
    p <- p + theme(legend.position = "none")
  else
    p <- p + guides(colour = guide_legend(
      nrow = 1,
      override.aes = list(size = 3.5, shape = 21, fill = "white")))
  p
}

# ==============================================================================
# Load coefficients from saved CSV tables (no model RDS required)
# ==============================================================================
message("Loading coefficients from CSV tables...")

load_coefs_csv <- function(path, panel_label) {
  dt <- fread(path)
  dt[panel == panel_label, .(var, archetype, coef, se, p)]
}

ct_adopt_A <- load_coefs_csv("output/tables/main/coefs_pc1_adoption.csv",        "Panel A")
ct_adopt_B <- load_coefs_csv("output/tables/main/coefs_pc1_adoption.csv",        "Panel B")
ct_aband_A <- load_coefs_csv("output/tables/abandonment/coefs_pc1_abandonment.csv", "Panel A")
ct_aband_B <- load_coefs_csv("output/tables/abandonment/coefs_pc1_abandonment.csv", "Panel B")

message("  Adoption Panel A: ", nrow(ct_adopt_A), " rows")
message("  Abandonment Panel A: ", nrow(ct_aband_A), " rows")

# Rug segments not available; plots render without marginal density strips
rug_adopt <- NULL
rug_aband <- NULL

# ==============================================================================
# Build 4 panels
# ==============================================================================
message("Building panels...")

pA <- build_panel(ct_adopt_A, rug_adopt,
                  title     = "Panel A \u2014 Source FE",
                  y_label   = "Log-hazard of skill adoption",
                  panel_tag = "(A)")

pB <- build_panel(ct_adopt_B, rug_adopt,
                  title     = "Panel B \u2014 Target FE",
                  y_label   = NULL,
                  panel_tag = "(B)")

pC <- build_panel(ct_aband_A, rug_aband,
                  title     = "Panel A \u2014 Source FE",
                  y_label   = "Log-hazard of skill abandonment",
                  panel_tag = "(C)")

pD <- build_panel(ct_aband_B, rug_aband,
                  title     = "Panel B \u2014 Target FE",
                  y_label   = NULL,
                  panel_tag = "(D)")

# ==============================================================================
# Ensamble
# ==============================================================================
make_row_label <- function(txt) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             angle = 90, size = 4.5, fontface = "bold",
             colour = "grey30", family = FONT) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

xlab_grob <- grid::textGrob(
  XLAB_SHARED, y = unit(1, "npc"), vjust = 1,
  gp = grid::gpar(fontsize = 13, fontfamily = FONT, lineheight = 1.1)
)

leg <- cowplot::get_legend(
  pA + theme(legend.position = "bottom") +
    guides(colour = guide_legend(
      nrow = 1,
      override.aes = list(size = 3.5, shape = 21, fill = "white"))))

row_adopt <- (make_row_label("Adoption")    | pA | pB) +
  plot_layout(widths = c(0.04, 1, 1))
row_aband <- (make_row_label("Abandonment") | pC | pD) +
  plot_layout(widths = c(0.04, 1, 1))

fig <- row_adopt / row_aband /
  patchwork::wrap_elements(full = xlab_grob, clip = FALSE) /
  patchwork::wrap_elements(full = leg) +
  plot_layout(heights = c(1, 1, 0.06, 0.07))

# ==============================================================================
# Save
# ==============================================================================
message("Saving figure...")

ggsave(file.path(out_figs, "fig2_regression.pdf"),
       fig, width = 11.0, height = 12.0,
       units = "in", device = cairo_pdf)
ggsave(file.path(out_figs, "fig2_regression.png"),
       fig, width = 11.0, height = 12.0,
       units = "in", dpi = 300)

message("  Saved: output/figures/main/fig2_regression.pdf")
message("  Saved: output/figures/main/fig2_regression.png")
message("\n>>> 08_fig2_regression.R complete.")
message("    Next: 09_fig3_projections.R (generated by 06_projections.R)")