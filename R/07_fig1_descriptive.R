# ==============================================================================
# 07_fig1_descriptive.R
#
# Fig. 1 — Adoption and abandonment rates as a function of:
#   (A) Skill profile distance — adoption
#   (B) Skill profile distance — abandonment
#   (C) Status gap PC1        — adoption
#   (D) Status gap PC1        — abandonment
#
# STATUS GAP:
#   status_pc1 per occupation from occ_status_scores.csv
#   G_ij = status_target - status_source
#
# MEMORY STRATEGY:
#   1. Load adoption  → keep needed cols → bin → free RAM
#   2. Load abandonment → idem → bin → free RAM
#   3. Build figure from binned objects only
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
#         output/tables/main/pca_status_decision.csv
#
# Output: output/figures/main/fig1_descriptive.pdf / .png
# ==============================================================================

rm(list = ls()); gc(); gc()
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(igraph)
  library(ggraph)
  library(tibble)
  library(grid)
})

if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

out_figs <- file.path("output", "figures", "main")
dir.create(out_figs, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Rutas de Datos
# ==============================================================================
adopt_rds   <- "data/derived/riskset_adoption.rds"
aband_rds   <- "data/derived/riskset_abandonment.rds"
scores_path <- "output/tables/main/occ_status_scores.csv"
dec_path    <- "output/tables/main/pca_status_decision.csv"

for (f in c(adopt_rds, aband_rds, scores_path, dec_path)) stopifnot(file.exists(f))

pca_dec     <- fread(dec_path)
pc1_pct_var <- round(pca_dec$pc1_pct_var, 1)
message(sprintf("PC1: %.1f%% var | flip: %s", pc1_pct_var, ifelse(pca_dec$needs_flip, "YES", "NO")))

# ==============================================================================
# 2. Constantes y Estilos
# ==============================================================================
CLASS_LEVELS  <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")
CLASS_LABELS <- c("Specialized socio-cognitive", "General socio-cognitive", "sensory-physical")
CLASS_SHORT <- c("Spec. SC", "Gen. SC", "sensory-physical") 

CLASS_COLOURS <- c("SC_Scaffolding" = "#3B4992", "SC_Specialized" = "#008280", "Physical_Terminal" = "#EE0000")
CLASS_COLOURS_LIGHT <- c("SC_Scaffolding" = "#8E99D2", "SC_Specialized" = "#66B2B2", "Physical_Terminal" = "#F56666")
CLASS_FILLS   <- CLASS_COLOURS

PT_AX     <- 9.5; PT_LAB <- 10.0; PT_TITLE <- 11.5; PT_ANN <- 8.5
LWD       <- 0.85; ALPHA_RIB <- 0.12; BORDER_LW <- 0.75
N_BINS    <- 25; N_BOOT <- 300; SET_SEED <- 42; NODE_SIZE <- 1.8
FONT      <- "sans"

# ==============================================================================
# 3. Core Helpers & Network Generation
# ==============================================================================
load_scores_with_quintiles <- function() {
  sc <- fread(scores_path)
  setDT(sc)
  sc[, occ := as.character(occ)]
  qs <- quantile(sc$status_pc1, probs = seq(0, 1, 0.2), na.rm = TRUE)
  sc[, q_status := as.integer(cut(status_pc1, breaks = qs, include.lowest = TRUE))]
  sc
}

add_status_and_quintiles <- function(dt, scores) {
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]
  dt[scores, on = .(source = occ), `:=`(status_source = i.status_pc1, q_src = i.q_status)]
  dt[scores, on = .(target = occ), `:=`(status_target = i.status_pc1, q_dest = i.q_status)]
  dt[, pc1_gap := status_target - status_source]
  dt[, c("status_source", "status_target") := NULL]
  invisible(dt)
}

make_skill_class <- function(dt) {
  cs_med <- dt[domain == "Cognitive", median(cs, na.rm = TRUE)]
  dt[, skill_class := fcase(
    domain == "Cognitive" & cs >= cs_med, "SC_Scaffolding",
    domain == "Cognitive" & cs <  cs_med, "SC_Specialized",
    domain == "Physical",                  "Physical_Terminal"
  )]
  invisible(dt)
}

binned_ci <- function(dt, xvar, yvar) {
  d <- dt[!is.na(get(xvar)) & !is.na(get(yvar)) & skill_class %in% CLASS_LEVELS,
          .SD, .SDcols = c(xvar, yvar, "skill_class")]
  breaks <- quantile(d[[xvar]], probs = seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE)
  d[, xbin := cut(get(xvar), breaks = breaks, include.lowest = TRUE, labels = FALSE)]
  mids <- d[, .(xmid = median(get(xvar), na.rm = TRUE)), by = xbin]
  
  set.seed(SET_SEED)
  boot_res <- d[, {
    v <- get(yvar)
    obs <- mean(v, na.rm = TRUE)
    bts <- replicate(N_BOOT, mean(sample(v, length(v), replace = TRUE), na.rm = TRUE))
    list(ymean = obs, ylo = quantile(bts, 0.025), yhi = quantile(bts, 0.975))
  }, by = .(xbin, skill_class)]
  
  pd <- merge(boot_res, mids, by = "xbin")
  pd <- pd[xmid >= quantile(pd$xmid, 0.01, na.rm = TRUE) & xmid <= quantile(pd$xmid, 0.99, na.rm = TRUE)]
  pd[, class_f := factor(skill_class, levels = CLASS_LEVELS, labels = CLASS_LABELS)]
  pd
}

build_network_flows <- function(dt, outcome_col) {
  er <- dt[!is.na(q_src) & !is.na(q_dest) & q_src != q_dest & skill_class %in% CLASS_LEVELS,
           .(rate = mean(get(outcome_col), na.rm = TRUE), n = .N),
           by = .(skill_class, pc1_q_src = q_src, pc1_q_dest = q_dest)]
  er <- er[n >= 50] 
  er[, pc1_direction := fcase(pc1_q_dest > pc1_q_src, "Up", pc1_q_dest < pc1_q_src, "Down")]
  er[, rate_norm := rate / max(rate, na.rm = TRUE), by = skill_class]
  
  dr <- dt[!is.na(q_src) & !is.na(q_dest) & q_src != q_dest & skill_class %in% CLASS_LEVELS,
           .(Up = mean(get(outcome_col)[q_dest > q_src], na.rm = TRUE),
             Down = mean(get(outcome_col)[q_dest < q_src], na.rm = TRUE)),
           by = skill_class]
  dr_long <- melt(dr, id.vars = "skill_class", variable.name = "pc1_direction", value.name = "rate")
  
  list(edges = er, rates = dr_long)
}

pad_range <- function(r, lo = 0.04, hi = 0.12) c(r[1] - lo * diff(r), r[2] + hi * diff(r))

# ==============================================================================
# 4. Full Data Processing (computationally intensive)
# ==============================================================================
message("\n>>> Loading status scores and computing global quintiles...")
scores <- load_scores_with_quintiles()

message("\n>>> [1/2] Processing adoption (curves + networks)...")
da <- readRDS(adopt_rds)
setDT(da)
da <- da[, .(source, target, diffusion, domain, cs, structural_distance)]
gc()
make_skill_class(da)
add_status_and_quintiles(da, scores)
da <- da[!is.na(pc1_gap) & !is.na(structural_distance) & skill_class %in% CLASS_LEVELS]
pd_A <- binned_ci(da, "structural_distance", "diffusion")
pd_C <- binned_ci(da, "pc1_gap", "diffusion")
net_adopt <- build_network_flows(da, "diffusion")
rm(da); gc(); gc()

message("\n>>> [2/2] Processing abandonment (curves + networks)...")
db <- readRDS(aband_rds)
setDT(db)
db <- db[, .(source, target, abandonment, domain, cs, structural_distance)]
gc()
make_skill_class(db)
add_status_and_quintiles(db, scores)
db <- db[!is.na(pc1_gap) & !is.na(structural_distance) & skill_class %in% CLASS_LEVELS]
pd_B <- binned_ci(db, "structural_distance", "abandonment")
pd_D <- binned_ci(db, "pc1_gap", "abandonment")
net_aband <- build_network_flows(db, "abandonment")
rm(db); gc(); gc()

# ==============================================================================
# 5. Plot Generation (base panels, no X-axis titles)
# ==============================================================================
message("\n>>> Computing axis limits...")
XLIM_DIST <- pad_range(range(c(pd_A$xmid, pd_B$xmid), na.rm = TRUE), lo = 0.02, hi = 0.02)
XLIM_PC1  <- pad_range(range(c(pd_C$xmid, pd_D$xmid), na.rm = TRUE), lo = 0.02, hi = 0.02)
all_y     <- c(pd_A$ylo, pd_A$yhi, pd_B$ylo, pd_B$yhi, pd_C$ylo, pd_C$yhi, pd_D$ylo, pd_D$yhi)
YLIM_ALL  <- pad_range(range(all_y, na.rm = TRUE), lo = 0.04, hi = 0.14)

YLIM_CD    <- c(YLIM_ALL[1], YLIM_ALL[1] + 1.15 * (YLIM_ALL[2] - YLIM_ALL[1]))
YBREAKS_CD <- scales::pretty_breaks(n = 5)(YLIM_ALL)

theme_grad <- function() {
  theme_classic(base_size = PT_AX) +
    theme(
      plot.tag = element_text(face = "bold", size = PT_TITLE),
      axis.title = element_text(size = PT_LAB), 
      axis.title.x = element_blank(), # X-axis title suppressed here; added globally below
      axis.text = element_text(size = PT_AX), axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = BORDER_LW),
      legend.position = "none", plot.margin = margin(8, 10, 6, 8), aspect.ratio = 1
    )
}

make_grad <- function(pd, xlim, ylim, y_lab, tag, vline=FALSE, arrows=FALSE, ybreaks=NULL, legend=FALSE) {
  xs <- diff(xlim); ya <- ylim[2] * 1.04
  p <- ggplot(pd, aes(x = xmid, y = ymean, colour = class_f, fill = class_f)) +
    geom_ribbon(aes(ymin = ylo, ymax = yhi), alpha = ALPHA_RIB, colour = NA) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.45, linewidth = LWD, fullrange = TRUE,
                method.args = list(control = loess.control(surface = "direct"))) +
    geom_point(size = NODE_SIZE, shape = 16) +
    scale_colour_manual(NULL, values = setNames(CLASS_COLOURS, CLASS_LABELS)) +
    scale_fill_manual(NULL,   values = setNames(CLASS_COLOURS, CLASS_LABELS)) +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(.02, .02))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = ylim, expand = expansion(mult = c(0, 0)), breaks = if (!is.null(ybreaks)) ybreaks else waiver()) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = y_lab, tag = tag) + theme_grad() # x = NULL

  if (vline) p <- p + geom_vline(xintercept = 0, colour = "grey55", linewidth = .4, linetype = "dashed")
  if (arrows) {
    p <- p +
      annotate("text", x = xlim[1] + xs * .02, y = ya, label = "\u2190 Downward", size = PT_ANN/.pt, colour = "grey45", hjust = 0) +
      annotate("text", x = xlim[2] - xs * .02, y = ya, label = "Upward \u2192", size = PT_ANN/.pt, colour = "grey45", hjust = 1)
  }
  if (legend) {
    p <- p + theme(legend.position = "bottom") + guides(colour = guide_legend(NULL, nrow = 1, override.aes = list(size = 3.2, shape = 16)), fill = "none")
  }
  p
}

message(">>> Building panels...")
pA <- make_grad(pd_A, XLIM_DIST, YLIM_ALL, "Adoption rate", "(A)", ybreaks = YBREAKS_CD)
pB <- make_grad(pd_B, XLIM_DIST, YLIM_ALL, "Abandonment rate", "(B)", ybreaks = YBREAKS_CD)
pC_base <- make_grad(pd_C, XLIM_PC1, YLIM_CD, "Adoption rate", "(C)", vline = TRUE, arrows = TRUE, ybreaks = YBREAKS_CD)
pD_leg  <- make_grad(pd_D, XLIM_PC1, YLIM_CD, "Abandonment rate", "(D)", vline = TRUE, arrows = TRUE, ybreaks = YBREAKS_CD, legend = TRUE)

leg <- cowplot::get_legend(pD_leg)
pD_base <- pD_leg + theme(legend.position = "none")

# ==============================================================================
# 6. Network Insets: Build and Integrate
# ==============================================================================
message("\n>>> Building network flow insets...")

theme_mf <- function(col) {
  theme_void(base_family = FONT) + theme(
    plot.background = element_rect(fill = NA, colour = NA),
    plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5, colour = col, lineheight = 1.0, margin = margin(t = 1, b = 0)),
    plot.subtitle = element_text(size = 5.8, hjust = 0.5, colour = "grey45", margin = margin(b = 1)),
    plot.margin = margin(1, 3, 1, 3)
  )
}

verts5 <- tibble::tibble(name = paste0("Q", 1:5), x = 1:5, y = 0)
arrow_mf <- grid::arrow(angle = 20, length = unit(1.3, "mm"), type = "open", ends = "last")

make_mini_flow <- function(net_obj, arch) {
  col <- CLASS_COLOURS[[arch]]; col_dn <- CLASS_COLOURS_LIGHT[[arch]]
  lbl <- CLASS_SHORT[match(arch, CLASS_LEVELS)]
  
  er <- net_obj$edges[skill_class == arch]
  if (!nrow(er)) return(ggplot() + theme_void())
  
  up_r <- net_obj$rates[skill_class == arch & pc1_direction == "Up", rate]
  dn_r <- net_obj$rates[skill_class == arch & pc1_direction == "Down", rate]
  
  g <- igraph::graph_from_data_frame(
    d = data.frame(from = paste0("Q", er$pc1_q_src), to = paste0("Q", er$pc1_q_dest), 
                   weight = er$rate_norm, direction = er$pc1_direction, stringsAsFactors = FALSE),
    vertices = verts5, directed = TRUE
  )
  
  ggraph(g, layout = "linear") +
    geom_edge_arc(aes(edge_alpha = weight, edge_width = weight, edge_colour = direction), strength = 0.28, arrow = arrow_mf, end_cap = circle(1.6, "mm"), show.legend = FALSE) +
    scale_edge_width(range = c(0.06, 0.90), guide = "none") + scale_edge_alpha(range = c(0.20, 1.0), guide = "none") +
    scale_edge_colour_manual(values = c(Up = col, Down = col_dn)) +
    geom_node_point(size = 3.2, colour = col, fill = "white", shape = 21, stroke = 0.60) +
    geom_node_text(aes(label = name), vjust = 0.45, size = 1.75, fontface = "bold", colour = "grey15", family = FONT) +
    scale_x_continuous(expand = expansion(add = 0.55)) + coord_cartesian(clip = "off") +
    labs(title = lbl, subtitle = sprintf("\u2191%.1f%%  \u2193%.1f%%", up_r * 100, dn_r * 100)) + theme_mf(col)
}

strip_theme <- theme(plot.background = element_rect(fill = alpha("white", 0.88), colour = "grey78", linewidth = 0.35), plot.margin = margin(2, 4, 2, 4))
strip_C <- (make_mini_flow(net_adopt, "SC_Scaffolding") | make_mini_flow(net_adopt, "SC_Specialized") | make_mini_flow(net_adopt, "Physical_Terminal")) + plot_annotation(theme = strip_theme)
strip_D <- (make_mini_flow(net_aband, "SC_Scaffolding") | make_mini_flow(net_aband, "SC_Specialized") | make_mini_flow(net_aband, "Physical_Terminal")) + plot_annotation(theme = strip_theme)

pC <- pC_base + inset_element(strip_C, left = 0.0, bottom = 0.74, right = 1.0, top = 1.00, align_to = "panel", clip = FALSE)
pD <- pD_base + inset_element(strip_D, left = 0.0, bottom = 0.74, right = 1.0, top = 1.00, align_to = "panel", clip = FALSE)

# ==============================================================================
# 7. Assembly and Export (shared X-axis titles)
# ==============================================================================
message("\n>>> Creating shared X-axis labels and assembling final figure...")

# Textos unificados
xlab_top <- textGrob("Skill profile distance by occupations", 
                     gp = gpar(fontsize = PT_LAB, fontfamily = FONT))
xlab_bot <- textGrob(sprintf("Status gap between occupations\n(target \u2212 source) \u2014 PC1 (%.1f%% variance)", pc1_pct_var), 
                     gp = gpar(fontsize = PT_LAB, fontfamily = FONT, lineheight = 1.0))

# Vertical assembly
fig <- (pA | pB) / 
       patchwork::wrap_elements(full = xlab_top, clip = FALSE) / 
       (pC | pD) / 
       patchwork::wrap_elements(full = xlab_bot, clip = FALSE) / 
       patchwork::wrap_elements(full = leg) + 
       plot_layout(heights = c(1, 0.05, 1, 0.08, 0.10))

ggsave(file.path(out_figs, "fig1_full_networks.pdf"), fig, width = 7.2, height = 7.8, units = "in", device = cairo_pdf)
ggsave(file.path(out_figs, "fig1_full_networks.png"), fig, width = 7.2, height = 7.8, units = "in", dpi = 300)

message(sprintf("  Saved PDF: %s", file.path(out_figs, "fig1_full_networks.pdf")))
message(">>> Script complete.")