# ==============================================================================
# SI_05_rca_threshold.R
# Section S3.3: RCA Threshold Sensitivity
#
# Produces Fig. S1: RCA threshold sensitivity across four specialization
# thresholds (0.90, 1.00, 1.10, 1.25), for both adoption and abandonment
# flows.
#
# Layout: 2 rows (Adoption | Abandonment) x 4 columns (thresholds)
#         Free Y axis per row to allow each flow's natural scale.
#
# The figure documents that the directional asymmetry is not sensitive to
# the choice of specialization threshold: beta_up_Cog > 0 and
# beta_up_Phy < 0 at all four thresholds in adoption; the reverse holds
# in abandonment. This confirms the directional asymmetry is not an
# artifact of the RCA = 1 cutpoint.
#
# Design: binned descriptive approach (observed adoption/abandonment rates
# by status gap decile). Full model estimates at alternative thresholds
# are reported in the coefficient table.
#
# Input:  data/derived/riskset_adoption.rds
#         data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
#         data/raw/onet/db_15_1/
#         data/raw/onet/db_29_2_text/
#         data/derived/skill_cs_scores.rds
#         data/crosswalk/2010_to_2019_Crosswalk.csv
#
# Output: output/figures/si/fig_S1_rca_threshold.pdf / .png
#         output/tables/si/table_rca_threshold_adoption_rates.csv
#         output/tables/si/table_rca_threshold_abandonment_rates.csv
# ==============================================================================

library(data.table)
library(ggplot2)
library(ggsci)
library(scales)

# ==============================================================================
# Constants
# ==============================================================================
THRESHOLDS <- c(0.90, 1.00, 1.10, 1.25)
N_BINS     <- 10L
N_BOOT     <- 500L
SEED_FIG   <- 42L

CLASS_LEVELS  <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")
CLASS_LABELS  <- c("Specialized socio-cognitive",
                   "General socio-cognitive",
                   "Physical-sensory")
CLASS_COLOURS <- c(
  "Specialized socio-cognitive" = "#3B4992",
  "General socio-cognitive"     = "#008280",
  "Physical-sensory"            = "#EE0000"
)

SI_TABLES <- "output/tables/si"
SI_FIGS   <- "output/figures/si"
dir.create(SI_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(SI_FIGS,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 1 — Load status scores and archetype lookup
# ==============================================================================
message(">>> Step 1: Loading status scores and base risk sets...")

scores <- fread("output/tables/main/occ_status_scores.csv")
setDT(scores)
scores[, occ := as.character(occ)]
scores_lkp <- scores[, .(occ, s_pc1 = status_pc1)]

# Archetype lookup from baseline adoption risk set
dt_base <- readRDS("data/derived/riskset_adoption.rds")
setDT(dt_base)
archetype_lkp <- unique(dt_base[, .(skill_name, atc_archetype)])
rm(dt_base); gc()

# Status gap decile breaks — fixed across thresholds and flows
# for visual comparability
dt_tmp <- readRDS("data/derived/riskset_adoption.rds")
setDT(dt_tmp)
dt_tmp[scores_lkp, on = .(source = occ), s_pc1 := i.s_pc1]
dt_tmp[scores_lkp, on = .(target = occ), t_pc1 := i.s_pc1]
dt_tmp <- dt_tmp[!is.na(s_pc1) & !is.na(t_pc1)]
gap_breaks <- unique(quantile(
  dt_tmp$t_pc1 - dt_tmp$s_pc1,
  probs = seq(0, 1, length.out = N_BINS + 1L),
  na.rm = TRUE
))
rm(dt_tmp); gc()
message(sprintf("  Gap breaks: [%.3f, %.3f] in %d bins",
                min(gap_breaks), max(gap_breaks), length(gap_breaks) - 1L))

# ==============================================================================
# Step 2 — Load O*NET importance data
# ==============================================================================
message(">>> Step 2: Loading O*NET importance data...")

read_onet_im <- function(folder) {
  files <- c("Skills.txt", "Abilities.txt", "Knowledge.txt")
  rbindlist(lapply(files, function(f) {
    path <- file.path(folder, f)
    if (!file.exists(path)) return(NULL)
    d <- fread(path, sep = "\t", quote = "",
               na.strings = c("NA", "n/a", "", "*"),
               showProgress = FALSE)
    setnames(d, names(d), make.names(names(d)))
    d <- d[Scale.ID == "IM"]
    clean_id <- function(x) {
      s <- gsub("-", "", as.character(x))
      trimws(gsub("\\.[0-9]+$", "", s))
    }
    data.table(
      soc        = clean_id(d$O.NET.SOC.Code),
      skill_name = d$Element.Name,
      importance = suppressWarnings(as.numeric(as.character(d$Data.Value)))
    )
  }), fill = TRUE)
}

path_t0 <- "data/raw/onet/db_15_1"
path_t1 <- "data/raw/onet/db_29_2_text"
if (file.exists("R/99_paths_local.R")) {
  source("R/99_paths_local.R")
  if (exists("PATH_ONET_2015")) path_t0 <- PATH_ONET_2015
  if (exists("PATH_ONET_2024")) path_t1 <- PATH_ONET_2024
}

im_t0 <- read_onet_im(path_t0)[!is.na(importance)]
im_t1 <- read_onet_im(path_t1)[!is.na(importance)]

# Apply SOC crosswalk (2010 to 2019)
cw_raw <- fread("data/crosswalk/2010_to_2019_Crosswalk.csv")
clean_id <- function(x) trimws(gsub("\\.[0-9]+$", "", gsub("-", "", as.character(x))))
cw <- unique(cw_raw[, .(
  soc10 = clean_id(`O*NET-SOC 2010 Code`),
  soc19 = clean_id(`O*NET-SOC 2019 Code`)
)])
rm(cw_raw)

apply_cw <- function(im) {
  im <- merge(im, cw, by.x = "soc", by.y = "soc10",
              all.x = FALSE, allow.cartesian = TRUE)
  im[, soc := soc19][, soc19 := NULL]
  im
}
im_t0 <- apply_cw(im_t0)
im_t1 <- apply_cw(im_t1)

# Compute RCA
compute_rca <- function(im) {
  total       <- im[, .(occ_total   = sum(importance, na.rm = TRUE)), by = soc]
  skill_total <- im[, .(skill_total = sum(importance, na.rm = TRUE)), by = skill_name]
  grand_total <- sum(im$importance, na.rm = TRUE)
  im2 <- merge(im,  total,       by = "soc")
  im2 <- merge(im2, skill_total, by = "skill_name")
  im2[, rca := (importance / occ_total) / (skill_total / grand_total)]
  im2[, .(soc, skill_name, rca)]
}

message("  Computing RCA at t0 and t1...")
rca_t0 <- compute_rca(im_t0)
rca_t1 <- compute_rca(im_t1)
rm(im_t0, im_t1, cw); gc()

# Attach archetype
rca_t0 <- merge(rca_t0, archetype_lkp, by = "skill_name", all.x = TRUE)
rca_t1 <- merge(rca_t1, archetype_lkp, by = "skill_name", all.x = TRUE)
rca_t0 <- rca_t0[atc_archetype %in% CLASS_LEVELS]
rca_t1 <- rca_t1[atc_archetype %in% CLASS_LEVELS]

message(sprintf("  RCA t0: %s rows | %d occupations | %d skills",
                format(nrow(rca_t0), big.mark = ","),
                uniqueN(rca_t0$soc),
                uniqueN(rca_t0$skill_name)))

# ==============================================================================
# Step 3 — Helper: binned diffusion rates with bootstrap CI
# ==============================================================================
bin_diffusion <- function(risk_dt, gap_breaks, n_boot, seed) {
  risk_dt[, gap_bin := cut(pc1_gap, breaks = gap_breaks,
                            include.lowest = TRUE, labels = FALSE)]
  risk_dt[, gap_mid := gap_breaks[gap_bin] +
              (gap_breaks[gap_bin + 1L] - gap_breaks[gap_bin]) / 2]
  set.seed(seed)
  risk_dt[!is.na(gap_bin), {
    v   <- diffusion
    obs <- mean(v, na.rm = TRUE)
    bts <- replicate(n_boot,
                     mean(sample(v, length(v), replace = TRUE), na.rm = TRUE))
    list(
      rate  = obs,
      ci_lo = quantile(bts, 0.025),
      ci_hi = quantile(bts, 0.975),
      n     = .N
    )
  }, by = .(gap_bin, gap_mid, atc_archetype)]
}

# ==============================================================================
# Step 4 — Main loop: iterate over thresholds
# ==============================================================================
message(">>> Step 4: Computing diffusion rates by threshold...")

results_list <- list()

for (thr in THRESHOLDS) {
  message(sprintf("\n  --- Threshold RCA >= %.2f ---", thr))

  spec_t0    <- rca_t0[rca >= thr, .(soc, skill_name, atc_archetype)]
  nonspec_t0 <- rca_t0[rca <  thr, .(soc, skill_name)]
  spec_t1    <- rca_t1[rca >= thr, .(soc, skill_name)]
  nonspec_t1 <- rca_t1[rca <  thr, .(soc, skill_name)]

  # --- Adoption ---
  message("    Building adoption risk set...")
  adopt <- merge(
    spec_t0[,    .(source = soc, skill_name, atc_archetype)],
    nonspec_t0[, .(target = soc, skill_name)],
    by = "skill_name", allow.cartesian = TRUE
  )
  adopt <- adopt[source != target]
  adopt[spec_t1, on = .(target = soc, skill_name), diffusion := 1L]
  adopt[is.na(diffusion), diffusion := 0L]
  adopt[scores_lkp, on = .(source = occ), s_pc1 := i.s_pc1]
  adopt[scores_lkp, on = .(target = occ), t_pc1 := i.s_pc1]
  adopt <- adopt[!is.na(s_pc1) & !is.na(t_pc1)]
  adopt[, pc1_gap := t_pc1 - s_pc1]
  message(sprintf("      n = %s | rate = %.4f",
                  format(nrow(adopt), big.mark = ","), mean(adopt$diffusion)))

  bin_a <- bin_diffusion(adopt, gap_breaks, N_BOOT, SEED_FIG)
  bin_a[, `:=`(threshold = thr, flow = "Adoption")]
  results_list[[paste0("adopt_", thr)]] <- bin_a
  rm(adopt, bin_a); gc()

  # --- Abandonment ---
  message("    Building abandonment risk set...")
  lost_skills <- merge(
    spec_t0[,    .(source = soc, skill_name, atc_archetype)],
    nonspec_t1[, .(soc,          skill_name)],
    by.x = c("source", "skill_name"),
    by.y = c("soc",    "skill_name")
  )
  had_t0 <- spec_t0[, .(target = soc, skill_name)]
  aband <- merge(
    lost_skills[, .(source, skill_name, atc_archetype)],
    had_t0[,      .(target, skill_name)],
    by = "skill_name", allow.cartesian = TRUE
  )
  aband <- aband[source != target]
  aband[nonspec_t1, on = .(target = soc, skill_name), diffusion := 1L]
  aband[is.na(diffusion), diffusion := 0L]
  aband[scores_lkp, on = .(source = occ), s_pc1 := i.s_pc1]
  aband[scores_lkp, on = .(target = occ), t_pc1 := i.s_pc1]
  aband <- aband[!is.na(s_pc1) & !is.na(t_pc1)]
  aband[, pc1_gap := t_pc1 - s_pc1]
  message(sprintf("      n = %s | rate = %.4f",
                  format(nrow(aband), big.mark = ","), mean(aband$diffusion)))

  bin_b <- bin_diffusion(aband, gap_breaks, N_BOOT, SEED_FIG)
  bin_b[, `:=`(threshold = thr, flow = "Abandonment")]
  results_list[[paste0("aband_", thr)]] <- bin_b
  rm(aband, bin_b, lost_skills, had_t0); gc()

  rm(spec_t0, nonspec_t0, spec_t1, nonspec_t1); gc()
}

all_res <- rbindlist(results_list)
rm(results_list, rca_t0, rca_t1); gc()

# Save intermediate CSVs
fwrite(all_res[flow == "Adoption"],
       file.path(SI_TABLES, "table_rca_threshold_adoption_rates.csv"))
fwrite(all_res[flow == "Abandonment"],
       file.path(SI_TABLES, "table_rca_threshold_abandonment_rates.csv"))
message("\n  Saved: table_rca_threshold_adoption_rates.csv")
message("  Saved: table_rca_threshold_abandonment_rates.csv")

# ==============================================================================
# Step 5 — Figure
# ==============================================================================
message(">>> Step 5: Building figure...")

thr_levels <- sprintf("RCA \u2265 %.2f", THRESHOLDS)
all_res[, atc_f   := factor(atc_archetype, CLASS_LEVELS, CLASS_LABELS)]
all_res[, thr_lab := factor(sprintf("RCA \u2265 %.2f", threshold),
                             levels = thr_levels)]
all_res[, flow_f  := factor(flow, c("Adoption", "Abandonment"))]

theme_si <- theme_classic(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12,
                                      margin = margin(t = 2, b = 2)),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 10, colour = "grey10"),
    panel.border       = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.7),
    axis.line          = element_blank(),
    axis.ticks         = element_line(linewidth = 0.7),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    legend.key.width   = unit(1.4, "lines"),
    panel.spacing.x    = unit(0.6, "lines"),
    panel.spacing.y    = unit(0.8, "lines"),
    plot.margin        = margin(6, 8, 4, 8)
  )

fig_S1 <- ggplot(all_res,
                  aes(x      = gap_mid,
                      y      = rate,
                      colour = atc_f,
                      fill   = atc_f)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.3, shape = 16) +
  geom_vline(xintercept = 0, colour = "grey50",
             linewidth = 0.35, linetype = "dashed") +
  # Rows = flow, columns = threshold, free Y per row
  facet_grid(flow_f ~ thr_lab, scales = "free_y") +
  scale_colour_manual(values = CLASS_COLOURS, name = NULL) +
  scale_fill_manual(values   = CLASS_COLOURS, name = NULL) +
  scale_x_continuous(expand  = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Status gap (target \u2212 source, PC1)",
    y = "Diffusion rate"
  ) +
  theme_si +
  guides(
    colour = guide_legend(nrow = 1,
                          override.aes = list(linewidth = 1.0)),
    fill   = "none"
  )

ggsave(file.path(SI_FIGS, "fig_S1_rca_threshold.pdf"),
       fig_S1, width = 12, height = 7, units = "in", device = cairo_pdf)
ggsave(file.path(SI_FIGS, "fig_S1_rca_threshold.png"),
       fig_S1, width = 12, height = 7, units = "in", dpi = 300)

message("  Saved: fig_S1_rca_threshold.pdf / .png")
message("\n>>> SI_S2_fig_S1_rca_threshold.R complete.")