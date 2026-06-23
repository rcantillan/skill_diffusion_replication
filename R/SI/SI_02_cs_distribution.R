# ==============================================================================
# SI_02_cs_distribution.R
#
# Section S3: Nestedness — Fig. S2
#
# Produces Fig. S2: Distribution of nestedness contributions (c_s) by skill
# domain. Kernel density estimates of the standardized leave-one-out nestedness
# contribution c_s for socio-cognitive (teal) and sensory/physical (grey)
# requirements. The vertical dashed line marks the within-domain median used
# to classify archetypes into SC_Scaffolding (c_s >= median within Cognitive)
# and SC_Specialized (c_s < median within Cognitive); all sensory/physical
# requirements map to Physical_Terminal regardless of c_s position.
#
# This figure is purely descriptive — no estimation required.
#
# Input:  data/derived/skill_cs_scores.rds
#         data/derived/riskset_adoption.rds (for domain lookup)
#
# Output: output/figures/si/fig_S2_cs_distribution.pdf / .png
#         output/tables/si/table_cs_distribution_stats.csv
# ==============================================================================

source("R/SI/00_setup_SI.R")
library(ggplot2)

# ==============================================================================
# Load c_s scores and domain labels
# ==============================================================================
message(">>> Loading c_s scores...")

cs <- readRDS("data/derived/skill_cs_scores.rds")
setDT(cs)
message(sprintf("  Skills with c_s: %d | NAs: %d",
                cs[!is.na(cs), .N], cs[is.na(cs), .N]))

# Domain lookup from the adoption risk set
dt <- readRDS("data/derived/riskset_adoption.rds")
setDT(dt)
domain_lkp <- unique(dt[, .(skill_name, domain)])[!is.na(domain)]
rm(dt); gc()

# Defensive domain-text normalization to prevent NA median
domain_lkp[grepl("(?i)cognitive", domain), domain := "Cognitive"]
domain_lkp[grepl("(?i)physical", domain), domain := "Physical"]

# Merge domain onto c_s scores
cs <- merge(cs, domain_lkp, by = "skill_name", all.x = TRUE)
cs <- cs[!is.na(domain) & !is.na(cs)]

message(sprintf("  Cognitive: %d skills | Physical: %d skills",
                cs[domain == "Cognitive", .N],
                cs[domain == "Physical",  .N]))

# ==============================================================================
# Archetype classification thresholds
# ==============================================================================
med_cog <- cs[domain == "Cognitive", median(cs, na.rm = TRUE)]
message(sprintf("  Within-domain median c_s (Cognitive): %.4f", med_cog))

cs[, archetype := fcase(
  domain == "Cognitive" & cs >= med_cog, "SC_Scaffolding",
  domain == "Cognitive" & cs <  med_cog, "SC_Specialized",
  domain == "Physical",                   "Physical_Terminal"
)]

# ==============================================================================
# Descriptive statistics
# ==============================================================================
stats_tbl <- cs[, .(
  n        = .N,
  mean_cs  = round(mean(cs,   na.rm = TRUE), 3),
  med_cs   = round(median(cs, na.rm = TRUE), 3),
  sd_cs    = round(sd(cs,     na.rm = TRUE), 3),
  min_cs   = round(min(cs,    na.rm = TRUE), 3),
  max_cs   = round(max(cs,    na.rm = TRUE), 3),
  pct_pos  = round(100 * mean(cs > 0, na.rm = TRUE), 1)
), by = domain]

fwrite(stats_tbl, file.path(SI_TABLES, "table_cs_distribution_stats.csv"))

message("\n  c_s statistics by domain:")
print(stats_tbl)

# Also report archetype breakdown
arch_tbl <- cs[, .N, by = .(domain, archetype)][order(domain, archetype)]
arch_tbl[, pct := round(100 * N / sum(N)), by = domain]
message("\n  Archetype breakdown:")
print(arch_tbl)

# ==============================================================================
# Figure S2 — c_s distribution by domain
# ==============================================================================
message("\n>>> Generating Fig. S2...")

# Domain colors: teal for Cognitive, grey for Physical (per SI convention)
DOM_COLS <- c("Cognitive" = "#008280", "Physical" = "#8c8c8c")
DOM_LBLS <- c("Cognitive" = "Socio-cognitive",
              "Physical"  = "Sensory/physical")

theme_si <- theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14, colour = "grey10"),
    panel.border     = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.8),
    axis.line        = element_blank(),
    axis.ticks       = element_line(linewidth = 0.8),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 14),
    legend.key.width = unit(1.4, "lines"),
    plot.margin      = margin(8, 10, 6, 8)
  )

fig_S2 <- ggplot(cs, aes(x = cs, fill = domain, colour = domain)) +
  geom_density(alpha = 0.35, linewidth = 0.8, bw = "nrd0") +
  
  geom_vline(xintercept = med_cog,
             colour = "grey30", linewidth = 0.6,
             linetype = "dashed") +
             
  annotate("text", x = med_cog - 0.5, y = Inf,
           label = "SC_specialized", hjust = 1, vjust = 1.5,
           size = 5.0, colour = "#008280", fontface = "italic",
           family = "Helvetica") +
  annotate("text", x = med_cog + 0.5, y = Inf,
           label = "SC_general", hjust = 0, vjust = 1.5,
           size = 5.0, colour = "#008280", fontface = "italic",
           family = "Helvetica") +
           
  # Center "Sensory_physical" label over the grey curve peak (approx at -5)
  annotate("text", x = -5.0, y = Inf,
           label = "Sensory_physical", hjust = 0.5, vjust = 1.5,
           size = 5.0, colour = "#8c8c8c", fontface = "italic",
           family = "Helvetica") +
           
  scale_fill_manual(  values = DOM_COLS, labels = DOM_LBLS, name = NULL) +
  scale_colour_manual(values = DOM_COLS, labels = DOM_LBLS, name = NULL) +
  
  # Limits set here so tails reach zero;
  # extended to 19 to leave room for the "SC_Scaffolding" label
  scale_x_continuous(
    limits = c(-12, 19), 
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  # Increased upper margin so the tall green curve does not
  # visually collide with the annotation text
  scale_y_continuous(
    expand = expansion(mult = c(0.00, 0.20))
  ) +
  
  # coord_cartesian removed — limits handled by scale_x_continuous above
  
  labs(
    x        = expression("Standardized nestedness contribution (" * italic(c)[s] * ")"),
    y        = "Density",
    tag      = NULL,
    title    = NULL,
    subtitle = NULL,
    caption  = NULL
  ) +
  theme_si +
  guides(fill   = guide_legend(nrow = 1),
         colour = guide_legend(nrow = 1))

ggsave(file.path(SI_FIGS, "fig_S2_cs_distribution.pdf"),
       fig_S2, width = 7, height = 5,
       units = "in", device = cairo_pdf)
ggsave(file.path(SI_FIGS, "fig_S2_cs_distribution.png"),
       fig_S2, width = 7, height = 5,
       units = "in", dpi = 300)

message("  Saved: fig_S2_cs_distribution.pdf / .png")
message("\n>>> SI_02_cs_distribution.R complete.")
message("    Next: SI_03_baseline.R")

print(stats_tbl)
