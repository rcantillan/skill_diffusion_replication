# ==============================================================================
# 05b_gravity_abandonment.R
#
# Cloglog gravity models for ABANDONMENT using the PC1 status index.
#
# DESIGN:
#   source i = co-specialized reference occupation (RCA>1 at t0)
#   target j = focal occupation that may lose the skill (RCA>1 at t0)
#   outcome  = abandonment_j = 1 if j loses RCA>1 between t0 and t1
#   gap      = G_ij = status_j - status_i (positive = j has higher status)
#
# STATUS:
#   status_pc1 per occupation from occ_status_scores.csv
#   G_ij = status_target - status_source
#
# IDENTIFICATION:
#   Panel A (Source FE)  = PRIMARY       -> within-source
#   Panel B (Target FE)  = CORROBORATION -> within-target
#
# DTC PREDICTIONS (abandonment):
#   SC archetypes:     beta_up < 0  (higher-status targets RETAIN cognitive)
#   Physical_Terminal: beta_up > 0  (higher-status targets SHED physical)
#
# Input:  data/derived/riskset_abandonment.rds
#         output/tables/main/occ_status_scores.csv
#         output/tables/main/pca_status_decision.csv
#
# Output: output/models/model_abandonment_panelA.rds
#         output/models/model_abandonment_panelB.rds
#         output/tables/abandonment/coefs_pc1_abandonment.csv
#
# Next: 06_projections.R
# ==============================================================================

gc()
library(data.table)
library(fixest)

if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

FORCE_REFIT <- TRUE

# Paths
dt_path     <- if (exists("PATH_ABANDONMENT")) PATH_ABANDONMENT else
               "data/derived/riskset_abandonment.rds"
scores_path <- "output/tables/main/occ_status_scores.csv"
dec_path    <- "output/tables/main/pca_status_decision.csv"
out_models  <- "output/models"
out_tables  <- "output/tables/abandonment"

for (d in c(out_models, out_tables))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

stopifnot("riskset_abandonment.rds not found"  = file.exists(dt_path))
stopifnot("occ_status_scores.csv not found"    = file.exists(scores_path))
stopifnot("pca_status_decision.csv not found"  = file.exists(dec_path))

pca_dec <- fread(dec_path)
message(sprintf("PC1: %.1f%% var | flip: %s | FORCE_REFIT: %s",
                pca_dec$pc1_pct_var,
                ifelse(pca_dec$needs_flip, "YES", "NO"),
                FORCE_REFIT))

MAIN_ARCHETYPES <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")

# ==============================================================================
# Step 1 — Load data and verify
# ==============================================================================
message("\n>>> Step 1: Loading data...")

dt <- readRDS(dt_path)
setDT(dt)
stopifnot("abandonment" %in% names(dt))
stopifnot("cs"          %in% names(dt))
stopifnot("domain"      %in% names(dt))
message(sprintf("  Rows: %s | Abandonment rate: %.4f",
                format(nrow(dt), big.mark = ","),
                mean(dt$abandonment)))

# atc_archetype comes from 03c_nestedness_merge_ab.R
if (!"atc_archetype" %in% names(dt)) {
  cs_med <- dt[domain == "Cognitive", median(cs, na.rm = TRUE)]
  dt[, atc_archetype := fcase(
    domain == "Cognitive" & cs >= cs_med, "SC_Scaffolding",
    domain == "Cognitive" & cs <  cs_med, "SC_Specialized",
    domain == "Physical",                  "Physical_Terminal"
  )]
}
dt[, atc_archetype := factor(atc_archetype, levels = MAIN_ARCHETYPES)]
dt <- dt[atc_archetype %in% MAIN_ARCHETYPES]
dt <- dt[, .(source, target, skill_name, abandonment,
             atc_archetype, structural_distance)]
dt <- dt[!is.na(structural_distance)]
gc()

message(sprintf("  Post-filter: %s rows", format(nrow(dt), big.mark = ",")))
print(dt[, .N, by = atc_archetype][order(atc_archetype)])

# ==============================================================================
# Step 2 — Merge PC1 status and build G_ij
# ==============================================================================
message("\n>>> Step 2: Building status gap G_ij...")

scores <- fread(scores_path)
setDT(scores)
scores[, occ := as.character(occ)]

dt[, source := as.character(source)]
dt[, target := as.character(target)]
dt[scores, on = .(source = occ), status_source := i.status_pc1]
dt[scores, on = .(target = occ), status_target := i.status_pc1]

n_na <- dt[is.na(status_source) | is.na(status_target), .N]
message(sprintf("  Dyads without status: %d (%.2f%%)", n_na, 100*n_na/nrow(dt)))
dt <- dt[!is.na(status_source) & !is.na(status_target)]

# G_ij = status_target - status_source
dt[, status_pc1 := status_target - status_source]
dt[, pc1_up     := pmax(0,  status_pc1)]
dt[, pc1_down   := pmin(0,  status_pc1)]
dt[, pc1_dummy  := fifelse(status_pc1 > 0, 1L, 0L)]
dt[, c("status_source", "status_target") := NULL]

message(sprintf("  G_ij range: [%.3f, %.3f]",
                min(dt$status_pc1), max(dt$status_pc1)))
message(sprintf("  Upward dyads (G_ij>0): %.1f%%", 100*mean(dt$pc1_dummy)))

# DTC diagnostic for abandonment
tab <- dt[, .(rate = mean(abandonment), n = .N),
          by = .(atc_archetype,
                 direction = fifelse(pc1_dummy == 1L, "up", "down"))]
setorder(tab, atc_archetype, direction)
message("\n  DTC diagnostic (abandonment):")
print(tab)
wide <- dcast(tab, atc_archetype ~ direction, value.var = "rate")
if (all(c("up","down") %in% names(wide))) {
  wide[, ratio := round(up/down, 3)]
  message("  Ratio up/down:")
  print(wide)
  message("  Expected — SC: ratio<1 | Physical: ratio>1")
}

dt[, status_pc1 := NULL]
gc()

# ==============================================================================
# Step 3 — Rug segments
# ==============================================================================
message("\n>>> Step 3: Rug segments...")

rug_full <- dt[, .(x = pc1_up + pc1_down,
                   archetype = as.character(atc_archetype))]
set.seed(42)
rug_samp <- rug_full[sample(.N, min(60000L, .N))]
rug_samp[, archetype := factor(archetype, levels = MAIN_ARCHETYPES)]
saveRDS(rug_samp, file.path(out_models, "rug_segs_abandonment.rds"))
message("  Saved: rug_segs_abandonment.rds")
rm(rug_full, rug_samp); gc()

# ==============================================================================
# Step 4 — Formula
# ==============================================================================
formula_m <- abandonment ~
  (pc1_dummy + pc1_up + pc1_down + structural_distance):atc_archetype

message(sprintf("\n  Formula: %s", deparse(formula_m)))
message("  Panel A: FE(source, skill_name) — PRIMARY")
message("  Panel B: FE(target, skill_name) — CORROBORATION")

# ==============================================================================
# Step 5 — Panel A: Source + Skill FE  [PRIMARY]
# ==============================================================================
message("\n>>> Step 5: Panel A — Source FE (PRIMARY)...")

path_A <- file.path(out_models, "model_abandonment_panelA.rds")

if (FORCE_REFIT || !file.exists(path_A)) {
  t0  <- proc.time()
  m_A <- feglm(
    formula_m,
    data      = dt,
    family    = binomial("cloglog"),
    fixef     = c("source", "skill_name"),
    cluster   = c("source", "target", "skill_name"),
    mem.clean = TRUE,
    nthreads  = 0,
    lean      = TRUE
  )
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message(sprintf("  Panel A: %.1f min", elapsed))
  saveRDS(m_A, path_A)
  rm(m_A)
} else {
  message("  Loading from disk (FORCE_REFIT=FALSE)")
}
gc()

# ==============================================================================
# Step 6 — Panel B: Target + Skill FE  [CORROBORATION]
# ==============================================================================
message("\n>>> Step 6: Panel B — Target FE (CORROBORATION)...")

path_B <- file.path(out_models, "model_abandonment_panelB.rds")

if (FORCE_REFIT || !file.exists(path_B)) {
  t0  <- proc.time()
  m_B <- feglm(
    formula_m,
    data      = dt,
    family    = binomial("cloglog"),
    fixef     = c("target", "skill_name"),
    cluster   = c("source", "target", "skill_name"),
    mem.clean = TRUE,
    nthreads  = 0,
    lean      = TRUE
  )
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message(sprintf("  Panel B: %.1f min", elapsed))
  saveRDS(m_B, path_B)
  rm(m_B)
} else {
  message("  Loading from disk (FORCE_REFIT=FALSE)")
}
gc()

# ==============================================================================
# Step 7 — Extract coefficients
# ==============================================================================
message("\n>>> Step 7: Extracting coefficients...")

extract_coefs <- function(model_path, panel_label) {
  m  <- readRDS(model_path)
  ct <- as.data.table(m$coeftable, keep.rownames = "term")
  rm(m); gc()
  nms <- names(ct)
  if ("Estimate"   %in% nms) setnames(ct, "Estimate",   "estimate")
  if ("Std. Error" %in% nms) setnames(ct, "Std. Error", "std_error")
  if ("z value"    %in% nms) setnames(ct, "z value",    "z")
  if ("Pr(>|z|)"   %in% nms) setnames(ct, "Pr(>|z|)",  "p")

  ct[, panel     := panel_label]
  ct[, archetype := fcase(
    grepl("SC_Scaffolding",    term), "SC_Scaffolding",
    grepl("SC_Specialized",    term), "SC_Specialized",
    grepl("Physical_Terminal", term), "Physical_Terminal",
    default = NA_character_
  )]
  ct[, var := fcase(
    grepl("pc1_dummy",           term) & !grepl("pc1_up|pc1_down", term), "kappa",
    grepl("pc1_up",              term),                                    "b_up",
    grepl("pc1_down",            term),                                    "b_dn",
    grepl("structural_distance", term),                                    "delta",
    default = "other"
  )]
  ct[, sig := fcase(
    p < 0.001, "***", p < 0.01, "**",
    p < 0.05,  "*",   p < 0.10, "~",
    default = "ns"
  )]
  ct[!is.na(archetype) & var != "other",
     .(panel, var, archetype,
       coef = round(estimate,  4),
       se   = round(std_error, 4),
       p    = round(p, 4), sig)]
}

coefs <- rbind(
  extract_coefs(path_A, "Panel A"),
  extract_coefs(path_B, "Panel B")
)

message("\n  Coefficients b_up / b_dn by panel and archetype:")
print(coefs[var %in% c("b_up","b_dn")][order(panel, archetype, var)])

atc_dt <- merge(
  coefs[var == "b_up", .(panel, archetype, b_up = coef)],
  coefs[var == "b_dn", .(panel, archetype, b_dn = coef)],
  by = c("panel","archetype")
)[, dtc_index := b_up - b_dn]

message("\n  DTC index (b_up - b_dn):")
print(atc_dt[order(panel, archetype)])

# Panel A verification
message("\n  CHECK Panel A:")
b_up_sc  <- coefs[panel=="Panel A" & var=="b_up" &
                  archetype %in% c("SC_Scaffolding","SC_Specialized"),
                  mean(coef)]
b_up_phy <- coefs[panel=="Panel A" & var=="b_up" &
                  archetype == "Physical_Terminal", coef[1]]
message(paste0("  beta_up SC mean  = ", round(b_up_sc,  4),
               "  [expected: < 0]  ",
               ifelse(b_up_sc  < 0, "[OK]", "[ATTENTION]")))
message(paste0("  beta_up Physical = ", round(b_up_phy, 4),
               "  [expected: > 0]  ",
               ifelse(b_up_phy > 0, "[OK]", "[ATTENTION]")))

fwrite(coefs,  file.path(out_tables, "coefs_pc1_abandonment.csv"))
fwrite(atc_dt, file.path(out_tables, "atc_pc1_abandonment.csv"))
message("\n  Saved: coefs_pc1_abandonment.csv | atc_pc1_abandonment.csv")

gc()
message("\n>>> 05b_gravity_abandonment.R complete.")
message("    Next: 06_projections.R")
