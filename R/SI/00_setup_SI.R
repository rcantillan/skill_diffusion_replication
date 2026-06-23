# ==============================================================================
# 00_setup_SI.R
#
# Common setup for all SI scripts.
# source("R/SI/00_setup_SI.R") at the start of each session.
#
# DESIGN:
#   Two-domain model (Cognitive / Physical).
#   Status gap: G_ij = status_target - status_source
#   50% subsample of source occupations (seed 42).
# ==============================================================================

library(data.table)
library(fixest)

SI_TABLES <- file.path("output", "tables",  "si")
SI_FIGS   <- file.path("output", "figures", "si")
SI_MODELS <- file.path("output", "models",  "si")
for (d in c(SI_TABLES, SI_FIGS, SI_MODELS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

SAMPLE_FRAC <- 0.50
SEED        <- 42L
B_PERM      <- 1000L
CUTOFF_GRID <- c(-1.00, -0.50, -0.25, 0.00, 0.25, 0.50, 1.00)

# ==============================================================================
# load_flow()
# ==============================================================================
load_flow <- function(flow = c("adoption", "abandonment"),
                      seed = SEED,
                      frac = SAMPLE_FRAC) {
  flow <- match.arg(flow)
  message(sprintf("\n>>> load_flow('%s', seed=%d, frac=%.0f%%)...",
                  flow, seed, frac * 100))
  rds_path    <- if (flow == "adoption")
                   "data/derived/riskset_adoption.rds" else
                   "data/derived/riskset_abandonment.rds"
  scores_path <- "output/tables/main/occ_status_scores.csv"
  outcome_col <- if (flow == "adoption") "diffusion" else "abandonment"
  stopifnot(file.exists(rds_path), file.exists(scores_path))

  dt <- readRDS(rds_path); setDT(dt)
  keep <- c("source","target","skill_name",outcome_col,
            "domain","structural_distance")
  dt <- dt[, ..keep]
  dt[, domain := factor(as.character(domain),
                        levels = c("Cognitive","Physical"))]
  dt <- dt[!is.na(domain) & !is.na(structural_distance)]
  gc()

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
  dt[, c("s_pc1","t_pc1","pc1_gap") := NULL]

  qs <- quantile(scores$status_pc1, probs = 0:5/5, na.rm = TRUE)
  scores[, status_q := as.integer(cut(status_pc1, breaks = qs,
                                       labels = 1:5,
                                       include.lowest = TRUE))]
  dt[scores, on = .(source = occ), status_q_source := i.status_q]
  dt <- dt[!is.na(status_q_source)]
  rm(scores); gc()

  ckpt <- file.path(SI_MODELS, sprintf("sources_seed%d.rds", seed))
  if (file.exists(ckpt)) {
    sources_sample <- readRDS(ckpt)
    message(sprintf("  Sources from disk: %d", length(sources_sample)))
  } else {
    set.seed(seed)
    all_src <- unique(dt$source)
    sources_sample <- sample(all_src, size = round(length(all_src) * frac))
    saveRDS(sources_sample, ckpt)
    message(sprintf("  Sources generated: %d of %d",
                    length(sources_sample), length(all_src)))
  }
  dt_s <- dt[source %in% sources_sample]
  message(sprintf("  Sample: %s triads | %d src | %d tgt | %d skills",
                  format(nrow(dt_s), big.mark=","),
                  uniqueN(dt_s$source),
                  uniqueN(dt_s$target),
                  uniqueN(dt_s$skill_name)))
  message(sprintf("  Rate %s: %.4f | Upward: %.1f%%",
                  outcome_col,
                  mean(dt_s[[outcome_col]], na.rm=TRUE),
                  100*mean(dt_s$up_dummy)))
  rm(dt); gc()

  fml <- as.formula(sprintf(
    "%s ~ (up_dummy + pc1_up + pc1_down + structural_distance) : domain",
    outcome_col))
  list(dt = dt_s, fml = fml, outcome = outcome_col, flow = flow)
}

# ==============================================================================
# extract_coefs_2d()
# Uses exact string matching (ct$term == p) — NOT grepl —
# to avoid regex interpretation of ":" in coefficient names.
# ==============================================================================
extract_coefs_2d <- function(m, panel_label, flow_label) {
  ct <- as.data.table(m$coeftable, keep.rownames = "term")
  # Normalize column names regardless of fixest version
  nm <- names(ct)
  if ("Estimate"   %in% nm) setnames(ct, "Estimate",   "estimate")
  if ("Std. Error" %in% nm) setnames(ct, "Std. Error", "std_error")
  if ("z value"    %in% nm) setnames(ct, "z value",    "z")
  if ("Pr(>|z|)"   %in% nm) setnames(ct, "Pr(>|z|)",  "p")

  # Exact match — avoids grepl regex on ":"
  find_coef <- function(patterns) {
    for (p in patterns) {
      idx <- which(ct[["term"]] == p)
      if (length(idx) > 0)
        return(c(ct[["estimate"]][idx[1L]],
                 ct[["std_error"]][idx[1L]]))
    }
    c(NA_real_, NA_real_)
  }

  params <- list(
    b_up_Cog  = find_coef(c("pc1_up:domainCognitive",
                             "domainCognitive:pc1_up")),
    b_dn_Cog  = find_coef(c("pc1_down:domainCognitive",
                             "domainCognitive:pc1_down")),
    kappa_Cog = find_coef(c("up_dummy:domainCognitive",
                             "domainCognitive:up_dummy")),
    delta_Cog = find_coef(c("structural_distance:domainCognitive",
                             "domainCognitive:structural_distance")),
    b_up_Phy  = find_coef(c("pc1_up:domainPhysical",
                             "domainPhysical:pc1_up")),
    b_dn_Phy  = find_coef(c("pc1_down:domainPhysical",
                             "domainPhysical:pc1_down")),
    kappa_Phy = find_coef(c("up_dummy:domainPhysical",
                             "domainPhysical:up_dummy")),
    delta_Phy = find_coef(c("structural_distance:domainPhysical",
                             "domainPhysical:structural_distance"))
  )

  out <- data.table(
    flow      = flow_label,
    panel     = panel_label,
    coef      = names(params),
    estimate  = vapply(params, function(x) x[1L], numeric(1L)),
    std_error = vapply(params, function(x) x[2L], numeric(1L))
  )
  out[, sig := fcase(
    is.na(estimate),                 "",
    abs(estimate/std_error) >= 3.29, "***",
    abs(estimate/std_error) >= 2.58, "**",
    abs(estimate/std_error) >= 1.96, "*",
    abs(estimate/std_error) >= 1.65, "~",
    default = "ns"
  )]
  out
}

message("[OK] 00_setup_SI.R ready.")
message(sprintf("     B_PERM=%d | SAMPLE_FRAC=%.0f%%", B_PERM, SAMPLE_FRAC*100))
