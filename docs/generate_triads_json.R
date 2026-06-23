# ==============================================================================
# generate_triads_json.R
#
# Generates triads_explorer.json for the interactive skill diffusion explorer.
#
# For each flow (adoption / abandonment) × distance range (D2D9 / D4D6):
#   - Computes lp (log partial hazard) per triad using Panel A coefficients
#     (source + skill FE) for both flows — Panel A is the correct specification
#     because TARGET is the focal unit in both adoption and abandonment risk sets.
#   - Assigns triads to 3×3 matrix cells:
#       rows    = lp tercile (High / Mid / Low) within pc1_zone
#       columns = pc1_zone T1 (downward) / T2 (lateral) / T3 (upward)
#   - Extracts N_TRIADS per cell using random and dyad-unique sampling
#   - Computes skill-level scatter (mean lp × mean gap × frequency)
#
# Output: R/SI/explorer/triads_explorer.json
# ==============================================================================

gc()
library(data.table)
library(jsonlite)

if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

# ==============================================================================
# Parameters
# ==============================================================================
N_TRIADS    <- 20L
MIN_N       <- 3L
SAMPLE_SEED <- 42L

# Fixed-effects panel for BOTH flows:
#   "A" = source + skill FE (preferred; TARGET is the focal unit that adopts
#         or abandons in both risk sets, so source FE is the correct absorber)
#   "B" = target + skill FE
PANEL <- "A"

DECILE_RANGES <- list(
  D2D9 = c(2L, 9L),   # "Full"
  D4D6 = c(4L, 6L)    # "Trimmed" (moderate structural distance only)
)

ARCH_3 <- c("SC_Scaffolding", "SC_Specialized", "Physical_Terminal")

out_dir  <- file.path("R", "SI", "explorer")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_json <- file.path(out_dir, "triads_explorer.json")

# ==============================================================================
# Helpers
# ==============================================================================
clean_id <- function(x) {
  s <- gsub("-", "", as.character(x))
  s <- gsub("\\.[0-9]+$", "", s)
  trimws(s)
}

trunc_name <- function(x, n = 32)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)

ntile_fn <- function(x, n) {
  r <- rank(x, ties.method = "first", na.last = "keep")
  ceiling(n * r / max(r, na.rm = TRUE))
}

# ==============================================================================
# O*NET occupation titles
# ==============================================================================
message(">>> Loading O*NET occupation titles...")
load_occ_titles <- function(folder) {
  f <- file.path(folder, "Occupation Data.txt")
  if (!file.exists(f)) return(NULL)
  d <- fread(f, sep = "\t", quote = "", showProgress = FALSE,
             na.strings = c("NA", "n/a", "", "*"))
  setnames(d, names(d), make.names(names(d)))
  cc <- grep("O.NET.SOC.Code|Code",  names(d), ignore.case=TRUE, value=TRUE)[1]
  ct <- grep("Title|Occupation",     names(d), ignore.case=TRUE, value=TRUE)[1]
  if (is.na(cc) || is.na(ct)) return(NULL)
  d <- d[, .(soc_code  = clean_id(get(cc)),
             occ_title = trimws(as.character(get(ct))))]
  d[!is.na(occ_title) & occ_title != "",
    .(occ_title = first(occ_title)), by = soc_code]
}

occ_titles <- NULL
for (p in c("data/raw/onet/db_29_2_text","data/raw/onet/db_15_1")) {
  if (dir.exists(p)) {
    occ_titles <- load_occ_titles(p)
    if (!is.null(occ_titles) && nrow(occ_titles) > 0) {
      message(sprintf("  %d titles from %s", nrow(occ_titles), p)); break
    }
  }
}
if (is.null(occ_titles)) stop("O*NET titles not found.")
setkey(occ_titles, soc_code)

# ==============================================================================
# Status scores
# ==============================================================================
scores <- fread("output/tables/main/occ_status_scores.csv")
setDT(scores); scores[, occ := as.character(occ)]

# ==============================================================================
# Model coefficients — Panel A (source + skill FE) for BOTH flows
# Panel A is correct for both because TARGET is the focal unit in both
# risk sets: target adopts in adoption, target abandons in abandonment.
# ==============================================================================
message(">>> Loading Panel A coefficients...")

get_panelA_coefs <- function(flow) {
  model_path <- sprintf("output/models/si/baseline_%s_%s.rds",
                        ifelse(flow=="adoption","adopt","aband"), PANEL)
  if (!file.exists(model_path)) {
    message(sprintf("  [!] %s not found", model_path))
    return(NULL)
  }
  m  <- readRDS(model_path)
  ct <- as.data.table(m$coeftable, keep.rownames="term")
  nm <- names(ct)
  if ("Estimate"   %in% nm) setnames(ct, "Estimate",   "estimate")
  if ("Std. Error" %in% nm) setnames(ct, "Std. Error", "std_error")
  ct
}

find_coef <- function(ct, patterns) {
  if (is.null(ct)) return(0)
  for (p in patterns) {
    idx <- which(ct[["term"]] == p)
    if (length(idx) > 0) return(ct[["estimate"]][idx[1]])
  }
  0
}

build_coef_tbl <- function(flow) {
  ct <- get_panelA_coefs(flow)
  rbindlist(lapply(ARCH_3, function(arch) {
    data.table(
      flow      = flow,
      archetype = arch,
      kappa = find_coef(ct, c(
        sprintf("up_dummy:domain%s",   ifelse(grepl("Phys",arch),"Physical","Cognitive")),
        sprintf("domain%s:up_dummy",   ifelse(grepl("Phys",arch),"Physical","Cognitive")))),
      b_up  = find_coef(ct, c(
        sprintf("pc1_up:domain%s",     ifelse(grepl("Phys",arch),"Physical","Cognitive")),
        sprintf("domain%s:pc1_up",     ifelse(grepl("Phys",arch),"Physical","Cognitive")))),
      b_dn  = find_coef(ct, c(
        sprintf("pc1_down:domain%s",   ifelse(grepl("Phys",arch),"Physical","Cognitive")),
        sprintf("domain%s:pc1_down",   ifelse(grepl("Phys",arch),"Physical","Cognitive"))))
    )
  }))
}

coef_tbl <- rbind(
  build_coef_tbl("adoption"),
  build_coef_tbl("abandonment")
)

message("  Coefficient table (Panel A, 2-domain):")
print(coef_tbl[, .(flow, archetype,
                   kappa = round(kappa,3),
                   b_up  = round(b_up, 3),
                   b_dn  = round(b_dn, 3))])

# ==============================================================================
# Process one flow
# ==============================================================================
process_flow <- function(flow) {
  message(sprintf("\n>>> Processing flow: %s", flow))

  rds_path    <- sprintf("data/derived/riskset_%s.rds", flow)
  outcome_col <- if (flow=="adoption") "diffusion" else "abandonment"

  dt <- readRDS(rds_path); setDT(dt)
  keep <- c("source","target","skill_name", outcome_col,
            "domain","atc_archetype","structural_distance")
  keep <- intersect(keep, names(dt))
  dt   <- dt[, ..keep]
  dt   <- dt[atc_archetype %in% ARCH_3 & !is.na(structural_distance)]
  if (outcome_col %in% names(dt)) setnames(dt, outcome_col, "event") else dt[, event := NA_integer_]
  dt[, source := as.character(source)]
  dt[, target := as.character(target)]

  # Status gap — keep raw s_pc1 / t_pc1 / pc1_gap throughout
  dt[scores, on=.(source=occ), s_pc1 := i.status_pc1]
  dt[scores, on=.(target=occ), t_pc1 := i.status_pc1]
  dt <- dt[!is.na(s_pc1) & !is.na(t_pc1)]
  dt[, pc1_gap   := t_pc1 - s_pc1]
  dt[, pc1_up    := pmax(0,  pc1_gap)]
  dt[, pc1_down  := pmin(0,  pc1_gap)]
  dt[, pc1_dummy := fifelse(pc1_gap > 0, 1L, 0L)]

  # Compute lp — Panel A 2-domain coefficients
  dt[, lp := NA_real_]
  .fl <- flow  # local copy: avoids collision with the `flow` column in coef_tbl
  for (arch in ARCH_3) {
    cc <- coef_tbl[flow == .fl & archetype == arch]
    if (nrow(cc) == 0) next
    dt[atc_archetype == arch,
       lp := cc$kappa * pc1_dummy + cc$b_up * pc1_up + cc$b_dn * pc1_down]
  }
  dt <- dt[!is.na(lp)]

  # Merge occupation titles
  dt[, sc := clean_id(source)]
  dt[, tc := clean_id(target)]
  dt <- merge(dt, occ_titles, by.x="sc", by.y="soc_code", all.x=TRUE)
  setnames(dt, "occ_title", "source_title")
  dt <- merge(dt, occ_titles, by.x="tc", by.y="soc_code", all.x=TRUE)
  setnames(dt, "occ_title", "target_title")
  dt[is.na(source_title), source_title := sc]
  dt[is.na(target_title), target_title := tc]
  dt[, c("sc","tc","source","target") := NULL]
  gc()

  message(sprintf("  Loaded: %s triads | lp [%.3f, %.3f] | gap [%.3f, %.3f]",
                  format(nrow(dt), big.mark=","),
                  min(dt$lp), max(dt$lp),
                  min(dt$pc1_gap), max(dt$pc1_gap)))

  # Global decile of structural distance
  dt[, dist_decile := ntile_fn(structural_distance, 10)]

  # PC1 zone using GLOBAL quantiles (consistent across decile ranges)
  pc1_q_global <- quantile(dt$pc1_gap, c(0, 1/3, 2/3, 1), na.rm=TRUE)
  dt[, pc1_zone := fcase(
    pc1_gap <= pc1_q_global[2], "T1",
    pc1_gap <= pc1_q_global[3], "T2",
    pc1_gap >  pc1_q_global[3], "T3",
    default = NA_character_
  )]
  dt <- dt[!is.na(pc1_zone)]

  result_list <- list()

  for (dr_name in names(DECILE_RANGES)) {
    dr    <- DECILE_RANGES[[dr_name]]
    dt_dr <- dt[dist_decile >= dr[1] & dist_decile <= dr[2]]
    message(sprintf("  %s: %s triads", dr_name,
                    format(nrow(dt_dr), big.mark=",")))
    if (nrow(dt_dr) == 0) next

    # ── SCATTER: mean lp AND mean gap per skill ────────────────────────────
    # Among REALIZED events only (event == 1). Over the full risk set the
    # per-skill mean gap is ~0 by construction (all directed dyads at risk),
    # which collapses the scatter at x = 0. Fallback: full risk set.
    dt_ev <- dt_dr[event == 1L]
    if (nrow(dt_ev) == 0) dt_ev <- dt_dr
    skill_scatter <- dt_ev[, .(
      lp_mean   = mean(lp,      na.rm=TRUE),
      gap_mean  = mean(pc1_gap, na.rm=TRUE),
      n_triads  = .N,
      archetype = as.character(first(atc_archetype)),
      domain    = as.character(first(domain))
    ), by=skill_name]

    # ── AGGREGATE per skill × zone ─────────────────────────────────────────
    skill_zone <- dt_dr[, {
      # Modal pair among realized events when available (matches the HTML
      # caption: "most frequent real observed occupational pair");
      # fallback to the full risk set if the skill x zone has no events.
      ev_idx <- which(event == 1L)
      pair_src <- if (length(ev_idx) > 0) ev_idx else seq_len(.N)
      tb <- sort(table(paste0(source_title[pair_src],"|||",target_title[pair_src])),
                 decreasing=TRUE)
      modal_pair <- names(tb)[1]
      list(
        lp_mean     = mean(lp, na.rm=TRUE),
        n_dyads     = .N,
        archetype   = as.character(first(atc_archetype)),
        domain      = as.character(first(domain)),
        source_title = trunc_name(trimws(sub("\\|\\|\\|.*","",modal_pair)),32),
        target_title = trunc_name(trimws(sub(".*\\|\\|\\|","",modal_pair)),32)
      )
    }, by=.(skill_name, pc1_zone)]
    skill_zone <- skill_zone[n_dyads >= MIN_N]

    # lp tercile at the SKILL x ZONE level (the unit actually displayed).
    # Computed on lp_mean within each zone, so every zone splits ~evenly
    # into High / Mid / Low and no cell is empty by construction.
    skill_zone[, lp_tercile := {
      r <- frank(lp_mean, ties.method="average", na.last="keep")
      m <- max(r, na.rm=TRUE)
      fcase(r/m > 2/3, "High", r/m > 1/3, "Mid", default="Low")
    }, by = pc1_zone]

    # ── CELL EXTRACTION ────────────────────────────────────────────────────
    cells_random      <- list()
    cells_dyad_unique <- list()

    for (lp_lev in c("High","Mid","Low")) {
      for (zone in c("T1","T2","T3")) {
        key <- paste0(lp_lev,"__",zone)
        sub <- skill_zone[lp_tercile==lp_lev & pc1_zone==zone]

        make_triads <- function(rows) {
          lapply(seq_len(nrow(rows)), function(i) list(
            skill     = rows$skill_name[i],
            lp        = round(rows$lp_mean[i], 3),
            archetype = rows$archetype[i],
            domain    = rows$domain[i],
            source    = rows$source_title[i],
            target    = rows$target_title[i],
            n_dyads   = rows$n_dyads[i]
          ))
        }

        # Random
        set.seed(SAMPLE_SEED)
        cells_random[[key]] <- if (nrow(sub)==0) list() else
          make_triads(setorder(sub[sample(.N, min(N_TRIADS,.N))], -lp_mean))

        # Dyad-unique
        set.seed(SAMPLE_SEED)
        sub_shuf   <- sub[sample(.N)]
        seen_dyads <- character(0)
        selected   <- integer(0)
        for (i in seq_len(nrow(sub_shuf))) {
          dk <- paste0(sub_shuf$source_title[i],"|||",sub_shuf$target_title[i])
          if (!dk %in% seen_dyads) {
            selected   <- c(selected, i)
            seen_dyads <- c(seen_dyads, dk)
          }
          if (length(selected) >= N_TRIADS) break
        }
        cells_dyad_unique[[key]] <- if (length(selected)==0) list() else
          make_triads(setorder(sub_shuf[selected], -lp_mean))
      }
    }

    result_list[[dr_name]] <- list(
      n_triads          = nrow(dt_dr),
      dist_range        = list(
        min = round(min(dt_dr$structural_distance), 3),
        max = round(max(dt_dr$structural_distance), 3)),
      pc1_cutpoints     = list(
        T1_max = round(pc1_q_global[2], 3),
        T2_max = round(pc1_q_global[3], 3)),
      scatter           = lapply(seq_len(nrow(skill_scatter)), function(i) list(
        skill     = skill_scatter$skill_name[i],
        lp        = round(skill_scatter$lp_mean[i],  3),
        gap       = round(skill_scatter$gap_mean[i],  3),
        n         = skill_scatter$n_triads[i],
        archetype = skill_scatter$archetype[i],
        domain    = skill_scatter$domain[i]
      )),
      cells_random      = cells_random,
      cells_dyad_unique = cells_dyad_unique
    )
  }

  result_list
}

# ==============================================================================
# Build and write JSON
# ==============================================================================
message("\n>>> Building JSON...")

full_data <- list(
  meta = list(
    generated         = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    panel             = PANEL,
    n_triads_per_cell = N_TRIADS,
    min_n             = MIN_N,
    seed              = SAMPLE_SEED,
    panel             = "Panel A (source + skill FE) for both flows",
    archetypes        = ARCH_3,
    decile_ranges     = names(DECILE_RANGES),
    decile_labels     = list(
      D2D9 = "Full",
      D4D6 = "Trimmed (D4\u2013D6)"
    ),
    flows             = c("adoption", "abandonment"),
    sampling          = c("random", "dyad_unique"),
    pc1_variance      = tryCatch({
      fread("output/tables/main/pca_status_decision.csv")$pc1_pct_var[1]
    }, error=function(e) 75.0)
  ),
  adoption    = process_flow("adoption"),
  abandonment = process_flow("abandonment")
)

message(">>> Writing JSON...")
json_str <- toJSON(full_data, auto_unbox=TRUE, pretty=FALSE,
                   null="null", na="null")
writeLines(json_str, out_json, useBytes=TRUE)
message(sprintf("  Saved: %s (%.2f MB)", out_json, file.size(out_json)/1e6))
message("\n>>> generate_triads_json.R complete.")
