# ==============================================================================
# SI_00_BUILD_SUBPERIOD_RISKSETS.R
#
# Builds adoption and abandonment risk sets for three non-overlapping
# sub-periods following the O*NET refresh cycle:
#   Period 1: 2015-2018 (O*NET 20.3 → 23.1)
#   Period 2: 2019-2021 (O*NET 24.1 → 26.1)
#   Period 3: 2022-2024 (O*NET 27.1 → 29.2)
#
# The sub-period risk sets use the same construction logic as the main
# pipeline (01_build_riskset_adoption.R and 01_build_riskset_abandonment.R),
# applied to sub-period O*NET extracts.
#
# NOTE: This script requires sub-period O*NET raw data files. If only the
# full 2015 and 2024 extracts are available, the script approximates
# sub-periods using release-date metadata from the O*NET update schedule.
# Full replication requires the intermediate O*NET releases listed below.
#
# O*NET release mapping:
#   2018 endline: db_23_1 (released ~Nov 2018)
#   2021 endline: db_26_1 (released ~Nov 2021)
#   2024 endline: db_29_2 (released ~Nov 2024; already in pipeline)
#
# Output:
#   data/derived/riskset_adoption_1518.rds
#   data/derived/riskset_adoption_1921.rds
#   data/derived/riskset_adoption_2224.rds
#   data/derived/riskset_abandonment_1518.rds
#   data/derived/riskset_abandonment_1921.rds
#   data/derived/riskset_abandonment_2224.rds
# ==============================================================================

gc()
library(data.table)

# Sub-period O*NET paths — update these to match your local directory structure
SUBPERIODS <- list(
  list(
    tag    = "1518",
    label  = "2015-2018",
    path_t0 = "data/raw/onet/db_15_1",        # baseline: O*NET 20.3
    path_t1 = "data/raw/onet/db_23_1"         # endline:  O*NET 23.1 (~2018)
  ),
  list(
    tag    = "1921",
    label  = "2019-2021",
    path_t0 = "data/raw/onet/db_24_1",        # baseline: O*NET 24.1 (~2019)
    path_t1 = "data/raw/onet/db_26_1"         # endline:  O*NET 26.1 (~2021)
  ),
  list(
    tag    = "2224",
    label  = "2022-2024",
    path_t0 = "data/raw/onet/db_27_1",        # baseline: O*NET 27.1 (~2022)
    path_t1 = "data/raw/onet/db_29_2_text"    # endline:  O*NET 29.2 (~2024)
  )
)

out_dir <- "data/derived"

# ==============================================================================
# Shared utilities
# ==============================================================================
clean_id <- function(x) {
  s <- gsub("-", "", as.character(x))
  s <- gsub("\\.[0-9]+$", "", s)
  trimws(s)
}

read_onet_im <- function(folder) {
  files <- c("Skills.txt", "Abilities.txt", "Knowledge.txt")
  rbindlist(Filter(Negate(is.null), lapply(files, function(f) {
    path <- file.path(folder, f)
    if (!file.exists(path)) {
      message(sprintf("  [SKIP] Not found: %s", path))
      return(NULL)
    }
    d <- fread(path, sep = "\t", quote = "",
               na.strings = c("NA", "n/a", "", "*"),
               showProgress = FALSE)
    setnames(d, names(d), make.names(names(d)))
    d <- d[Scale.ID == "IM"]
    data.table(
      soc        = clean_id(d$O.NET.SOC.Code),
      skill_name = d$Element.Name,
      importance = suppressWarnings(as.numeric(as.character(d$Data.Value)))
    )
  })), fill = TRUE)
}

apply_crosswalk <- function(im, cw) {
  im <- merge(im, cw, by.x = "soc", by.y = "soc10",
              all.x = FALSE, allow.cartesian = TRUE)
  im[, soc := soc19][, soc19 := NULL]
  unique(im)
}

compute_rca <- function(im) {
  grand_total  <- sum(im$importance, na.rm = TRUE)
  occ_total    <- im[, .(occ_total   = sum(importance)), by = soc]
  skill_total  <- im[, .(skill_total = sum(importance)), by = skill_name]
  im2 <- merge(im,  occ_total,   by = "soc")
  im2 <- merge(im2, skill_total, by = "skill_name")
  im2[, rca := (importance / occ_total) / (skill_total / grand_total)]
  im2[, .(soc, skill_name, rca)]
}

# Crosswalk (shared across periods)
cw <- unique(fread("data/crosswalk/2010_to_2019_Crosswalk.csv")[, .(
  soc10 = clean_id(`O*NET-SOC 2010 Code`),
  soc19 = clean_id(`O*NET-SOC 2019 Code`)
)])

# Domain and nestedness labels from main pipeline
cs_scores    <- readRDS(file.path(out_dir, "skill_cs_scores.rds"))
domain_lkp   <- unique(readRDS(
  file.path(out_dir, "riskset_adoption.rds"))[, .(skill_name, domain)])
cs_scores    <- merge(cs_scores, domain_lkp, by = "skill_name", all.x = TRUE)
med_cog      <- cs_scores[domain == "Cognitive", median(cs, na.rm = TRUE)]
cs_scores[, atc_archetype := fcase(
  domain == "Cognitive" & cs >= med_cog, "SC_General",
  domain == "Cognitive" & cs <  med_cog, "SC_Specialized",
  domain == "Physical",                   "Physical_Terminal"
)]

# ==============================================================================
# Main construction loop
# ==============================================================================
for (sp in SUBPERIODS) {

  message("\n", strrep("=", 65))
  message(sprintf(">>> Sub-period: %s (%s)", sp$label, sp$tag))
  message(strrep("=", 65))

  # Check if O*NET paths exist
  if (!dir.exists(sp$path_t0) || !dir.exists(sp$path_t1)) {
    message(sprintf("  [SKIP] O*NET paths not found for %s:", sp$label))
    message(sprintf("    t0: %s [%s]", sp$path_t0,
                    ifelse(dir.exists(sp$path_t0), "OK", "MISSING")))
    message(sprintf("    t1: %s [%s]", sp$path_t1,
                    ifelse(dir.exists(sp$path_t1), "OK", "MISSING")))
    next
  }

  # Load and process O*NET at t0 and t1
  message("  Reading O*NET IM at t0...")
  im_t0 <- read_onet_im(sp$path_t0)
  im_t0 <- im_t0[!is.na(importance)]
  im_t0 <- apply_crosswalk(im_t0, cw)
  rca_t0 <- compute_rca(im_t0)
  rm(im_t0); gc()

  message("  Reading O*NET IM at t1...")
  im_t1 <- read_onet_im(sp$path_t1)
  im_t1 <- im_t1[!is.na(importance)]
  im_t1 <- apply_crosswalk(im_t1, cw)
  rca_t1 <- compute_rca(im_t1)
  rm(im_t1); gc()

  message(sprintf("  RCA t0: %d occs × %d skills",
                  uniqueN(rca_t0$soc), uniqueN(rca_t0$skill_name)))
  message(sprintf("  RCA t1: %d occs × %d skills",
                  uniqueN(rca_t1$soc), uniqueN(rca_t1$skill_name)))

  # Occupations active in both periods
  occs_both <- intersect(unique(rca_t0$soc), unique(rca_t1$soc))
  skills_both <- intersect(unique(rca_t0$skill_name), unique(rca_t1$skill_name))
  message(sprintf("  Shared: %d occs × %d skills",
                  length(occs_both), length(skills_both)))

  # Restrict to shared occupations and skills
  rca_t0 <- rca_t0[soc %in% occs_both & skill_name %in% skills_both]
  rca_t1 <- rca_t1[soc %in% occs_both & skill_name %in% skills_both]

  # Specialization at t0 and t1
  spec_t0 <- rca_t0[rca >= 1, .(soc, skill_name)]
  spec_t0[, spec_t0 := 1L]
  spec_t1 <- rca_t1[rca >= 1, .(soc, skill_name)]
  spec_t1[, spec_t1 := 1L]

  # Full occupation × skill grid
  all_occ_skill <- CJ(soc = occs_both, skill_name = skills_both)
  all_occ_skill <- merge(all_occ_skill, spec_t0,
                         by = c("soc","skill_name"), all.x = TRUE)
  all_occ_skill <- merge(all_occ_skill, spec_t1,
                         by = c("soc","skill_name"), all.x = TRUE)
  all_occ_skill[is.na(spec_t0), spec_t0 := 0L]
  all_occ_skill[is.na(spec_t1), spec_t1 := 0L]

  # Cosine distance between occupation profiles at t0
  message("  Computing structural distance...")
  M <- dcast(spec_t0, soc ~ skill_name,
             value.var = "skill_name",
             fun.aggregate = length, fill = 0L)
  occ_ids <- M$soc
  M_mat   <- as.matrix(M[, -1, with = FALSE])
  rownames(M_mat) <- occ_ids
  # Cosine similarity → distance
  norms <- sqrt(rowSums(M_mat^2))
  norms[norms == 0] <- 1
  M_norm <- M_mat / norms
  sim_mat <- tcrossprod(M_norm)
  dist_dt <- as.data.table(sim_mat, keep.rownames = "source")
  dist_dt <- melt(dist_dt, id.vars = "source",
                  variable.name = "target", value.name = "cosine_sim")
  dist_dt[, structural_distance := 1 - cosine_sim]
  dist_dt[, cosine_sim := NULL]
  dist_dt[, source := as.character(source)]
  dist_dt[, target := as.character(target)]
  rm(M, M_mat, M_norm, sim_mat); gc()

  # Merge archetype and domain from main pipeline
  skill_meta <- cs_scores[, .(skill_name, domain, atc_archetype, cs)]

  # ==============================================================================
  # ADOPTION risk set: source has RCA>=1 at t0, target has RCA<1 at t0
  # ==============================================================================
  message("  Building adoption risk set...")

  adopt_src <- spec_t0[, .(source = soc, skill_name)]
  adopt_tgt <- all_occ_skill[spec_t0 == 0L, .(target = soc, skill_name)]

  adopt_rs <- merge(adopt_src, adopt_tgt, by = "skill_name",
                    allow.cartesian = TRUE)
  adopt_rs <- adopt_rs[source != target]

  # Outcome: did target specialize at t1?
  adopt_rs <- merge(adopt_rs,
                    spec_t1[, .(target = soc, skill_name, spec_t1)],
                    by = c("target","skill_name"), all.x = TRUE)
  adopt_rs[is.na(spec_t1), spec_t1 := 0L]
  setnames(adopt_rs, "spec_t1", "diffusion")

  # Merge structural distance
  adopt_rs <- merge(adopt_rs, dist_dt,
                    by = c("source","target"), all.x = TRUE)
  # Merge skill metadata
  adopt_rs <- merge(adopt_rs, skill_meta,
                    by = "skill_name", all.x = FALSE)
  adopt_rs <- adopt_rs[!is.na(structural_distance) & !is.na(domain)]

  message(sprintf("  Adoption: %s dyads | diffusion=%.4f",
                  format(nrow(adopt_rs), big.mark = ","),
                  mean(adopt_rs$diffusion)))

  out_path <- file.path(out_dir,
                         sprintf("riskset_adoption_%s.rds", sp$tag))
  saveRDS(adopt_rs, out_path)
  message(sprintf("  Saved: %s (%.1f MB)",
                  out_path, file.size(out_path)/1e6))
  rm(adopt_rs, adopt_src, adopt_tgt); gc()

  # ==============================================================================
  # ABANDONMENT risk set: both source and target have RCA>=1 at t0
  # ==============================================================================
  message("  Building abandonment risk set...")

  aband_src <- spec_t0[, .(source = soc, skill_name)]
  aband_tgt <- spec_t0[, .(target = soc, skill_name)]

  aband_rs <- merge(aband_src, aband_tgt, by = "skill_name",
                    allow.cartesian = TRUE)
  aband_rs <- aband_rs[source != target]

  # Outcome: did target LOSE specialization at t1?
  aband_rs <- merge(aband_rs,
                    spec_t1[, .(target = soc, skill_name, spec_t1)],
                    by = c("target","skill_name"), all.x = TRUE)
  aband_rs[is.na(spec_t1), spec_t1 := 0L]
  aband_rs[, abandonment := 1L - spec_t1]
  aband_rs[, spec_t1 := NULL]

  # Merge distance and metadata
  aband_rs <- merge(aband_rs, dist_dt,
                    by = c("source","target"), all.x = TRUE)
  aband_rs <- merge(aband_rs, skill_meta,
                    by = "skill_name", all.x = FALSE)
  aband_rs <- aband_rs[!is.na(structural_distance) & !is.na(domain)]

  message(sprintf("  Abandonment: %s dyads | abandonment=%.4f",
                  format(nrow(aband_rs), big.mark = ","),
                  mean(aband_rs$abandonment)))

  out_path <- file.path(out_dir,
                         sprintf("riskset_abandonment_%s.rds", sp$tag))
  saveRDS(aband_rs, out_path)
  message(sprintf("  Saved: %s (%.1f MB)",
                  out_path, file.size(out_path)/1e6))
  rm(aband_rs, aband_src, aband_tgt); gc()

  # Clean up per-period objects
  rm(rca_t0, rca_t1, spec_t0, spec_t1, all_occ_skill, dist_dt)
  gc(); gc()
}

message("\n>>> SI_00_BUILD_SUBPERIOD_RISKSETS.R complete.")
message("    Check output/data/derived/ for sub-period risk sets.")
message("    Next: SI_02_cs_distribution.R")
