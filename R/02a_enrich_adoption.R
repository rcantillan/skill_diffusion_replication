# ==============================================================================
# 02a_enrich_adoption.R
#
# Enriches the adoption risk set with dyadic covariates.
#
# Input:  data/derived/all_events_adoption.rds  (output of 01a_risk_set_adoption.R)
#         data/raw/onet/db_15_1/               (education, cog_score)
#         data/raw/bls/national_M2015_dl.xlsx  (BLS wages)
#         data/crosswalk/2010_to_2019_Crosswalk.csv
#
# Output: data/derived/all_events_adoption_enriched.rds
#         data/derived/cog_lkp.rds  (socio-cognitive fraction [0,1])
#
# Covariates added:
#   wage: s_wage, t_wage, log_s_wage, log_t_wage,
#         wage_gap, wage_up, wage_down, up_dummy
#   edu:  s_edu, t_edu, log_s_edu, log_t_edu,
#         edu_gap, edu_up, edu_down, edu_up_dummy
#   cog:  s_cog, t_cog, cog_gap, cog_up, cog_down, cog_up_dummy
#
# cog_score = socio-cognitive fraction (Alabdulkareem et al. 2018, eq. 3)
#   cog_j = sum_{s in Cognitive} onet(j,s) / sum_{s in S} onet(j,s)
#   Range [0,1] — direct measure independent of wage and education
#
# edu = Required Level of Education (O*NET Scale RL), weighted category
#
# Intermediates: deleted before final save
# Next:          03_nestedness.R
# ==============================================================================

gc()
library(data.table)
library(readxl)

# ==============================================================================
# Paths
# ==============================================================================
if (file.exists("R/99_paths_local.R")) source("R/99_paths_local.R")

path_2015   <- if (exists("PATH_ONET_2015")) PATH_ONET_2015 else
               "data/raw/onet/db_15_1"
path_bls    <- if (exists("PATH_BLS_2015")) PATH_BLS_2015 else
               "data/raw/bls/national_M2015_dl.xlsx"
path_cw     <- "data/crosswalk/2010_to_2019_Crosswalk.csv"
out_dir     <- "data/derived"
input_file  <- file.path(out_dir, "all_events_adoption.rds")

stopifnot("Input not found — run 01a_risk_set_adoption.R first" =
            file.exists(input_file))
stopifnot("O*NET 2015 not found" = dir.exists(path_2015))
stopifnot("BLS 2015 not found"   = file.exists(path_bls))
stopifnot("Crosswalk not found"  = file.exists(path_cw))
message("[OK] All inputs exist")

# ==============================================================================
# clean_id — identical to 01a_risk_set_adoption.R
# ==============================================================================
clean_id <- function(x) {
  s <- as.character(x)
  s <- gsub("-", "", s)
  s <- gsub("\\.[0-9]+$", "", s)
  trimws(s)
}
stopifnot(clean_id("15-1212.01") == "151212")
message("[OK] clean_id verified")

# ==============================================================================
# Step 1 — Load base risk set
# ==============================================================================
message("\n>>> Step 1: Loading adoption risk set...")

dt <- readRDS(input_file)
setDT(dt)

message(sprintf("  Rows:       %s", format(nrow(dt), big.mark = ",")))
message(sprintf("  Source occ: %d", uniqueN(dt$source)))
message(sprintf("  Target occ: %d", uniqueN(dt$target)))
message(sprintf("  Skills:     %d", uniqueN(dt$skill_name)))
message(sprintf("  Diffusion:  %.4f", mean(dt$diffusion)))

occ_in_risk_set <- union(unique(dt$source), unique(dt$target))

# Extract domain lookup from the dataset
domain_lkp <- unique(dt[, .(skill_name, domain)])
message("  Domains:")
print(domain_lkp[, .N, by = domain][order(-N)])

gc()

# ==============================================================================
# Step 2 — Wages (BLS OEWS 2015, A_MEDIAN)
# ==============================================================================
message("\n>>> Step 2: BLS 2015 wages...")

bls <- as.data.table(read_excel(path_bls, sheet = 1))
setnames(bls, names(bls), toupper(names(bls)))

# Keep only detailed level
if ("OCC_GROUP" %in% names(bls))
  bls <- bls[OCC_GROUP == "detailed"]

bls[, soc_code := clean_id(OCC_CODE)]
bls[, wage := suppressWarnings(as.numeric(A_MEDIAN))]
bls <- bls[!is.na(wage) & nchar(soc_code) == 6]

wages <- bls[, .(wage = mean(wage, na.rm = TRUE)), by = soc_code]
wages[, log_wage := log(wage)]
message(sprintf("  BLS occupations: %d", nrow(wages)))
message(sprintf("  Wage: median=$%.0f | range [$%.0f, $%.0f]",
                median(wages$wage), min(wages$wage), max(wages$wage)))
rm(bls); gc()

# Merge wage to dataset
setkey(wages, soc_code)
dt[wages, on = .(source = soc_code),
   `:=`(s_wage = i.wage, log_s_wage = i.log_wage)]
dt[wages, on = .(target = soc_code),
   `:=`(t_wage = i.wage, log_t_wage = i.log_wage)]
rm(wages); gc()

dt[, wage_gap  := log_t_wage - log_s_wage]
dt[, wage_up   := pmax(0, wage_gap)]
dt[, wage_down := pmin(0, wage_gap)]
dt[, up_dummy  := fifelse(!is.na(wage_gap) & wage_gap > 0, 1L, 0L)]

pct_wage <- 100 * mean(!is.na(dt$wage_gap))
message(sprintf("  wage_gap coverage: %.1f%%", pct_wage))
if (pct_wage < 80)
  warning("wage_gap coverage < 80% — check BLS-SOC crosswalk")

# ==============================================================================
# Step 3 — Education (O*NET 2015, Scale RL = Required Level of Education)
# ==============================================================================
message("\n>>> Step 3: Education O*NET 2015 (Scale RL)...")

edu_file <- file.path(path_2015,
                      "Education, Training, and Experience.txt")
stopifnot("O*NET education file not found" = file.exists(edu_file))

edu_raw <- fread(edu_file, sep = "\t", quote = "",
                 na.strings = c("NA", "n/a", "", "*"),
                 showProgress = FALSE)
setnames(edu_raw, names(edu_raw), make.names(names(edu_raw)))

edu_raw <- edu_raw[
  Element.Name == "Required Level of Education" &
  Scale.ID == "RL"
]
edu_raw[, soc_code := clean_id(O.NET.SOC.Code)]
edu_raw[, cat := suppressWarnings(as.numeric(Category))]
edu_raw[, val := suppressWarnings(as.numeric(as.character(Data.Value)))]
edu_raw <- edu_raw[!is.na(cat) & !is.na(val) & val > 0]
message(sprintf("  RL rows: %d | Occupations: %d",
                nrow(edu_raw), uniqueN(edu_raw$soc_code)))

# Crosswalk 2010 -> 2019
cw_e <- unique(fread(path_cw)[, .(
  soc10 = clean_id(`O*NET-SOC 2010 Code`),
  soc19 = clean_id(`O*NET-SOC 2019 Code`)
)])
edu_raw <- merge(edu_raw, cw_e,
                 by.x = "soc_code", by.y = "soc10",
                 all.x = FALSE, allow.cartesian = TRUE)
edu_raw[, soc_code := soc19]
rm(cw_e)

# Weighted expected category (weighted mean by response proportion)
edu <- edu_raw[, .(
  edu_expected = sum(cat * val, na.rm = TRUE) / sum(val, na.rm = TRUE)
), by = soc_code]
edu <- edu[!is.na(edu_expected)]
edu[, log_edu := log(edu_expected)]
rm(edu_raw); gc()

message(sprintf("  Occupations with edu: %d", nrow(edu)))
message(sprintf("  edu: median=%.2f | range [%.2f, %.2f]",
                median(edu$edu_expected),
                min(edu$edu_expected),
                max(edu$edu_expected)))

setkey(edu, soc_code)
dt[edu, on = .(source = soc_code),
   `:=`(s_edu = i.edu_expected, log_s_edu = i.log_edu)]
dt[edu, on = .(target = soc_code),
   `:=`(t_edu = i.edu_expected, log_t_edu = i.log_edu)]
rm(edu); gc()

dt[, edu_gap      := log_t_edu - log_s_edu]
dt[, edu_up       := pmax(0,  edu_gap)]
dt[, edu_down     := pmin(0,  edu_gap)]
dt[, edu_up_dummy := fifelse(!is.na(edu_gap) & edu_gap > 0, 1L, 0L)]

message(sprintf("  edu_gap coverage: %.1f%%",
                100 * mean(!is.na(dt$edu_gap))))

# ==============================================================================
# Step 4 — Cognitive score (Alabdulkareem et al. 2018, eq. 3)
#
# cog_j = sum_{s in Cognitive} onet(j,s) / sum_{s in S} onet(j,s)
# where onet(j,s) = (IM - 1) / 4  [normalized to [0,1]]
#
# Direct measure independent of wage and education -> clean composite PCA
# ==============================================================================
message("\n>>> Step 4: Cognitive score (Alabdulkareem eq. 3)...")

read_im <- function(folder, filename) {
  f <- file.path(folder, filename)
  if (!file.exists(f)) { message("  WARNING: not found — ", filename); return(NULL) }
  d <- fread(f, sep = "\t", quote = "",
             na.strings = c("NA", "n/a", "", "*"),
             showProgress = FALSE)
  setnames(d, names(d), make.names(names(d)))
  d <- d[Scale.ID == "IM"]
  d[, .(
    soc_code   = clean_id(O.NET.SOC.Code),
    skill_name = Element.Name,
    im_raw     = suppressWarnings(as.numeric(as.character(Data.Value)))
  )]
}

im <- rbindlist(Filter(Negate(is.null), list(
  read_im(path_2015, "Skills.txt"),
  read_im(path_2015, "Abilities.txt"),
  read_im(path_2015, "Knowledge.txt")
)), fill = TRUE)

im <- im[!is.na(im_raw) & soc_code %in% occ_in_risk_set]
im[, onet := (im_raw - 1) / 4]  # normalize to [0,1]
im[, im_raw := NULL]

message(sprintf("  IM rows: %s | Occupations: %d | Skills: %d",
                format(nrow(im), big.mark = ","),
                uniqueN(im$soc_code),
                uniqueN(im$skill_name)))

# Merge with domain (Cognitive / Physical)
im <- merge(im, domain_lkp, by = "skill_name", all.x = TRUE)
im <- im[!is.na(domain)]

# Detect cognitive domain label
cog_domain <- unique(im$domain)[grepl("cog", unique(im$domain),
                                      ignore.case = TRUE)][1]
message(sprintf("  Cognitive domain detected: '%s'", cog_domain))

# cog_score per occupation
cog_lkp <- im[, .(
  sum_cog   = sum(onet[domain == cog_domain], na.rm = TRUE),
  sum_total = sum(onet, na.rm = TRUE)
), by = soc_code]
cog_lkp[, cog_score := fifelse(sum_total > 0,
                                sum_cog / sum_total, NA_real_)]
cog_lkp <- cog_lkp[!is.na(cog_score), .(soc_code, cog_score)]
rm(im); gc()

stopifnot("cog_score outside [0,1]" =
            min(cog_lkp$cog_score) >= 0 &
            max(cog_lkp$cog_score) <= 1)
message(sprintf("  Occupations: %d | mean=%.4f | range [%.4f, %.4f]",
                nrow(cog_lkp),
                mean(cog_lkp$cog_score),
                min(cog_lkp$cog_score),
                max(cog_lkp$cog_score)))
message("  [OK] Range [0,1] verified")

# Save cog_lkp for use in 02b_enrich_abandonment.R and 04_status_pca.R
saveRDS(cog_lkp, file.path(out_dir, "cog_lkp.rds"))
message("  [OK] Saved: data/derived/cog_lkp.rds")

# Merge to dataset
setkey(cog_lkp, soc_code)
dt[cog_lkp, on = .(source = soc_code), s_cog := i.cog_score]
dt[cog_lkp, on = .(target = soc_code), t_cog := i.cog_score]
rm(cog_lkp); gc()

dt[, cog_gap      := t_cog - s_cog]
dt[, cog_up       := pmax(0,  cog_gap)]
dt[, cog_down     := pmin(0,  cog_gap)]
dt[, cog_up_dummy := fifelse(!is.na(cog_gap) & cog_gap > 0, 1L, 0L)]

message(sprintf("  cog_gap coverage: %.1f%%",
                100 * mean(!is.na(dt$cog_gap))))
message(sprintf("  cog_gap: sd=%.4f | range [%.4f, %.4f]",
                sd(dt$cog_gap, na.rm = TRUE),
                min(dt$cog_gap, na.rm = TRUE),
                max(dt$cog_gap, na.rm = TRUE)))

# Diagnostic: correlation cog_gap ~ wage_gap (expected > 0.25)
r_cog_wage <- cor(dt$cog_gap, dt$wage_gap, use = "complete.obs")
message(sprintf("  Correlation cog_gap ~ wage_gap: %.3f", r_cog_wage))
if (r_cog_wage < 0.25)
  warning("Correlation < 0.25 — composite PCA may be weak")

# ==============================================================================
# Step 5 — Coverage checks
# ==============================================================================
message("\n>>> Step 5: Coverage checks...")

check_cov <- function(col, label, thr = 0.85) {
  pct <- mean(!is.na(dt[[col]]))
  status <- ifelse(pct >= thr, "[OK]", "[CHECK]")
  message(sprintf("  %s %-25s %.1f%%", status, label, 100 * pct))
}
check_cov("wage_gap",            "wage_gap")
check_cov("edu_gap",             "edu_gap",  thr = 0.70)
check_cov("cog_gap",             "cog_gap")
check_cov("structural_distance", "structural_distance")
stopifnot(sum(is.na(dt$structural_distance)) == 0)

# ==============================================================================
# Step 6 — Preliminary DTC diagnostic
# ==============================================================================
message("\n>>> Step 6: Preliminary DTC diagnostic...")

dt[, direction := fcase(
  wage_gap > 0,     "upward",
  wage_gap < 0,     "downward",
  !is.na(wage_gap), "lateral"
)]
tab <- dt[!is.na(direction) & direction != "lateral",
          .(n = .N, rate = mean(diffusion)),
          by = .(domain, direction)]
setorder(tab, domain, direction)
print(tab)
wide <- dcast(tab, domain ~ direction, value.var = "rate")
if (all(c("upward", "downward") %in% names(wide))) {
  wide[, ratio := round(upward / downward, 3)]
  message("  Ratio upward/downward:")
  print(wide)
  message("  Expected — Cognitive: ratio > 1 | Physical: ratio < 1")
}
dt[, direction := NULL]

# ==============================================================================
# Step 7 — Save and clean up
# ==============================================================================
message("\n>>> Step 7: Saving...")

message(sprintf("  Rows:    %s", format(nrow(dt), big.mark = ",")))
message(sprintf("  Columns: %s", paste(sort(names(dt)), collapse = ", ")))
message(sprintf("  NAs wage_gap: %d | edu_gap: %d | cog_gap: %d",
                sum(is.na(dt$wage_gap)),
                sum(is.na(dt$edu_gap)),
                sum(is.na(dt$cog_gap))))

out_path <- file.path(out_dir, "all_events_adoption_enriched.rds")
saveRDS(dt, out_path)
message(sprintf("  Saved: %s", out_path))
message(sprintf("  Size:  %.1f MB", file.size(out_path) / 1e6))

# Remove intermediate file
file.remove(file.path(out_dir, "all_events_adoption.rds"))
message("  [OK] Removed: all_events_adoption.rds")

rm(dt); gc()
message("\n>>> 02a_enrich_adoption.R complete.")
message("    Next: 02b_enrich_abandonment.R")
