# =============================================================================
# 99_paths_local.R
# Local paths — in .gitignore, never committed to the repository
# =============================================================================

# Canonical risk sets (Layer 2+)
PATH_ADOPTION    <- "data/derived/riskset_adoption.rds"
PATH_ABANDONMENT <- "data/derived/riskset_abandonment.rds"

# Raw O*NET data (Layer 1 — risk set construction)
PATH_ONET_2015 <- "data/raw/onet/db_15_1"
PATH_ONET_2024 <- "data/raw/onet/db_29_2_text"

# BLS OEWS wages
PATH_BLS_2015  <- "data/raw/bls/national_M2015_dl.xlsx"

# ACS education
PATH_ACS       <- "data/raw/acs"

# Crosswalk
PATH_CROSSWALK <- "data/crosswalk/onet_soc_crosswalk.xlsx"
