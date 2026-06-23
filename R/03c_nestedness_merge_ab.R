# ==============================================================================
# 03c_nestedness_merge_ab.R
#
# Merges all_events_abandonment_enriched.rds + skill_cs_scores.rds
# Final output: data/derived/riskset_abandonment.rds
#
# Deletes all_events_abandonment_enriched.rds on completion.
# Next: 04_status_pca.R
# ==============================================================================

gc()
library(data.table)

out_dir     <- "data/derived"
input_file  <- file.path(out_dir, "all_events_abandonment_enriched.rds")
cs_file     <- file.path(out_dir, "skill_cs_scores.rds")
output_file <- file.path(out_dir, "riskset_abandonment.rds")

stopifnot("Input not found — run 02b_enrich_abandonment.R first" =
            file.exists(input_file))
stopifnot("skill_cs_scores.rds not found — run 03_nestedness.R first" =
            file.exists(cs_file))

message(">>> Loading data...")
dt <- readRDS(input_file)
setDT(dt)
message(sprintf("  Rows: %s | Skills: %d",
                format(nrow(dt), big.mark = ","),
                uniqueN(dt$skill_name)))

cs <- readRDS(cs_file)
setDT(cs)

# Merge cs and nested
setkey(cs, skill_name)
dt[cs, on = .(skill_name), `:=`(cs = i.cs, nested = i.nested)]

n_na_cs <- sum(is.na(dt$cs))
pct_cs  <- 100 * mean(!is.na(dt$cs))
message(sprintf("  cs coverage: %.1f%% | NAs: %d", pct_cs, n_na_cs))
stopifnot("cs coverage < 99%" = pct_cs >= 99)
message("  [OK] cs coverage >= 99%")

# Archetypes — same classification as adoption
# Median computed from the same cs_file
med_cog <- median(dt[domain == "Cognitive", cs], na.rm = TRUE)
message(sprintf("  Median cs (Cognitive): %.3f", med_cog))

dt[, atc_archetype := fcase(
  domain == "Cognitive" & cs >= med_cog, "SC_Scaffolding",
  domain == "Cognitive" & cs <  med_cog, "SC_Specialized",
  domain == "Physical",                   "Physical_Terminal"
)]
message("  Archetype distribution:")
print(dt[, .(n_dyads = .N), by = atc_archetype][order(atc_archetype)])

message(sprintf("\n  Rows: %s", format(nrow(dt), big.mark = ",")))

saveRDS(dt, output_file)
message(sprintf("  Saved: %s", output_file))
message(sprintf("  Size:  %.1f MB", file.size(output_file) / 1e6))

# Remove intermediate file
file.remove(input_file)
message("  [OK] Removed: all_events_abandonment_enriched.rds")

rm(dt, cs); gc()
message("\n>>> 03c_nestedness_merge_ab.R complete.")
message("    Next: 04_status_pca.R")
