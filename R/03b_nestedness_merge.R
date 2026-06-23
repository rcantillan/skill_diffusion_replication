# ==============================================================================
# 03b_nestedness_merge.R
#
# Merges all_events_adoption_enriched.rds + skill_cs_scores.rds
# Final output: data/derived/riskset_adoption.rds
#
# Deletes all_events_adoption_enriched.rds on completion.
# Next: 03c_nestedness_merge_ab.R
# ==============================================================================

gc()
library(data.table)

out_dir     <- "data/derived"
input_file  <- file.path(out_dir, "all_events_adoption_enriched.rds")
cs_file     <- file.path(out_dir, "skill_cs_scores.rds")
output_file <- file.path(out_dir, "riskset_adoption.rds")

stopifnot("Input not found" = file.exists(input_file))
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
message(sprintf("  cs scores: %d skills | NAs: %d",
                nrow(cs), sum(is.na(cs$cs))))

# Merge cs and nested into dataset
setkey(cs, skill_name)
dt[cs, on = .(skill_name), `:=`(cs = i.cs, nested = i.nested)]

# Check coverage
n_na_cs <- sum(is.na(dt$cs))
pct_cs  <- 100 * mean(!is.na(dt$cs))
message(sprintf("  cs coverage: %.1f%% | NAs: %d", pct_cs, n_na_cs))
stopifnot("cs coverage < 99%" = pct_cs >= 99)
message("  [OK] cs coverage >= 99%")

# Classify ATC archetypes
# SC_Scaffolding   = Cognitive, cs >= median(cs | Cognitive)  [deeply nested]
# SC_Specialized   = Cognitive, cs <  median(cs | Cognitive)  [modular]
# Physical_Terminal = Physical (all cs levels)
med_cog <- median(dt[domain == "Cognitive", cs], na.rm = TRUE)
message(sprintf("  Median cs (Cognitive): %.3f", med_cog))

dt[, atc_archetype := fcase(
  domain == "Cognitive" & cs >= med_cog, "SC_Scaffolding",
  domain == "Cognitive" & cs <  med_cog, "SC_Specialized",
  domain == "Physical",                   "Physical_Terminal"
)]
message("  Archetype distribution:")
print(dt[, .(n_dyads = .N), by = atc_archetype][order(atc_archetype)])

# Final report
message(sprintf("\n  Final columns: %s",
                paste(sort(names(dt)), collapse = ", ")))
message(sprintf("  Rows: %s", format(nrow(dt), big.mark = ",")))

saveRDS(dt, output_file)
message(sprintf("  Saved: %s", output_file))
message(sprintf("  Size:  %.1f MB", file.size(output_file) / 1e6))

# Remove intermediate file
file.remove(input_file)
message("  [OK] Removed: all_events_adoption_enriched.rds")

rm(dt, cs); gc()
message("\n>>> 03b_nestedness_merge.R complete.")
message("    Next: 03c_nestedness_merge_ab.R")
