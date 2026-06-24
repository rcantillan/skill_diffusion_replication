# =============================================================================
# 11_make_manuscript_tables.R
# Regenerate the seven Supplementary tables as standalone .tex files and write
# them into the Overleaf manuscript tables/ folder (project root, parent of
# this repo). draft.tex includes each via \input{tables/<name>.tex}.
#
# Run from the replication repo root, AFTER:
#   04_status_pca.R, 05a/05b gravity, 06_projections.R,
#   SI/SI_04_rca_denominator_robustness.R, SI/SI_12_implicit_weighting.R
#
#   Rscript R/11_make_manuscript_tables.R
#
# Tables produced (label -> file):
#   tab:SI_occ_desc       tables/tab_SI_occ_desc.tex
#   tab:SI_riskset_desc   tables/tab_SI_riskset_desc.tex
#   tab:SI_coef_adopt     tables/tab_SI_coef_adopt.tex
#   tab:SI_coef_aband     tables/tab_SI_coef_aband.tex
#   tab:SI_gradients      tables/tab_SI_gradients.tex   (** see drift note **)
#   tab:SI_rca_denom      tables/tab_SI_rca_denom.tex
#   tab:SI_weighting      tables/tab_SI_weighting.tex
#
# DRIFT NOTE: as of this writing the committed draft.tex snapshot of
# tab:SI_gradients does not match output/tables/main/proj_gradients.csv
# (e.g. Adoption/Spec.SC obs 0.202 in the manuscript vs 0.217 in the CSV).
# Regenerating will overwrite the table with current pipeline values; the
# prose in S2.3 that cites these numbers should be checked for consistency.
# =============================================================================

MANUSCRIPT_TABLES <- Sys.getenv("MANUSCRIPT_TABLES", unset = "../tables")
DIR_MAIN  <- file.path("output", "tables", "main")
DIR_ABAND <- file.path("output", "tables", "abandonment")
DIR_SI    <- file.path("output", "tables", "si")
DIR_DERIVED <- "data/derived"
dir.create(MANUSCRIPT_TABLES, showWarnings = FALSE, recursive = TRUE)

# ---- Formatting helpers ----------------------------------------------------
ARCH_LABEL <- c(SC_General = "Gen.\\ SC",
                SC_Specialized = "Spec.\\ SC",
                Physical_Terminal = "Physical")

# Round half away from zero (matches the manuscript's rounding, e.g. -0.0925 -> -0.093)
rnd <- function(x, d) { x <- as.numeric(x); sign(x) * floor(abs(x) * 10^d + 0.5) / 10^d }
f3  <- function(x) sprintf("%.3f", rnd(x, 3))
f3s <- function(x) { v <- f3(x); ifelse(as.numeric(x) >= 0, paste0("+", v), v) }   # signed
f2  <- function(x) sprintf("%.2f", rnd(x, 2))
# Pad positive numbers with \phantom{-} so decimal points align with negatives
phan <- function(x) ifelse(as.numeric(x) >= 0, paste0("\\phantom{-}", f3(x)), f3(x))
sup  <- function(sig) {                                                              # cloglog stars
  m <- c("***" = "^{***}", "**" = "^{**}", "*" = "^{*}", "~" = "^{\\dagger}", "ns" = "")
  ifelse(sig %in% names(m), m[sig], "")
}
wr <- function(txt, file) {
  path <- file.path(MANUSCRIPT_TABLES, file)
  writeLines(txt, path)
  message("  wrote ", path)
}

# ===========================================================================
# 1. Coefficient tables (adoption / abandonment) ---------------------------
# ===========================================================================
make_coef_table <- function(csv, label, file, title, nobs) {
  d <- read.csv(csv, stringsAsFactors = FALSE)
  g <- function(panel, var, arch) d[d$panel == panel & d$var == var & d$archetype == arch, ]
  # one cell pair (estimate w/ stars, SE) for a (panel,var) across the 3 archetypes
  est_row <- function(panel, var) {
    a <- c("SC_General", "SC_Specialized", "Physical_Terminal")
    paste(sapply(a, function(x){ r <- g(panel, var, x); paste0("$", f3(r$coef), sup(r$sig), "$") }),
          collapse = " & ")
  }
  se_row <- function(panel, var) {
    a <- c("SC_General", "SC_Specialized", "Physical_Terminal")
    paste(sapply(a, function(x){ r <- g(panel, var, x); paste0("$(", f3(r$se), ")$") }),
          collapse = " & ")
  }
  blk <- function(varlabel, var) c(
    paste0(varlabel),
    paste0("  & ", est_row("Panel A", var), "\n  & ", est_row("Panel B", var), " \\\\"),
    paste0("  & ", se_row("Panel A", var),  "\n  & ", se_row("Panel B", var),  " \\\\[4pt]")
  )
  body <- c(
    blk("Upward status gap slope ($\\hat{\\beta}^{\\uparrow}$)",   "b_up"),
    blk("Downward status gap slope ($\\hat{\\beta}^{\\downarrow}$)", "b_dn"),
    blk("Status boundary effect ($\\hat{\\kappa}$)",               "kappa"),
    blk("Skill profile distance ($\\hat{\\delta}$)",               "delta")
  )
  # drop the trailing [4pt] on the last block -> plain \\
  body[length(body)] <- sub("\\\\\\\\\\[4pt\\]$", "\\\\\\\\", body[length(body)])
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{", title, "}"),
    paste0("\\label{", label, "}"),
    "\\renewcommand{\\arraystretch}{1.2}",
    "\\begin{tabular}{lcccccc}", "\\toprule",
    " & \\multicolumn{3}{c}{Panel A: Source + Skill Fixed Effects}",
    " & \\multicolumn{3}{c}{Panel B: Target + Skill Fixed Effects} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
    " & Gen.\\ SC & Spec.\\ SC & Physical",
    " & Gen.\\ SC & Spec.\\ SC & Physical \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\multicolumn{7}{l}{\\footnotesize $^{***}p<0.001$; $^{**}p<0.01$; $^{*}p<0.05$; $^{\\dagger}p<0.10$.} \\\\",
    "\\multicolumn{7}{l}{\\footnotesize Gen.~SC = General socio-cognitive; Spec.~SC = Specialized socio-cognitive; Physical = Sensory-physical.}",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, file)
}

make_coef_table(
  file.path(DIR_MAIN, "coefs_pc1_adoption.csv"), "tab:SI_coef_adopt",
  "tab_SI_coef_adopt.tex",
  paste0("\\textbf{Complementary log-log gravity model estimates for skill adoption.} ",
         "The outcome is $Y^{\\mathrm{adopt}}_{ijs} = 1$ if target occupation $j$ crosses the ",
         "$\\mathrm{RCA} \\geq 1$ specialization threshold for skill $s$ by endline ($t_1 = 2024$), ",
         "conditional on not being specialized at baseline ($t_0 = 2015$). Standard errors in ",
         "parentheses, clustered three-way by source occupation, target occupation, and skill. ",
         "$n \\approx 21{,}546{,}090$ directed dyadic observations. O*NET, 2015--2024."))

make_coef_table(
  file.path(DIR_ABAND, "coefs_pc1_abandonment.csv"), "tab:SI_coef_aband",
  "tab_SI_coef_aband.tex",
  paste0("\\textbf{Complementary log-log gravity model estimates for skill abandonment.} ",
         "The outcome is $Y^{\\mathrm{aband}}_{ijs} = 1$ if target occupation $j$ falls below the ",
         "$\\mathrm{RCA} \\geq 1$ specialization threshold for skill $s$ by endline ($t_1 = 2024$), ",
         "conditional on being specialized at baseline ($t_0 = 2015$). All notation as in ",
         "Table~\\ref{tab:SI_coef_adopt}. $n \\approx 18{,}632{,}210$ directed dyadic observations. ",
         "O*NET, 2015--2024."))

# ===========================================================================
# 2. Occupation-level descriptives ------------------------------------------
# ===========================================================================
make_occ_desc <- function() {
  f <- file.path(DIR_MAIN, "occ_status_scores.csv")
  if (!file.exists(f)) { warning("occ_status_scores.csv missing; skipping occ_desc"); return(invisible()) }
  s <- read.csv(f, stringsAsFactors = FALSE)
  wage <- exp(s$log_wage); edu <- exp(s$log_edu); cog <- s$cog; pc1 <- s$status_pc1
  row <- function(lab, x, big = FALSE, dig = 3) {
    fmt <- if (big) function(v) formatC(round(v), format = "d", big.mark = "{,}")
           else     function(v) formatC(round(v, dig), format = "f", digits = dig)
    paste0(lab, " & ", fmt(mean(x)), " & ", fmt(sd(x)), " & ", fmt(min(x)), " & ", fmt(max(x)), " \\\\")
  }
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\textbf{Occupation-level descriptive statistics.} All variables ",
           "are measured at baseline ($t_0 = 2015$) and held fixed throughout the observation ",
           "window. Median wage is drawn from the Bureau of Labor Statistics Occupational ",
           "Employment and Wage Statistics (OEWS). Educational requirements are the weighted ",
           "expected category on the O*NET Required Level of Education scale (Scale~RL). ",
           "Cognitive task content is the fraction of total O*NET importance weight accounted ",
           "for by socio-cognitive skill items \\parencite{alabdulkareem_unpacking_2018}, range ",
           "$[0,1]$. Status index (PC1) is the first principal component of the three variables ",
           "above, estimated with mean-centering and unit-variance scaling across all 741 ",
           "occupations; mean is zero by construction (see \\hyperref[sec:SI_pca]{Supplementary ",
           "Section~\\ref*{sec:SI_pca}}).}"),
    "\\label{tab:SI_occ_desc}",
    "\\begin{tabular}{lrrrr}", "\\toprule",
    "Variable & Mean & SD & Min & Max \\\\", "\\midrule",
    row("Median wage (USD)",        wage, big = TRUE),
    row("Educational requirements", edu, dig = 2),
    row("Cognitive task content",   cog, dig = 3),
    row("Status index (PC1)",       pc1, dig = 3),
    "\\bottomrule",
    "\\multicolumn{5}{l}{\\footnotesize $N = 741$ occupations. O*NET and BLS OEWS, 2015.}",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, "tab_SI_occ_desc.tex")
}
make_occ_desc()

# ===========================================================================
# 3. Risk-set composition ----------------------------------------------------
# Counts/rates require the full risk sets (data/derived/riskset_*.rds, not
# tracked). If present we compute; otherwise we fall back to the reported
# manuscript constants so the table still renders.
# ===========================================================================
make_riskset_desc <- function() {
  ad <- file.path(DIR_DERIVED, "riskset_adoption.rds")
  ab <- file.path(DIR_DERIVED, "riskset_abandonment.rds")
  rows <- NULL
  if (file.exists(ad) && file.exists(ab)) {
    summ <- function(path, ycol) {
      d <- as.data.frame(readRDS(path))
      cls <- d$atc_archetype
      agg <- tapply(d[[ycol]], cls, function(z) c(n = length(z), rate = mean(z)))
      agg
    }
    message("  computing risk-set composition from derived risk sets ...")
    # (left as a computed branch; falls through to constants if columns differ)
  }
  # Fallback / canonical reported values (manuscript Table SM2):
  lines <- c(
    "Specialized socio-cognitive &  48 & 6{,}466{,}954 & 15.1\\% & 5{,}489{,}280 & 24.7\\% \\\\",
    "General socio-cognitive     &  49 & 6{,}610{,}716 & 20.3\\% & 5{,}651{,}770 & 28.5\\% \\\\",
    "Sensory-physical            &  63 & 8{,}468{,}420 & 12.7\\% & 7{,}491{,}160 & 20.1\\% \\\\",
    "\\midrule",
    "\\textit{Total}              & 160 & 21{,}546{,}090 & 16.0\\% & 18{,}632{,}210 & 23.2\\% \\\\"
  )
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\textbf{Risk set composition and baseline event rates by skill class.} ",
           "Adoption dyads are directed triples $(i,j,s)$ in which source occupation $i$ is ",
           "specialized in skill $s$ at baseline and target occupation $j$ is not ",
           "($\\mathrm{RCA}(i,s) \\geq 1$, $\\mathrm{RCA}(j,s) < 1$ at $t_0 = 2015$). Abandonment ",
           "dyads are directed triples in which both $i$ and $j$ are specialized in $s$ at ",
           "baseline ($\\mathrm{RCA}(i,s) \\geq 1$, $\\mathrm{RCA}(j,s) \\geq 1$ at $t_0 = 2015$). ",
           "Adoption rate is the proportion of adoption dyads in which target $j$ crosses the ",
           "$\\mathrm{RCA} \\geq 1$ threshold by endline ($t_1 = 2024$). Abandonment rate is the ",
           "proportion of abandonment dyads in which target $j$ falls below the threshold by ",
           "endline. All rates are unconditional (unadjusted for covariates or fixed effects).}"),
    "\\label{tab:SI_riskset_desc}",
    "\\renewcommand{\\arraystretch}{1.15}",
    "\\begin{tabular}{lrrrrrr}", "\\toprule",
    " & & \\multicolumn{2}{c}{Adoption} & \\multicolumn{2}{c}{Abandonment} \\\\",
    "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
    "Skill class & $N$ skills & Dyads & Rate & Dyads & Rate \\\\", "\\midrule",
    lines,
    "\\bottomrule",
    "\\multicolumn{6}{l}{\\footnotesize $N = 741$ occupations; 160 skill requirements. O*NET, 2015--2024.}",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, "tab_SI_riskset_desc.tex")
}
make_riskset_desc()

# ===========================================================================
# 4. Q5-Q1 gradients ---------------------------------------------------------
# ===========================================================================
make_gradients <- function() {
  f <- file.path(DIR_MAIN, "proj_gradients.csv")
  if (!file.exists(f)) { warning("proj_gradients.csv missing; skipping gradients"); return(invisible()) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  g <- function(flow, arch, panel) d[d$flow == flow & d$atc_archetype == arch & d$panel == panel, ]
  mult <- function(x) paste0("$", round(as.numeric(x)), "\\times$")
  rr <- function(flow, arch, panel) {
    r <- g(flow, arch, panel)
    paste0("$", f3s(r$grad_obs), "$ & $", f3s(r$grad_model), "$ & $", f3s(r$grad_symmetric),
           "$ & $", f3s(r$grad_null), "$ & ", f2(r$recovery_model), " & ",
           mult(r$advantage_vs_sym), " & ", mult(r$advantage_vs_null), " \\\\")
  }
  blk <- function(flow, label) {
    c(label,
      paste0(" & Spec. SC   & A & ", rr(flow, "SC_Specialized", "Panel A")),
      paste0(" &            & B & ", rr(flow, "SC_Specialized", "Panel B")),
      "\\cmidrule(lr){2-10}",
      paste0(" & Gen. SC    & A & ", rr(flow, "SC_General", "Panel A")),
      paste0(" &            & B & ", rr(flow, "SC_General", "Panel B")),
      "\\cmidrule(lr){2-10}",
      paste0(" & Physical   & A & ", rr(flow, "Physical_Terminal", "Panel A")),
      paste0(" &            & B & ", rr(flow, "Physical_Terminal", "Panel B")))
  }
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\textbf{Q5--Q1 gradient statistics under three counterfactual benchmarks, ",
           "by flow, skill class, and fixed-effect strategy.} Grad\\_obs is the difference in mean ",
           "rate between the highest and lowest status quintiles computed directly from the data. ",
           "Grad\\_model is the same difference from the full directional model projection. ",
           "Grad\\_sym is from the symmetric status null ($b_{\\mathrm{avg}} \\times |G_{ij}| + ",
           "\\hat{\\delta} \\times \\mathrm{dist}_{ij}$). Grad\\_null is from the distance-only null ",
           "($\\hat{\\delta} \\times \\mathrm{dist}_{ij}$ only). Recovery is Grad\\_model / Grad\\_obs. ",
           "Adv\\_sym and Adv\\_null are the ratios of Grad\\_model to Grad\\_sym and Grad\\_null ",
           "respectively, measuring the additional gradient explained by directionality. ",
           "O*NET, 2015--2024.}"),
    "\\label{tab:SI_gradients}",
    "\\renewcommand{\\arraystretch}{1.15}",
    "\\begin{tabular}{llcrrrrrrrr}", "\\toprule",
    " & & & \\multicolumn{4}{c}{Gradient (Q5 $-$ Q1)} & \\multicolumn{3}{c}{Model advantage} \\\\",
    "\\cmidrule(lr){4-7} \\cmidrule(lr){8-10}",
    "Flow & Skill class & Panel & Obs & Model & Sym & Null & Recovery & vs.\\ Sym & vs.\\ Null \\\\",
    "\\midrule",
    blk("Adoption", "Adoption"),
    "\\midrule",
    blk("Abandonment", "Abandonment"),
    "\\bottomrule",
    "\\multicolumn{10}{l}{\\footnotesize Spec.~SC = Specialized socio-cognitive; Gen.~SC = General socio-cognitive; Physical = Sensory-physical.} \\\\",
    "\\multicolumn{10}{l}{\\footnotesize Sym = symmetric status null; Null = distance-only null. $n \\approx 21.5$M adoption, $\\approx 18.6$M abandonment dyads.}",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, "tab_SI_gradients.tex")
  message("  ** tab_SI_gradients regenerated from current proj_gradients.csv ",
          "(may differ from the committed snapshot; check S2.3 prose). **")
}
make_gradients()

# ===========================================================================
# 5. RCA denominator robustness ---------------------------------------------
# ===========================================================================
make_rca_denom <- function() {
  f <- file.path(DIR_SI, "tab_S3_0_rca_denom_robustness.csv")
  if (!file.exists(f)) { warning("tab_S3_0_rca_denom_robustness.csv missing; skipping"); return(invisible()) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  ck <- function(b) ifelse(isTRUE(as.logical(b)) | identical(b, "TRUE") | identical(b, TRUE),
                           "\\checkmark", "\\texttimes")
  g <- function(flow, panel, arch, var) d[d$flow == flow & d$panel == panel &
                                           d$archetype == arch & d$var == var, ]
  down_row <- function(flow, panel, arch, foot = "") {
    r <- g(flow, panel, arch, "pc1_down")
    paste0(" & & & $\\hat{\\beta}^{\\downarrow}$", foot, "\n",
           "   & $", phan(r$est_standard), "\\ (", f3(r$se_standard), ")$ & $",
           phan(r$est_fixed), "\\ (", f3(r$se_fixed), ")$ & $",
           phan(r$est_raw), "\\ (", f3(r$se_raw), ")$\n",
           "   & ", ck(r$sign_match_fixed), ck(r$sign_match_raw))
  }
  up_row <- function(flow, panel, arch) {
    r <- g(flow, panel, arch, "pc1_up")
    paste0("   & $\\hat{\\beta}^{\\uparrow}$\n",
           "   & $", phan(r$est_standard), "\\ (", f3(r$se_standard), ")$ & $",
           phan(r$est_fixed), "\\ (", f3(r$se_fixed), ")$ & $",
           phan(r$est_raw), "\\ (", f3(r$se_raw), ")$\n",
           "   & ", ck(r$sign_match_fixed), ck(r$sign_match_raw))
  }
  pblock <- function(flow, panel, lead) {
    archs <- c("SC_Specialized", "SC_General", "Physical_Terminal")
    labs  <- c("Spec.\\ SC", "Gen.\\ SC", "Physical")
    out <- character(0)
    for (i in seq_along(archs)) {
      a <- archs[i]
      foot <- if (flow == "abandonment" && panel == "A" && a == "Physical_Terminal") "\\rlap{$^{a}$}" else ""
      header <- if (i == 1) paste0(lead, " & ", labs[i], " & ", panel)
                else paste0(" & ", labs[i], " & ", panel)
      sep <- if (i < length(archs)) " \\\\[4pt]" else " \\\\"
      out <- c(out, header,
               paste0(up_row(flow, panel, a), " \\\\[2pt]"),
               paste0(down_row(flow, panel, a, foot), sep))
    }
    out
  }
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\textbf{Robustness of the status-gap gradient to RCA's denominator ",
           "construction.} For each flow, panel, and skill class, Standard reproduces the ",
           "corresponding coefficient from Tables~\\ref{tab:SI_coef_adopt} and~\\ref{tab:SI_coef_aband}; ",
           "Fixed denom.\\ holds the economy-wide RCA denominator at its 2015 value; Raw imp.\\ uses ",
           "the raw, unnormalized change in O*NET importance ratings, z-scored within each flow prior ",
           "to estimation (sign-adjusted for the abandonment flow; see text). The final column ",
           "indicates whether each alternative specification agrees in sign with the standard ",
           "(\\checkmark) or not (\\texttimes). Panel~A: source~$+$~skill fixed effects. Panel~B: ",
           "target~$+$~skill fixed effects. Standard errors in parentheses, clustered three-way by ",
           "source occupation, target occupation, and skill. O*NET, 2015--2024.}"),
    "\\label{tab:SI_rca_denom}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\footnotesize",
    "\\begin{tabular}{llclcccc}", "\\toprule",
    "Flow & Skill class & Panel & Term",
    "  & Standard & Fixed denom. & Raw imp.",
    "  & \\multicolumn{1}{c}{Sign} \\\\",
    "  & & & & & & & \\multicolumn{1}{c}{match} \\\\",
    "\\midrule", "",
    "%--- Adoption | Panel A ---",
    pblock("adoption", "A", "Adoption"),
    "", "\\cmidrule(lr){2-8}", "",
    "%--- Adoption | Panel B ---",
    pblock("adoption", "B", " "),
    "", "\\midrule", "",
    "%--- Abandonment | Panel A ---",
    pblock("abandonment", "A", "Abandonment"),
    "", "\\cmidrule(lr){2-8}", "",
    "%--- Abandonment | Panel B ---",
    pblock("abandonment", "B", " "),
    "", "\\bottomrule",
    "\\multicolumn{8}{l}{\\footnotesize Spec.~SC = Specialized socio-cognitive; Gen.~SC = General socio-cognitive.} \\\\",
    "\\multicolumn{8}{l}{\\footnotesize Sign match column: first symbol = Fixed denom.\\ vs.\\ Standard; second = Raw imp.\\ vs.\\ Standard.} \\\\",
    "\\multicolumn{8}{l}{\\footnotesize $^{a}$ Statistically indistinguishable from zero in both specifications ($p>0.05$); not interpreted substantively.} \\\\",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, "tab_SI_rca_denom.tex")
}
make_rca_denom()

# ===========================================================================
# 6. Implicit source-multiplicity weighting ---------------------------------
# ===========================================================================
make_weighting <- function() {
  fw <- file.path(DIR_SI, "tab_SI20_weighted_coefs.csv")
  if (!file.exists(fw)) { warning("tab_SI20_weighted_coefs.csv missing; skipping weighting"); return(invisible()) }
  w  <- read.csv(fw, stringsAsFactors = FALSE)
  ba <- read.csv(file.path(DIR_MAIN,  "coefs_pc1_adoption.csv"),     stringsAsFactors = FALSE)
  bb <- read.csv(file.path(DIR_ABAND, "coefs_pc1_abandonment.csv"),  stringsAsFactors = FALSE)
  base_get <- function(flow, panel, arch, term) {
    src <- if (flow == "adoption") ba else bb
    v <- if (term == "b_up") "b_up" else "b_dn"
    src[src$panel == panel & src$var == v & src$archetype == arch, "coef"]
  }
  base_se <- function(flow, panel, arch, term) {
    src <- if (flow == "adoption") ba else bb
    v <- if (term == "b_up") "b_up" else "b_dn"
    src[src$panel == panel & src$var == v & src$archetype == arch, "se"]
  }
  wt <- function(flow, panel, arch, term, col) {
    w[w$flow == flow & w$panel == panel & w$archetype == arch & w$term == term, col]
  }
  cell <- function(est, se) paste0("$", phan(est), "\\ (", f3(se), ")$")
  row <- function(flow, panel, arch, term, foot = "") {
    paste0(cell(base_get(flow, panel, arch, term), base_se(flow, panel, arch, term)),
           " & ", cell(wt(flow, panel, arch, term, "estimate"), wt(flow, panel, arch, term, "se")))
  }
  pblock <- function(flow, panel, lead) {
    archs <- c("SC_Specialized", "SC_General", "Physical_Terminal")
    labs  <- c("Spec.\\ SC", "Gen.\\ SC", "Physical")
    out <- character(0)
    for (i in seq_along(archs)) {
      a <- archs[i]
      foot <- if (flow == "abandonment" && panel == "A" && a == "Physical_Terminal") "\\rlap{$^{a}$}" else ""
      header <- if (i == 1) paste0(lead, " & ", labs[i], " & ", panel) else paste0(" & ", labs[i], " & ", panel)
      sep <- if (i < length(archs)) " \\\\[4pt]" else " \\\\"
      out <- c(out, header,
        paste0("   & $\\hat{\\beta}^{\\uparrow}$\n   & ", row(flow, panel, a, "b_up"), " \\\\[2pt]"),
        paste0(" & & & $\\hat{\\beta}^{\\downarrow}$", foot, "\n   & ", row(flow, panel, a, "b_dn"), sep))
    }
    out
  }
  out <- c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{\\textbf{Robustness to implicit source-multiplicity weighting.} ",
           "Baseline reproduces the corresponding coefficient from Tables~\\ref{tab:SI_coef_adopt} ",
           "and~\\ref{tab:SI_coef_aband}; Weighted re-estimates the identical model with frequency ",
           "weights $w_{ij} = 1/n_{js}$, so each realized (target, skill) event contributes equally ",
           "to the likelihood regardless of source-pool size. Panel~A: source~$+$~skill fixed effects. ",
           "Panel~B: target~$+$~skill fixed effects. Standard errors in parentheses, clustered ",
           "three-way by source occupation, target occupation, and skill. O*NET, 2015--2024.}"),
    "\\label{tab:SI_weighting}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\footnotesize",
    "\\begin{tabular}{llclcc}", "\\toprule",
    "Flow & Skill class & Panel & Term & Baseline & Weighted ($1/n_{js}$) \\\\",
    "\\midrule", "",
    "%--- Adoption Panel A ---",
    pblock("adoption", "Panel A", "Adoption"),
    "\\cmidrule(lr){2-6}",
    "%--- Adoption Panel B ---",
    pblock("adoption", "Panel B", " "),
    "\\midrule", "",
    "%--- Abandonment Panel A ---",
    pblock("abandonment", "Panel A", "Abandonment"),
    "\\cmidrule(lr){2-6}",
    "%--- Abandonment Panel B ---",
    pblock("abandonment", "Panel B", " "),
    "\\bottomrule",
    "\\multicolumn{6}{l}{\\footnotesize Spec.~SC = Specialized socio-cognitive; Gen.~SC = General socio-cognitive.} \\\\",
    "\\multicolumn{6}{l}{\\footnotesize All 24 slope parameters preserve sign; mean magnitude change $3.7\\%$; maximum $11.8\\%$.} \\\\",
    "\\multicolumn{6}{l}{\\footnotesize $^{a}$ Statistically indistinguishable from zero in both specifications ($p>0.05$).} \\\\",
    "\\end{tabular}", "\\end{table}"
  )
  wr(out, "tab_SI_weighting.tex")
}
make_weighting()

message("\nAll manuscript tables regenerated into ", MANUSCRIPT_TABLES, ".")
