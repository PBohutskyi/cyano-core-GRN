# 15_crossvalidation.R
# Cross-validation of high-influence core regulators against the three species
# networks (Figure 5), and its sensitivity to the rank cutoffs.
# A regulator is cross-validated when it reaches the core top k on one measure
# and its ortholog reaches the species top n on the same measure. Measures:
# degree, k-core, betweenness, stress, eigenvector and IC; group labels are for
# reporting only. Ties at a cutoff (frequent for k-core) are broken by IC.
#   S6-6  cross-validated regulators and enrichment of stress-related ones
#   S6-7  sensitivity grid over core and species cutoffs
#
# Inputs : outputs/tables/centralities_{core,s7942,s6803,s7002}.csv  (05_centralities.R)
#          curation/S5_curated_regulators.xlsx   (S5-4 ortholog map)
#          curation/core_regulator_curation.csv
# Outputs: outputs/tables/S6-6_xval_enrichment.csv
#          outputs/tables/S6-7_xval_sensitivity.csv
#          outputs/tables/xval_matches.csv    (every regulator x species x measure match)

source("scripts/common.R")
CORE_K <- 12; SPECIES_N <- 10          # baseline cutoffs
GRID_K <- c(8, 10, 12, 15, 19); GRID_N <- c(5, 8, 10, 12, 15)

XVAL_MEASURES <- c(degree = "local (degree/k-core)",
                   kcore  = "local (degree/k-core)",
                   betweenness = "global (betweenness/stress)",
                   stress      = "global (betweenness/stress)",
                   eigenvector = "community (eigenvector)",
                   IC          = "integrated (IC)")

rd <- function(nm) read.csv(o(paste0("centralities_", nm, ".csv")), stringsAsFactors = FALSE)
core <- rd("core"); sp <- list("7942" = rd("s7942"), "6803" = rd("s6803"), "7002" = rd("s7002"))
ortho <- ortholog_map()

s6 <- core_curation()
sc <- as.character(s6$Stress_category); sc[is.na(sc)] <- ""
stress <- setNames(trimws(sc) != "", s6$locus_tag)
tfname <- setNames(as.character(s6$TF_name), s6$locus_tag)
scat   <- setNames(sc, s6$locus_tag)
N <- nrow(core); K <- sum(stress[core$locus_tag], na.rm = TRUE)

topby <- function(d, m, n) {
  x <- d[[m]]; x[!is.finite(x)] <- 0
  ## ties broken by IC
  ic <- d$IC; ic[!is.finite(ic)] <- 0
  d$locus_tag[order(-x, -ic)][seq_len(min(n, nrow(d)))]
}

## every (regulator, species, group) match at a given pair of cutoffs
instances <- function(k, n) {
  out <- list()
  for (m in names(XVAL_MEASURES)) {
    ck <- topby(core, m, k)
    for (spn in names(sp)) {
      st  <- topby(sp[[spn]], m, n)
      hit <- ck[!is.na(ortho[[spn]][ck]) & ortho[[spn]][ck] %in% st]
      for (h in hit) out[[length(out) + 1]] <-
        data.frame(regulator = unname(tfname[h]), locus_tag = h, species = spn,
                   measure = m, group = unname(XVAL_MEASURES[m]),
                   stringsAsFactors = FALSE)
    }
  }
  if (length(out)) { out <- do.call(rbind, out); return(unique(out)) }
  data.frame()
}

base <- instances(CORE_K, SPECIES_N)
write.csv(base, o("xval_matches.csv"), row.names = FALSE)
distinct <- unique(base$locus_tag)
k_stress <- sum(stress[distinct], na.rm = TRUE)
cat("baseline: core top ", CORE_K, " x species top ", SPECIES_N, " -> ",
    nrow(base), " matches, ", length(distinct), " distinct regulators, ",
    k_stress, " stress-related\n", sep = "")

## ---------------------------------------------------------------- S6-6
rows <- list(); add <- function(...) rows[[length(rows)+1]] <<- c(list(...), rep(list(""),7))[1:7]
add("Enrichment of stress-related regulators in the cross-validated set")
add(sprintf("Baseline: core top %d per measure, species top %d per measure",
            CORE_K, SPECIES_N))
add(sprintf("%d matches (regulator x species x measure), %d distinct regulators",
            nrow(base), length(distinct)))
add("")
add("Cross-validated regulators")
add("regulator", "locus tag", "species", "n species", "measure groups",
    "stress category", "stress-related")
for (g in distinct[order(-sapply(distinct, function(x)
      length(unique(base$species[base$locus_tag == x]))), unname(tfname[distinct]))]) {
  b <- base[base$locus_tag == g, ]
  add(unname(tfname[g]), g, paste(sort(unique(b$species)), collapse = ", "),
      length(unique(b$species)), paste(sort(unique(b$group)), collapse = "; "),
      unname(scat[g]), ifelse(isTRUE(stress[[g]]), "Y", "N"))
}
add("")
add(sprintf("Hypergeometric test: draw n=%d from a pool of %d with %d stress-related",
            length(distinct), N, K))
add("K stress in pool", "k stress in cross-validated set", "expected", "fold", "p", "significant")
p <- phyper(k_stress - 1, K, N - K, length(distinct), lower.tail = FALSE)
add(K, k_stress, round(length(distinct) * K / N, 1),
    round((k_stress / length(distinct)) / (K / N), 2), signif(p, 4),
    ifelse(p < 0.05, "yes", "no"))
add("")
add("Sensitivity to the classification of individual regulators")
for (dk in c(-1, 1)) {
  Kx <- K + dk; kx <- max(0, min(k_stress + dk, length(distinct)))
  px <- phyper(kx - 1, Kx, N - Kx, length(distinct), lower.tail = FALSE)
  add(Kx, kx, "", "", signif(px, 4), ifelse(px < 0.05, "yes", "no"))
}
write.csv(do.call(rbind, lapply(rows, function(r)
  setNames(as.data.frame(r, stringsAsFactors = FALSE), paste0("V", 1:7)))),
  o("S6-6_xval_enrichment.csv"), row.names = FALSE)

## ---------------------------------------------------------------- S6-7
rows <- list()
add("Threshold sensitivity of the cross-validation")
add("A regulator is cross-validated when it ranks in the core top-k on any one")
add("measure and its ortholog ranks in the species top-n on that same measure.")
add(sprintf("Baseline: core top %d, species top %d", CORE_K, SPECIES_N))
add(sprintf("Pool: %d regulators, %d stress-related", N, K))
add("")
add("A. Sensitivity grid")
add("core top-k", "species top-n", "matches", "distinct regulators",
    "retained of baseline", "stress-related", "hypergeometric p")
for (k in GRID_K) for (n in GRID_N) {
  x  <- instances(k, n); dd <- unique(x$locus_tag)
  ks <- sum(stress[dd], na.rm = TRUE)
  pp <- phyper(ks - 1, K, N - K, length(dd), lower.tail = FALSE)
  add(ifelse(k == CORE_K && n == SPECIES_N, paste0(k, " (baseline)"), k), n,
      nrow(x), length(dd),
      paste0(length(intersect(dd, distinct)), "/", length(distinct)),
      paste0(ks, "/", length(dd)), signif(pp, 3))
}
write.csv(do.call(rbind, lapply(rows, function(r)
  setNames(as.data.frame(r, stringsAsFactors = FALSE), paste0("V", 1:7)))),
  o("S6-7_xval_sensitivity.csv"), row.names = FALSE)
