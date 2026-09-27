# 11_IC_enrichment_robustness.R
# Enrichment of stress-coupled regulators among the top-ranked core regulators
# (Dataset S6-5) and robustness of the IC ranking (Dataset S6-4).
#   S6-4  leave-one-measure-out, correlation of each measure with IC, dominance
#         after max-normalization, alternative normalizations
#   S6-5  hypergeometric test for the top 5 to top 20, a cutoff-free rank test,
#         and sensitivity to the classification of single regulators
#
# Inputs : outputs/tables/centralities_core.csv   (05_centralities.R)
#          curation/core_regulator_curation.csv
# Outputs: outputs/tables/S6-4_IC_robustness.csv
#          outputs/tables/S6-5_IC_enrichment.csv
#          outputs/tables/core_top20_IC.csv

source("scripts/common.R")
TOPN_MAIN  <- 15                     # cutoff used for the headline test
TOPN_RANGE <- 5:20                   # cutoff sensitivity

d <- read.csv(o("centralities_core.csv"), stringsAsFactors = FALSE)
s6 <- core_curation()

## Stress-coupled flag: Stress_category is set. Held in its own vector because
## "stress" is also the name of a centrality column in d.
scat_v  <- as.character(s6$Stress_category); scat_v[is.na(scat_v)] <- ""
is_stress <- setNames(trimws(scat_v) != "", s6$locus_tag)
d$stress_coupled <- unname(is_stress[d$locus_tag])
d$stress_coupled[is.na(d$stress_coupled)] <- FALSE
stopifnot(is.numeric(d$stress))            # the centrality column must survive
N <- nrow(d); K <- sum(d$stress_coupled)
cat("pool ", N, " regulators, ", K, " stress-coupled (", round(100*K/N), "%)\n", sep = "")
cat("check: stress centrality max = ", max(d$stress), "\n", sep = "")

maxnorm  <- function(x) { x[!is.finite(x)] <- 0; m <- max(x); if (m > 0) x / m else x * 0 }
ranknorm <- function(x) rank(x, ties.method = "average") / length(x)
zscore   <- function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else x * 0
minmax   <- function(x) { r <- range(x); if (diff(r) > 0) (x - r[1]) / diff(r) else x * 0 }
skewness <- function(x) { m <- mean(x); s <- sd(x); if (s == 0) 0 else mean((x - m)^3) / s^3 }
hyper    <- function(k, n) phyper(k - 1, K, N - K, n, lower.tail = FALSE)
nstress  <- function(tp) sum(is_stress[tp], na.rm = TRUE)
topset   <- function(sc, n) d$locus_tag[order(-sc)][seq_len(n)]

ic_full  <- rowSums(sapply(MEASURES, function(m) maxnorm(d[[m]])))
top_full <- topset(ic_full, TOPN_MAIN)

## ---------------------------------------------------- S6-4, four sections
rows <- list()
add  <- function(...) rows[[length(rows) + 1]] <<- c(list(...), rep(list(""), 8))[1:8]

add("Robustness of the integrated centrality score")
add(paste0("IC = sum over ", length(MEASURES),
           " measures of (measure / max measure): ", paste(MEASURES, collapse = ", ")))
add(sprintf("Background: %d of %d regulators are stress-coupled (%.0f%%)", K, N, 100*K/N))
add("")
add("1) Leave-one-measure-out")
add("measure dropped", "Spearman vs full IC", paste0("top-", TOPN_MAIN, " retained"),
    paste0("stress in top ", TOPN_MAIN), "hypergeometric p")
add("(none)", 1, paste0(TOPN_MAIN, "/", TOPN_MAIN),
    paste0(nstress(top_full), "/", TOPN_MAIN),
    signif(hyper(nstress(top_full), TOPN_MAIN), 3))
for (m in MEASURES) {
  ic <- rowSums(sapply(setdiff(MEASURES, m), function(x) maxnorm(d[[x]])))
  tp <- topset(ic, TOPN_MAIN); k <- nstress(tp)
  add(m, round(cor(ic_full, ic, method = "spearman"), 3),
      paste0(length(intersect(tp, top_full)), "/", TOPN_MAIN),
      paste0(k, "/", TOPN_MAIN), signif(hyper(k, TOPN_MAIN), 3))
}
add("")
add("2) Correlation of each measure with the full IC")
add("measure", "rho vs IC")
for (m in MEASURES)
  add(m, round(cor(maxnorm(d[[m]]), ic_full, method = "spearman"), 3))
add("")
add("3) Outlier and dominance after max-normalization")
add("measure", "max / 2nd-max", "top-node share of sum", "skewness")
for (m in MEASURES) {
  v <- sort(maxnorm(d[[m]]), decreasing = TRUE)
  add(m, ifelse(v[2] > 0, round(v[1] / v[2], 2), NA),
      round(v[1] / sum(v), 3), round(skewness(v), 2))
}
add("")
add("4) Ranking stability under alternative normalization")
add("normalization", "Spearman vs max-norm IC", paste0("top-", TOPN_MAIN, " retained"),
    paste0("stress in top ", TOPN_MAIN))
for (nmz in c("rank-sum", "z-score", "min-max")) {
  fn <- switch(nmz, "rank-sum" = ranknorm, "z-score" = zscore, "min-max" = minmax)
  ic <- rowSums(sapply(MEASURES, function(x) fn(d[[x]])))
  tp <- topset(ic, TOPN_MAIN)
  add(nmz, round(cor(ic_full, ic, method = "spearman"), 3),
      paste0(length(intersect(tp, top_full)), "/", TOPN_MAIN),
      paste0(nstress(tp), "/", TOPN_MAIN))
}
write.csv(do.call(rbind, lapply(rows, function(r)
  setNames(as.data.frame(r, stringsAsFactors = FALSE), paste0("V", 1:8)))),
  o("S6-4_IC_robustness.csv"), row.names = FALSE)

## ---------------------------------------------------- S6-5, enrichment
rows <- list()
add("Enrichment of stress-coupled regulators among top-ranked regulators")
add("Hypergeometric test, P(X >= k), drawing n regulators without replacement")
add(sprintf("Pool: %d regulators, %d stress-coupled", N, K))
add("")
add("Cutoff sensitivity")
add("top n", "stress-coupled in top n", "expected", "fold", "p", "significant at 0.05")
for (n in TOPN_RANGE) {
  tp <- topset(ic_full, n); k <- nstress(tp)
  p  <- hyper(k, n)
  add(n, k, round(n * K / N, 1), round((k / n) / (K / N), 2), signif(p, 3),
      ifelse(p < 0.05, "yes", "no"))
}
add("")
## Cutoff-free test: IC of all stress-coupled regulators against all others,
## Mann-Whitney, two-sided.
add("Cutoff-free rank test")
add("Mann-Whitney (Wilcoxon rank sum), two-sided, IC of stress-coupled vs other regulators")
ic_s <- ic_full[d$stress_coupled]
ic_o <- ic_full[!d$stress_coupled]
wt   <- suppressWarnings(wilcox.test(ic_s, ic_o))
add("group", "n", "median IC", "", "", "")
add("stress-coupled", length(ic_s), round(median(ic_s), 2), "", "", "")
add("other",          length(ic_o), round(median(ic_o), 2), "", "", "")
add("U", "p", "", "", "", "")
add(unname(wt$statistic), signif(wt$p.value, 3), "", "", "", "")
cat(sprintf("cutoff-free rank test: median IC %.2f vs %.2f, U = %g, p = %.3g\n",
            median(ic_s), median(ic_o), unname(wt$statistic), wt$p.value))
add("")
add("Sensitivity to the classification of individual regulators")
add("The test repeated with the stress-coupled pool increased or decreased by one regulator.")
add("pool K", paste0("stress in top ", TOPN_MAIN), "p", "significant at 0.05")
k0 <- nstress(top_full)
for (dk in c(-1, 0, 1)) {
  Kx <- K + dk; kx <- max(0, min(k0 + dk, TOPN_MAIN))
  p  <- phyper(kx - 1, Kx, N - Kx, TOPN_MAIN, lower.tail = FALSE)
  add(Kx, kx, signif(p, 3), ifelse(p < 0.05, "yes", "no"))
}
write.csv(do.call(rbind, lapply(rows, function(r)
  setNames(as.data.frame(r, stringsAsFactors = FALSE), paste0("V", 1:8)))),
  o("S6-5_IC_enrichment.csv"), row.names = FALSE)

## ---------------------------------------------------- core top-20 ranking
tf <- setNames(as.character(s6$TF_name), s6$locus_tag)
ord <- order(-ic_full)
top20 <- data.frame(
  rank        = seq_len(min(20, nrow(d))),
  locus_tag   = d$locus_tag[ord][1:20],
  TF_name     = unname(tf[d$locus_tag[ord][1:20]]),
  stress_coupled = ifelse(d$stress_coupled[ord][1:20], "yes", "no"),
  stress_category = unname(scat_v[match(d$locus_tag[ord][1:20], s6$locus_tag)]),
  IC          = round(ic_full[ord][1:20], 3),
  IC_seedmean = round(d$IC_seedmean[ord][1:20], 3),
  IC_seedsd   = round(d$IC_seedsd[ord][1:20], 3),
  stringsAsFactors = FALSE)
for (m in MEASURES) top20[[m]] <- round(d[[m]][ord][1:20], 4)
write.csv(top20, o("core_top20_IC.csv"), row.names = FALSE)
cat("\ncore top 20 by IC:\n")
print(top20[, c("rank","TF_name","stress_coupled","IC","IC_seedmean","IC_seedsd")], row.names = FALSE)

