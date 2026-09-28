# build_dataset_S8.R
# Dataset S8: reproducibility of the four networks across ten seeded GENIE3
# runs. Every seeded network is pruned exactly as its representative network, and
# every value is read from the tables behind Figure 4 and Supplementary
# Figure S2, so the dataset and the figures carry identical numbers.
#
# Inputs : outputs/tables/fig4_A_edge_overlap.csv, figS2_A_edge_overlap.csv
#          outputs/tables/seed_outdegree_core.csv, seed_outdegree_species.csv
#          outputs/tables/seed_centralities_core.csv, fig4_D_top15_recurrence.csv,
#          fig4_E_stress_by_seed.csv, fig4_B_outdegree.csv, centralities_core.csv
#          (scripts 05, 12, 13, 14)
#          outputs/tables/network_recurrence_summary.csv (script 14)
# Output : datasets/Dataset_S8_Seeded_run_stability.xlsx

source("scripts/common.R")
suppressPackageStartupMessages(library(openxlsx))

OUT  <- "datasets/Dataset_S8_Seeded_run_stability.xlsx"
TOPN <- 15
SHEETS <- c(legend = "S8-0. Legend",
            ov     = "S8-1. Edge overlap",
            od     = "S8-2. Regulator out-degree",
            cent   = "S8-3. Core per-seed centrality",
            rec    = "S8-4. Core top-15 recurrence",
            enr    = "S8-5. Core stress enrichment",
            nrec   = "S8-6. Network recurrence")
## Excel rejects sheet names longer than 31 characters
stopifnot(all(nchar(SHEETS) <= 31))

NET_LAB <- sapply(NETWORKS, function(n) n$label)

need <- function(f) { if (!file.exists(f)) stop("missing ", f); f }
for (f in c("fig4_A_edge_overlap.csv", "figS2_A_edge_overlap.csv",
            "seed_outdegree_core.csv", "seed_outdegree_species.csv",
            "fig4_B_outdegree.csv", "fig4_D_top15_recurrence.csv",
            "fig4_E_stress_by_seed.csv", "seed_centralities_core.csv",
            "centralities_core.csv", "network_recurrence_summary.csv")) need(o(f))

## ---------------------------------------------------------------------------
## S8-1. Edge overlap between runs (Figure 4A, Figure S2A)
## ---------------------------------------------------------------------------
A4 <- read.csv(o("fig4_A_edge_overlap.csv"), stringsAsFactors = FALSE)
A4$network <- "core"
AS <- read.csv(o("figS2_A_edge_overlap.csv"), stringsAsFactors = FALSE)
A  <- rbind(A4[, c("network", "group", "run_a", "run_b", "overlap")],
            AS[, c("network", "group", "run_a", "run_b", "overlap")])
ov <- data.frame(
  Network    = unname(NET_LAB[A$network]),
  Comparison = ifelse(A$group == "rep_vs_seed", "representative vs seeded", "seeded vs seeded"),
  Run_A      = A$run_a,
  Run_B      = A$run_b,
  Shared_edges_percent = round(A$overlap, 2),
  stringsAsFactors = FALSE)
run_num <- function(x) { v <- suppressWarnings(as.integer(sub("seed", "", x, fixed = TRUE))); v[is.na(v)] <- 0L; v }
ov <- ov[c_order(match(ov$Network, NET_LAB), ov$Comparison, run_num(ov$Run_A), run_num(ov$Run_B)), ]
stopifnot(all(table(ov$Network, ov$Comparison)[, "seeded vs seeded"] == choose(length(SEEDS), 2)))
med <- tapply(ov$Shared_edges_percent, paste(ov$Network, ov$Comparison, sep = " | "), median)
cat("median shared edges (%):\n"); print(round(med, 1))

## ---------------------------------------------------------------------------
## S8-2. Regulator out-degree across seeded runs (Figure 4B, Figure S2B)
## ---------------------------------------------------------------------------
odc <- read.csv(o("seed_outdegree_core.csv"), stringsAsFactors = FALSE); odc$network <- "core"
ods <- read.csv(o("seed_outdegree_species.csv"), stringsAsFactors = FALSE)
odl <- rbind(odc[, c("network", "seed", "regulator", "out_degree")],
             ods[, c("network", "seed", "regulator", "out_degree")])
od_rows <- list()
for (nm in names(NETWORKS)) {
  d <- odl[odl$network == nm, c("regulator", "seed", "out_degree")]
  w <- reshape(d, idvar = "regulator", timevar = "seed", direction = "wide")
  names(w) <- sub("^out_degree\\.", "seed", names(w))
  sc <- paste0("seed", SEEDS)
  for (s in sc) if (!s %in% names(w)) w[[s]] <- NA
  v <- as.matrix(w[, sc])
  od_rows[[nm]] <- data.frame(
    Network   = unname(NET_LAB[nm]),
    locus_tag = w$regulator,
    w[, sc],
    Seed_mean = round(rowMeans(v, na.rm = TRUE), 2),
    Seed_sd   = round(apply(v, 1, sd, na.rm = TRUE), 2),
    Seeds_present = rowSums(!is.na(v)),
    stringsAsFactors = FALSE, check.names = FALSE)
}
od <- do.call(rbind, od_rows)
od <- od[c_order(match(od$Network, NET_LAB), -od$Seed_mean, od$locus_tag), ]
## the core rows must agree with Figure 4B
B <- read.csv(o("fig4_B_outdegree.csv"), stringsAsFactors = FALSE)
odcore <- od[od$Network == NET_LAB[["core"]], ]
i <- match(B$locus_tag, odcore$locus_tag)
stopifnot(!any(is.na(i)), max(abs(odcore$Seed_mean[i] - round(B$seed_mean, 2))) < 1e-9)
cat(sprintf("out-degree: %d regulators across %d networks\n", nrow(od), length(NETWORKS)))

## ---------------------------------------------------------------------------
## S8-3. Per-seed centralities, core GRN
## ---------------------------------------------------------------------------
## Computed in 12_figure4_core_stability.R with the functions used for the
## representative networks, so S8-3 and S7-1 are directly comparable.
cent <- read.csv(o("seed_centralities_core.csv"), stringsAsFactors = FALSE,
                 check.names = FALSE)
names(cent) <- sub("^kcore", "k_core", names(cent))
## same stress labels as S8-4
if ("stress_related" %in% names(cent)) {
  cent$stress_related <- ifelse(cent$stress_related == "yes", "Y", "N")
  names(cent)[names(cent) == "stress_related"] <- "Stress_coupled"
  names(cent)[names(cent) == "stress_category"] <- "Stress_category"
}
names(cent)[names(cent) == "seed"] <- "Seed"
front <- c("Seed", "locus_tag", "TF_name", "IC", "IC_rank")
front <- front[front %in% names(cent)]
cent <- cent[, c(front, setdiff(names(cent), front))]
cent <- cent[order(cent$Seed, cent$IC_rank), ]
stopifnot(length(unique(cent$Seed)) == length(SEEDS))
cat(sprintf("per-seed centralities: %d rows, %d to %d regulators per seed\n",
            nrow(cent), min(table(cent$Seed)), max(table(cent$Seed))))

## ---------------------------------------------------------------------------
## S8-4. Core top-15 recurrence
## ---------------------------------------------------------------------------
rec <- read.csv(o("fig4_D_top15_recurrence.csv"), stringsAsFactors = FALSE)
names(rec) <- c("locus_tag", "TF_name", "Representative_IC_rank",
                "Seeds_in_top15", "Seed_rank_min", "Seed_rank_median",
                "Seed_rank_max", "Stress_coupled")
rec$Stress_coupled <- ifelse(rec$Stress_coupled == "yes", "Y", "N")
rec <- rec[order(-rec$Seeds_in_top15, rec$Representative_IC_rank), ]
n_all  <- sum(rec$Seeds_in_top15 == length(SEEDS))
n_most <- sum(rec$Seeds_in_top15 >= 8)
cat(sprintf("top-%d recurrence: %d regulators in all %d runs, %d in at least 8 (%d stress-coupled)\n",
            TOPN, n_all, length(SEEDS), n_most,
            sum(rec$Seeds_in_top15 >= 8 & rec$Stress_coupled == "Y")))

## ---------------------------------------------------------------------------
## S8-5. Stress enrichment per seeded run, core GRN
## ---------------------------------------------------------------------------
cc <- read.csv(o("centralities_core.csv"), stringsAsFactors = FALSE)
N <- nrow(cc); K <- sum(cc$stress_related == "yes")
E <- read.csv(o("fig4_E_stress_by_seed.csv"), stringsAsFactors = FALSE)
enr <- data.frame(
  Run = sprintf("seed%d", E$seed),
  Stress_coupled_in_top15 = E$k_stress,
  Expected_at_random = round(TOPN * K / N, 2),
  Hypergeometric_p = signif(E$p, 3),
  Significant_at_0.05 = ifelse(E$p < 0.05, "yes", "no"),
  stringsAsFactors = FALSE)
k_rep <- sum(cc$stress_related[cc$IC_rank <= TOPN] == "yes")
p_rep <- phyper(k_rep - 1, K, N - K, TOPN, lower.tail = FALSE)
enr <- rbind(data.frame(Run = "representative",
                        Stress_coupled_in_top15 = k_rep,
                        Expected_at_random = round(TOPN * K / N, 2),
                        Hypergeometric_p = signif(p_rep, 3),
                        Significant_at_0.05 = ifelse(p_rep < 0.05, "yes", "no"),
                        stringsAsFactors = FALSE), enr)
## p-values recomputed here must match those behind Figure 4E
stopifnot(max(abs(signif(phyper(E$k_stress - 1, K, N - K, TOPN,
                               lower.tail = FALSE), 3) - enr$Hypergeometric_p[-1])) == 0)
cat(sprintf("stress enrichment: representative %d of %d (p = %s), seeded median %d, %d of %d runs significant\n",
            k_rep, TOPN, signif(p_rep, 3), median(E$k_stress),
            sum(E$p < 0.05), length(SEEDS)))

## ---------------------------------------------------------------------------
## S8-6. Recurrence of each representative network in its seeded runs
## ---------------------------------------------------------------------------
R6 <- read.csv(o("network_recurrence_summary.csv"), stringsAsFactors = FALSE)
nrec <- data.frame(
  Network = R6$network,
  Edges_in_representative_network = R6$representative_edges,
  Edges_recurring_in_at_least_1_run_percent = R6$edges_in_at_least_1_run_percent,
  Edges_recurring_in_at_least_5_runs_percent = R6$edges_in_at_least_5_runs_percent,
  Edges_recurring_in_all_10_runs_percent = R6$edges_in_all_10_runs_percent,
  Gene_overlap_median_percent = R6$node_overlap_median_percent,
  Regulator_overlap_median_percent = R6$regulator_overlap_median_percent,
  stringsAsFactors = FALSE)
stopifnot(nrow(nrec) == length(NETWORKS))

## ---------------------------------------------------------------------------
## Legend
## ---------------------------------------------------------------------------
legend <- rbind(
  c("Dataset S8. Reproducibility of the gene regulatory networks across ten seeded GENIE3 runs", ""),
  c("", ""),
  c(sprintf("Each network was rebuilt %d times with fixed random seeds, using the same input expression compendium and the same pruning as the representative network. Dataset S6 describes the core GRN and Dataset S7 the four representative networks.", length(SEEDS)), ""),
  c("", ""),
  c("Sheet", "Contents"),
  c(SHEETS[["ov"]],   sprintf("Percentage of edges shared between runs, for every pair of the %d seeded runs of each network and between the representative network and each seeded run. Shared edges are expressed as a percentage of the mean size of the two pruned edge sets (Figure 4A, Figure S2A).", length(SEEDS))),
  c(SHEETS[["od"]],   sprintf("Out-degree of every regulator in each of the %d pruned seeded runs of each network, with the mean and standard deviation across the runs in which it has targets (Figure 4B, Figure S2B).", length(SEEDS))),
  c(SHEETS[["cent"]], sprintf("Centrality measures, normalized values, Integrated Centrality and IC rank for the core GRN regulators in each of the %d seeded runs. Computed with the same functions as Dataset S7, so the two are directly comparable.", length(SEEDS))),
  c(SHEETS[["rec"]],  sprintf("How many of the %d seeded runs place each core regulator in the top %d by IC, with the range of its rank across runs.", length(SEEDS), TOPN)),
  c(SHEETS[["enr"]],  sprintf("Number of stress-coupled regulators among the top %d by IC in each seeded run and in the representative network, with the hypergeometric test against a pool of %d regulators of which %d are stress-coupled.", TOPN, N, K)),
  c(SHEETS[["nrec"]], sprintf("How much of each representative network recurs in its %d seeded runs: the percentage of its edges found in at least 1, at least 5 and all %d runs, and the median overlap of its genes and of its regulators with each seeded run.", length(SEEDS), length(SEEDS))),
  c("", ""),
  c("Column", "Definition"),
  c("Network", "Core GRN or the species network the row refers to."),
  c("Run_A, Run_B, Run", "Seeded run, identified by its seed, or the unseeded representative network."),
  c("Comparison", "Whether the row compares two seeded runs or the representative network against a seeded run (S8-1)."),
  c("Shared_edges_percent", "Edges present in both runs, as a percentage of the mean size of the two edge sets."),
  c("seed1 to seed10", "Value in the run with that seed."),
  c("Seed_mean, Seed_sd, Seeds_present", "Mean and standard deviation across seeded runs, and the number of runs in which the regulator was present."),
  c("degree, k_core, betweenness, stress, eigenvector", "Centrality measures, computed as described in Dataset S6."),
  c("*_norm, IC, IC_rank", "Measure divided by its maximum across the ranked regulators of that run; their sum (Eq. 2); and the rank by IC within that run."),
  c("Representative_IC_rank", "Rank of the regulator by IC in the representative network."),
  c("Seeds_in_top15", sprintf("Number of seeded runs placing the regulator in the top %d by IC.", TOPN)),
  c("Seed_rank_min, Seed_rank_median, Seed_rank_max", "Range and median of the regulator's IC rank across the seeded runs in which it was ranked."),
  c("Stress_coupled", "Y if a stress role is assigned to the regulator in Dataset S6, N otherwise."),
  c("Stress_category", "Stress role assigned to the regulator (Dataset S6-2)."),
  c("Expected_at_random", sprintf("Stress-coupled regulators expected among %d drawn at random from the pool, that is %d x %d / %d.", TOPN, TOPN, K, N)),
  c("Hypergeometric_p", "One-sided probability of observing at least this many stress-coupled regulators by chance."),
  c("Edges_recurring_in_...", "Edges of the representative network (regulator-target pairs) that are also present in the given number of pruned seeded runs, as a percentage of its edges."),
  c("Gene_overlap_median_percent, Regulator_overlap_median_percent", "Genes, or regulators (nodes with targets), present in both the representative network and a seeded run, as a percentage of the mean count of the two; median across the seeded runs.")
)
legend <- as.data.frame(legend, stringsAsFactors = FALSE); names(legend) <- NULL

## ---------------------------------------------------------------------------
## Workbook
## ---------------------------------------------------------------------------
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold")
add_table <- function(name, d) {
  addWorksheet(wb, name)
  writeData(wb, name, d, colNames = TRUE, rowNames = FALSE, headerStyle = hs)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = 2)
  setColWidths(wb, name, cols = seq_len(ncol(d)), widths = "auto")
  cat(sprintf("  %-34s %5d rows, %2d columns\n", name, nrow(d), ncol(d)))
}
addWorksheet(wb, SHEETS[["legend"]])
writeData(wb, SHEETS[["legend"]], legend, colNames = FALSE, rowNames = FALSE)
addStyle(wb, SHEETS[["legend"]], hs, rows = c(1, 5, 13), cols = 1:2, gridExpand = TRUE)
setColWidths(wb, SHEETS[["legend"]], cols = 1:2, widths = c(40, 120))
add_table(SHEETS[["ov"]],   ov)
add_table(SHEETS[["od"]],   od)
add_table(SHEETS[["cent"]], cent)
add_table(SHEETS[["rec"]],  rec)
add_table(SHEETS[["enr"]],  enr)
add_table(SHEETS[["nrec"]], nrec)
saveWorkbook(wb, OUT, overwrite = TRUE)

stopifnot(nrow(read.xlsx(OUT, sheet = SHEETS[["ov"]]))   == nrow(ov),
          nrow(read.xlsx(OUT, sheet = SHEETS[["cent"]])) == nrow(cent),
          nrow(read.xlsx(OUT, sheet = SHEETS[["enr"]]))  == length(SEEDS) + 1)
cat("written ", OUT, "\n", sep = "")
