# build_dataset_S7.R
# Dataset S7: centrality measures, Integrated Centrality and seeded-run
# statistics for the ranked regulators of the four representative networks.
#
# Inputs : outputs/tables/centralities_{core,s7942,s6803,s7002}.csv  (05_centralities.R)
# Output : datasets/Dataset_S7_Centrality_and_stability.xlsx

source("scripts/common.R")
suppressPackageStartupMessages(library(openxlsx))

OUT  <- "datasets/Dataset_S7_Centrality_and_stability.xlsx"
NETS <- c(core = "S7-1. Core GRN", s7942 = "S7-2. PCC 7942",
          s6803 = "S7-3. PCC 6803", s7002 = "S7-4. PCC 7002")

tabs <- lapply(names(NETS), function(nm) {
  f <- o(paste0("centralities_", nm, ".csv"))
  if (!file.exists(f)) stop("missing ", f)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
})
names(tabs) <- names(NETS)

legend <- rbind(
  c("Dataset S7. Centrality measures and Integrated Centrality of the regulators in the core GRN and the three species GRNs", ""),
  c("", ""),
  c("Sheet", "Contents"),
  c(NETS[["core"]],  sprintf("The %d ranked regulators of the core GRN.", nrow(tabs$core))),
  c(NETS[["s7942"]], sprintf("The %d ranked regulators of the S. elongatus PCC 7942 GRN.", nrow(tabs$s7942))),
  c(NETS[["s6803"]], sprintf("The %d ranked regulators of the Synechocystis sp. PCC 6803 GRN.", nrow(tabs$s6803))),
  c(NETS[["s7002"]], sprintf("The %d ranked regulators of the Picosynechococcus sp. PCC 7002 GRN.", nrow(tabs$s7002))),
  c("", ""),
  c("Column", "Definition"),
  c("locus_tag", "RefSeq locus tag."),
  c("TF_name", "Regulator name, where one is assigned (Dataset S5, Dataset S6-2)."),
  c("stress_related, stress_category", "Core GRN only: yes if a stress role is assigned, and that role (Dataset S6-2)."),
  c("degree", "Number of edges touching the regulator."),
  c("kcore", "Largest k such that the regulator belongs to the k-core of the undirected network."),
  c("betweenness", "Number of shortest directed paths through the regulator, each weighted by 1 / number of equivalent shortest paths."),
  c("eigenvector", "Eigenvector centrality on the undirected network, scaled to a maximum of 1."),
  c("stress", "Number of shortest directed paths through the regulator."),
  c("*_norm", "Measure divided by its maximum across the ranked regulators of that network."),
  c("IC", "Integrated Centrality, Eq. 2: sum of the five *_norm values."),
  c("*_seedmean, *_seedsd", "Mean and standard deviation of the value across ten GENIE3 runs with fixed seeds, pruned as the representative network; a regulator absent from a run counts as 0."),
  c("IC_rank", "Rank by IC; 1 = highest."),
  c("IC_z_vs_seeds", "(IC - IC_seedmean) / IC_seedsd.")
)
legend <- as.data.frame(legend, stringsAsFactors = FALSE); names(legend) <- NULL

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold")
addWorksheet(wb, "S7-0. Legend")
writeData(wb, "S7-0. Legend", legend, colNames = FALSE, rowNames = FALSE)
addStyle(wb, "S7-0. Legend", hs, rows = c(1, 3, 9), cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "S7-0. Legend", cols = 1:2, widths = c(34, 120))
for (nm in names(NETS)) {
  addWorksheet(wb, NETS[[nm]])
  writeData(wb, NETS[[nm]], tabs[[nm]], colNames = TRUE, rowNames = FALSE, headerStyle = hs)
  freezePane(wb, NETS[[nm]], firstActiveRow = 2, firstActiveCol = 2)
  cat(sprintf("  %-14s %3d regulators, %2d columns\n", NETS[[nm]], nrow(tabs[[nm]]), ncol(tabs[[nm]])))
}
saveWorkbook(wb, OUT, overwrite = TRUE)
stopifnot(nrow(read.xlsx(OUT, sheet = NETS[["core"]])) == 38)
cat("written ", OUT, "\n", sep = "")
