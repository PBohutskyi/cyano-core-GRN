# build_dataset_S6.R
# Dataset S6: core GRN nodes, regulators, edges, and the robustness,
# enrichment and cross-validation tables. Node measures are recomputed from the
# core network with the functions of common.R; every other number is read from
# a pipeline output.
#
# Inputs : networks/core/network_of_record.graphml
#          curation/core_regulator_curation.csv   names and classification
#          curation/S6_source_annotation.xlsx     gene and edge annotation
#          outputs/tables/centralities_core.csv, S6-4 to S6-7 tables,
#          xval_matches.csv                     (scripts 05, 11, 15)
# Output : datasets/Dataset_S6_Core_GRN_nodes_TFs_and_edges.xlsx

source("scripts/common.R")
suppressPackageStartupMessages(library(openxlsx))

OUT <- "datasets/Dataset_S6_Core_GRN_nodes_TFs_and_edges.xlsx"
## Figure numbers cited in the legend sheet. Kept here so a renumbering of the
## manuscript figures is a one-line change rather than a search through prose.
FIG_XVAL     <- 5    # cross-validation of core and species high-influence regulators
FIG_CLUSTERS <- 3    # core GRN with Louvain communities

SHEETS <- c(legend = "S6-0. Legend",
            nodes  = "S6-1. All nodes in network",
            tfs    = "S6-2. Transcription regulators",
            edges  = "S6-3. List of TF to gene edges",
            s4     = "S6-4. IC robustness",
            s5     = "S6-5. IC enrichment tests",
            s6     = "S6-6. Xval enrichment",
            s7     = "S6-7. Xval sensitivity",
            s8     = "S6-8. Xval matches")

need <- function(f) { if (!file.exists(f)) stop("missing ", f); f }
SRC <- need(CUR_S6)
IN  <- TAB
for (f in c("centralities_core.csv", "S6-4_IC_robustness.csv",
            "S6-5_IC_enrichment.csv", "S6-6_xval_enrichment.csv",
            "S6-7_xval_sensitivity.csv", "xval_matches.csv"))
  need(file.path(IN, f))

rd <- function(sheet) as.data.frame(readxl::read_excel(SRC, sheet = sheet,
                                                       .name_repair = "minimal"))

## ---------------------------------------------------------------------------
## 1. Core graph and node measures (same functions as 05_centralities.R)
## ---------------------------------------------------------------------------
g  <- read_network_of_record(NETWORKS$core)
nm <- node_names(g)
stopifnot(vcount(g) == 889, ecount(g) == 1100, !any(duplicated(nm)))

meas <- node_measures(g)                 # all 889 nodes
meas$indegree <- degree(g, mode = "in")
names(meas)[names(meas) == "kcore"] <- "k_core"
meas <- meas[, c("locus_tag", "degree", "indegree", "outdegree",
                 "k_core", "betweenness", "stress", "eigenvector")]
cat("graph: ", vcount(g), " nodes, ", ecount(g), " edges, ",
    sum(meas$outdegree > 0), " nodes with outgoing edges\n", sep = "")

## ---------------------------------------------------------------------------
## 2. S6-1 all nodes: annotation + recomputed measures
## ---------------------------------------------------------------------------
ann_nodes <- rd(SHEETS[["nodes"]])
stopifnot(nrow(ann_nodes) == 889, setequal(ann_nodes$locus_tag, nm))
ANNOT <- c("locus_tag", "start", "end", "strand", "feature_interval_length",
           "product_length", "product_accession", "name", "symbol", "GeneID",
           "old_locus_tag")
nodes <- ann_nodes[, c(ANNOT, "Louvain_cluster", "TF")]
nodes$TF <- ifelse(nodes$TF %in% 1, "Y", "N")
names(nodes)[names(nodes) == "TF"] <- "Regulator"
nodes <- merge(nodes, meas, by = "locus_tag", sort = FALSE)
nodes <- nodes[order(nodes$locus_tag), ]
stopifnot(nrow(nodes) == 889, sum(nodes$Regulator == "Y") == 38,
          all(nodes$outdegree[nodes$Regulator == "N"] == 0))
## degree from the graph must agree with the degree in the annotation
stopifnot(all(nodes$degree == ann_nodes$degree[match(nodes$locus_tag, ann_nodes$locus_tag)]))

## ---------------------------------------------------------------------------
## 3. S6-2 regulators: curated classification + centralities_core.csv
## ---------------------------------------------------------------------------
ann_tfs <- rd(SHEETS[["tfs"]])
stopifnot(nrow(ann_tfs) == 38, all(ann_tfs$locus_tag %in% nodes$locus_tag[nodes$Regulator == "Y"]))

tfs <- ann_tfs[, c(ANNOT, "Original description", "class", "type",
               "P2TF description", "Louvain_cluster")]
cur <- core_curation()
stopifnot(setequal(cur$locus_tag, tfs$locus_tag))
tfs <- merge(tfs, cur[, c("locus_tag", "TF_name", "Stress_category",
                          "Non_stress_category", "Assignment_confidence",
                          "Assignment_rationale")],
             by = "locus_tag", sort = FALSE)
sc <- trimws(ifelse(is.na(tfs$Stress_category), "", tfs$Stress_category))
tfs$Stress_category <- ifelse(sc == "", NA, sc)
tfs$Stress_coupled  <- ifelse(sc == "", "N", "Y")

stopifnot(all(tfs$Assignment_confidence %in% c("low", "med", "high")))

cc <- read.csv(file.path(IN, "centralities_core.csv"), stringsAsFactors = FALSE,
               check.names = FALSE)
stopifnot(nrow(cc) == nrow(tfs), setequal(cc$locus_tag, tfs$locus_tag))
## classification in the CSV (via annotate_regulators) must match the curation
m <- match(cc$locus_tag, tfs$locus_tag)
stopifnot(all(cc$stress_related == ifelse(tfs$Stress_coupled[m] == "Y", "yes", "no")))
## IC must equal the sum of the five normalized measures, Eq. 2
stopifnot(max(abs(rowSums(cc[, paste0(MEASURES, "_norm")]) - cc$IC)) < 1e-9)
stopifnot(all(sort(cc$IC_rank) == seq_len(nrow(cc))))

cc <- cc[, setdiff(names(cc), c("TF_name", "stress_related", "stress_category"))]
names(cc) <- sub("^kcore", "k_core", names(cc))
## measures come from the graph for every regulator; the CSV supplies the
## normalized columns, the seeded statistics and the rank
cc <- cc[, setdiff(names(cc), c("degree", "k_core", "betweenness",
                                "eigenvector", "stress"))]
tfs <- merge(tfs, meas, by = "locus_tag", all.x = TRUE, sort = FALSE)
tfs <- merge(tfs, cc,   by = "locus_tag", all.x = TRUE, sort = FALSE)
stopifnot(!any(is.na(tfs$degree)), !any(is.na(tfs$IC)),
          all(sort(tfs$IC_rank) == seq_len(nrow(tfs))))
front <- c(ANNOT, "Original description", "class", "type", "P2TF description",
           "TF_name", "Louvain_cluster", "Stress_category", "Non_stress_category",
           "Stress_coupled", "Assignment_confidence", "Assignment_rationale",
           "indegree", "outdegree")
tfs <- tfs[, c(front, setdiff(names(tfs), front))]
tfs <- tfs[order(tfs$IC_rank), ]

## headline numbers, recomputed here and printed for the log
pool <- tfs
K <- sum(pool$Stress_coupled == "Y"); N <- nrow(pool)
top15 <- pool[order(pool$IC_rank), ][1:15, ]
k15 <- sum(top15$Stress_coupled == "Y")
p15 <- phyper(k15 - 1, K, N - K, 15, lower.tail = FALSE)
cat(sprintf("regulators: %d ranked, %d with outgoing edges, %d stress-coupled (%.0f%%)\n",
            N, sum(tfs$outdegree > 0), K, 100 * K / N))
cat(sprintf("top 15 by IC: %d stress-coupled, hypergeometric p = %.4f\n", k15, p15))

## ---------------------------------------------------------------------------
## 4. S6-3 edges: annotation order and target annotation, checked against the graph
## ---------------------------------------------------------------------------
ann_edges <- rd(SHEETS[["edges"]])
stopifnot(nrow(ann_edges) == 1100)
el <- as_edgelist(g, names = FALSE)
gedges <- paste(nm[el[, 1]], nm[el[, 2]])
vedges <- paste(ann_edges[["Regulator (TF) locus_tag"]], ann_edges[["Gene Target"]])
stopifnot(setequal(gedges, vedges), !any(duplicated(vedges)))

edges <- data.frame(
  Edge                 = ann_edges$Edge,
  Regulator_locus_tag  = ann_edges[["Regulator (TF) locus_tag"]],
  Regulator_name       = tfs$TF_name[match(ann_edges[["Regulator (TF) locus_tag"]], tfs$locus_tag)],
  Regulator_Louvain_cluster = nodes$Louvain_cluster[match(ann_edges[["Regulator (TF) locus_tag"]], nodes$locus_tag)],
  Target_locus_tag     = ann_edges[["Gene Target"]],
  Target_gene_name     = ann_edges[["Targeted gene name"]],
  Target_description_eggNOG = ann_edges[["Gene description Eggnog"]],
  stringsAsFactors = FALSE)
stopifnot(!any(is.na(edges$Regulator_name)), !any(is.na(edges$Regulator_Louvain_cluster)),
          length(unique(edges$Regulator_locus_tag)) == 36)

## ---------------------------------------------------------------------------
## 5. S6-4 to S6-7: tables from 11_IC_enrichment_robustness.R and 15_crossvalidation.R
## ---------------------------------------------------------------------------
read_block <- function(f) {
  d <- read.csv(file.path(IN, f), stringsAsFactors = FALSE, check.names = FALSE,
                colClasses = "character", na.strings = character(0))
  d[is.na(d)] <- ""
  keep <- vapply(d, function(x) any(x != ""), logical(1))   # drop empty columns
  d <- d[, keep, drop = FALSE]
  names(d) <- NULL
  as.matrix(d)
}
s4 <- read_block("S6-4_IC_robustness.csv")
s5 <- read_block("S6-5_IC_enrichment.csv")
s6 <- read_block("S6-6_xval_enrichment.csv")
s7 <- read_block("S6-7_xval_sensitivity.csv")

## cross-check the headline row of S6-5 against the recomputation above
r15 <- s5[s5[, 1] == "15", , drop = FALSE]
stopifnot(nrow(r15) == 1, as.integer(r15[1, 2]) == k15,
          abs(as.numeric(r15[1, 5]) - p15) < 5e-5)
## S6-4 to S6-7 must have been regenerated against the same regulator pool.
## Stale files from an earlier pool would otherwise ship silently.
for (nmb in c("S6-4", "S6-5", "S6-6", "S6-7")) {
  b <- get(tolower(sub("S6-", "s", nmb)))
  hit <- grep(sprintf("(pool|Background|Pool):? .*\\b%d\\b|pool of %d", N, N), b, value = TRUE)
  if (!length(grep(sprintf("\\b%d\\b", N), b[grep("pool|Pool|Background", b)])))
    stop(nmb, " does not match a pool of ", N, " regulators")
}

## cross-check S6-6 header against xval_matches.csv
xv <- read.csv(file.path(IN, "xval_matches.csv"), stringsAsFactors = FALSE)
hdr6 <- s6[3, 1]
stopifnot(grepl(sprintf("^%d matches", nrow(xv)), hdr6),
          grepl(sprintf("%d distinct", length(unique(xv$regulator))), hdr6))

## S6-8: the individual cross-validation matches behind S6-6. Each row is one
## regulator recovered in one species network by one measure. S6-6 gives the
## counts; this sheet shows which regulator matched where.
SPECIES_LAB <- c("7942" = "S. elongatus PCC 7942",
                 "6803" = "Synechocystis sp. PCC 6803",
                 "7002" = "Picosynechococcus sp. PCC 7002")
xi <- xv[, c("locus_tag", "regulator", "species", "measure", "group")]
names(xi) <- c("locus_tag", "TF_name", "Species", "Measure", "Measure_group")
xi$Species  <- unname(SPECIES_LAB[as.character(xi$Species)])
xi$Measure  <- ifelse(xi$Measure == "IC", "IC", sub("^kcore$", "k_core", xi$Measure))
xi$Stress_coupled <- tfs$Stress_coupled[match(xi$locus_tag, tfs$locus_tag)]
xi$Core_IC_rank   <- tfs$IC_rank[match(xi$locus_tag, tfs$locus_tag)]
stopifnot(!any(is.na(xi$Stress_coupled)), !any(is.na(xi$Species)))
MEAS_ORDER <- c("degree", "k_core", "betweenness", "stress", "eigenvector", "IC")
xi <- xi[c_order(xi$Core_IC_rank, xi$Species, match(xi$Measure, MEAS_ORDER)), ]
cat(sprintf("cross-validation: %d matches, %d distinct regulators, %d stress-coupled\n",
            nrow(xi), length(unique(xi$locus_tag)),
            length(unique(xi$locus_tag[xi$Stress_coupled == "Y"]))))

## ---------------------------------------------------------------------------
## 6. Legend sheet
## ---------------------------------------------------------------------------
legend <- rbind(
  c("Dataset S6. Conserved core gene regulatory network (GRN): nodes, transcriptional regulators, edges, and robustness tests", ""),
  c("", ""),
  c("Sheet", "Contents"),
  c(SHEETS[["nodes"]], sprintf("All %d nodes of the pruned core GRN with RefSeq annotation for S. elongatus PCC 7942, Louvain cluster, and centrality measures computed with igraph on the directed network (eigenvector on the undirected network).", nrow(nodes))),
  c(SHEETS[["tfs"]], sprintf("The %d transcriptional regulators in the pruned core GRN with curated name, functional category, stress classification, centrality measures, Integrated Centrality (IC, Eq. 2) and IC rank. All %d are ranked; %d of them carry outgoing edges after pruning, and the remaining %d retained only an incoming edge and therefore rank last.", nrow(tfs), N, sum(tfs$outdegree > 0), N - sum(tfs$outdegree > 0))),
  c(SHEETS[["edges"]], sprintf("The %d regulator-to-gene edges of the pruned core GRN with target gene name and eggNOG description.", nrow(edges))),
  c(SHEETS[["s4"]],  "Robustness of the IC ranking: leave-one-measure-out, correlation of each measure with IC, dominance after max-normalization, and alternative normalizations."),
  c(SHEETS[["s5"]],  "Enrichment of stress-coupled regulators among the top-ranked regulators by IC, from the top 5 to the top 20, and sensitivity to the classification of individual regulators."),
  c(SHEETS[["s6"]],  sprintf("Regulators cross-validated between the core GRN and the species GRNs (Figure %d), with enrichment of stress-coupled regulators in that set.", FIG_XVAL)),
  c(SHEETS[["s7"]],  "Sensitivity of the cross-validation to the core and species rank cutoffs."),
  c(SHEETS[["s8"]],  sprintf("The individual matches behind the cross-validation: %d matches of a core regulator recovered in a species network by one measure, covering %d distinct regulators. S6-6 gives the counts; this sheet gives the matches.", nrow(xi), length(unique(xi$locus_tag)))),
  c("", ""),
  c("Column", "Definition"),
  c("locus_tag, old_locus_tag", "Current RefSeq and previous locus tags in PCC 7942."),
  c("Louvain_cluster", sprintf("Cluster identifier from Louvain clustering (Methods); cluster 0 is the group shown as n/a in Figure %dB.", FIG_CLUSTERS)),
  c("Regulator", "Y if the node is one of the 38 curated regulators, N otherwise (S6-1)."),
  c("TF_name", "Standardized regulator name used in figures and tables."),
  c("Stress_category / Non_stress_category", "Stress role and primary non-stress role of the regulator. A regulator is stress-coupled (Stress_coupled = Y) when a stress role is assigned."),
  c("Assignment_confidence, Assignment_rationale", "Confidence (low, med, high) and one-line basis for the functional assignment."),
  c("degree, indegree, outdegree", "Number of edges touching the node; incoming; outgoing."),
  c("k_core", "Largest k such that the node belongs to the k-core of the undirected network."),
  c("betweenness", "Number of shortest directed paths through the node, each path weighted by 1 / number of equivalent shortest paths."),
  c("stress", "Number of shortest directed paths through the node."),
  c("eigenvector", "Eigenvector centrality on the undirected network, scaled to a maximum of 1."),
  c("*_norm", "Measure divided by its maximum across the ranked regulators (S6-2)."),
  c("IC", "Integrated Centrality, Eq. 2: sum of degree_norm, k_core_norm, betweenness_norm, eigenvector_norm and stress_norm."),
  c("*_seedmean, *_seedsd", "Mean and standard deviation of the value across ten GENIE3 runs with fixed seeds."),
  c("IC_rank", "Rank by IC among the ranked regulators; 1 = highest."),
  c("IC_z_vs_seeds", "(IC - IC_seedmean) / IC_seedsd for the representative network."),
  c("Regulator_Louvain_cluster", "Louvain cluster of the regulator (S6-3)."),
  c("Species, Measure, Measure_group", "Species network in which the regulator was recovered, the measure that recovered it, and the scope of that measure (S6-8)."),
  c("Core_IC_rank", "Rank of the regulator by IC in the core GRN (S6-8).")
)
legend <- as.data.frame(legend, stringsAsFactors = FALSE); names(legend) <- NULL

## ---------------------------------------------------------------------------
## 7. Workbook
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
add_block <- function(name, m) {
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(m, stringsAsFactors = FALSE),
            colNames = FALSE, rowNames = FALSE)
  addStyle(wb, name, hs, rows = 1, cols = 1)
  setColWidths(wb, name, cols = seq_len(ncol(m)), widths = "auto")
  cat(sprintf("  %-34s %5d rows, %2d columns\n", name, nrow(m), ncol(m)))
}
addWorksheet(wb, SHEETS[["legend"]])
writeData(wb, SHEETS[["legend"]], legend, colNames = FALSE, rowNames = FALSE)
addStyle(wb, SHEETS[["legend"]], hs, rows = c(1, 3, 12), cols = 1:2, gridExpand = TRUE)
setColWidths(wb, SHEETS[["legend"]], cols = 1:2, widths = c(40, 120))
add_table(SHEETS[["nodes"]], nodes)
add_table(SHEETS[["tfs"]],   tfs)
add_table(SHEETS[["edges"]], edges)
add_block(SHEETS[["s4"]], s4)
add_block(SHEETS[["s5"]], s5)
add_block(SHEETS[["s6"]], s6)
add_block(SHEETS[["s7"]], s7)
add_table(SHEETS[["s8"]], xi)
saveWorkbook(wb, OUT, overwrite = TRUE)

## verify what landed on disk
chk <- read.xlsx(OUT, sheet = SHEETS[["tfs"]])
stopifnot(nrow(chk) == 38, "IC_rank" %in% names(chk),
          nrow(read.xlsx(OUT, sheet = SHEETS[["nodes"]])) == 889,
          nrow(read.xlsx(OUT, sheet = SHEETS[["edges"]])) == 1100,
          nrow(read.xlsx(OUT, sheet = SHEETS[["s8"]])) == nrow(xi))
cat("written ", OUT, "\n", sep = "")
