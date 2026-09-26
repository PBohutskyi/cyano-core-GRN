# 06_table3_network_properties.R
# Topological properties of the four networks of record: Table 3 for the
# species GRNs and the core GRN values given in Results.
#
# Inputs : networks/<network>/network_of_record.graphml
#          networks/<network>/regulators_GENIE3_input.csv
#          datasets/Dataset_S3_Core_genome_homologs.xlsx   core gene set size
#          reference/genome_gene_lists.csv                 genome sizes
#          curation/core_regulator_curation.csv
# Output : outputs/tables/table3_network_properties.csv
#
# Definitions (igraph, unweighted):
#   regulators_input        regulators supplied to GENIE3 for this network
#   regulators_in_network   of those, present as nodes after pruning
#   regulators_with_targets of those, with at least one outgoing edge
#   regulators_ranked       the ranked pool used for IC (ranked_pool() in common.R)
#   genes_pct_of_gene_set   nodes / genes in the reference set x 100; the
#                           reference set is the tri-homologous core genome for
#                           the core GRN and the protein-coding genome otherwise
#   edges_per_node          edges / nodes
#   density                 edge_density() of the directed network
#   clustering              mean local clustering of the undirected, collapsed
#                           network; nodes with fewer than two neighbours count 0
#   transitivity_global     global transitivity of the same network
#   median_regulator_degree median total degree of the ranked regulators
#   median_regulator_betweenness
#                           median directed betweenness of the ranked
#                           regulators, divided by the largest of them
#   top10_edge_share_pct    edges leaving the ten regulators with the highest
#                           out-degree, as % of all edges

source("scripts/common.R")

core_genes <- nrow(as.data.frame(read_excel("datasets/Dataset_S3_Core_genome_homologs.xlsx",
                                            sheet = "S3-2. Tri-homologous core")))
gl <- read.csv("reference/genome_gene_lists.csv", stringsAsFactors = FALSE)
n_genes <- function(s) sum(!is.na(gl[[s]]) & trimws(gl[[s]]) != "")
GENE_SET <- c(core = core_genes, s7942 = n_genes("PCC7942"),
              s6803 = n_genes("PCC6803"), s7002 = n_genes("PCC7002"))
print(GENE_SET)

read_tf_list <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  v <- unique(trimws(as.character(d$locus_tag)))
  v[!is.na(v) & v != ""]
}

rows <- list()
for (nm in names(NETWORKS)) {
  net <- NETWORKS[[nm]]
  g <- read_network_of_record(net)
  nodes <- node_names(g)
  say(nm, ": ", vcount(g), " nodes, ", ecount(g), " edges")

  tf_in    <- read_tf_list(tf_input_file(net))
  in_net   <- intersect(tf_in, nodes)
  with_tgt <- intersect(in_net, regulators(g))
  pool     <- ranked_pool(g, net)
  stopifnot(length(setdiff(regulators(g), tf_in)) == 0)

  deg  <- setNames(degree(g, mode = "all"), nodes)
  outd <- setNames(degree(g, mode = "out"), nodes)
  btw  <- setNames(betweenness(g, directed = TRUE, weights = NA), nodes)
  u    <- as_undirected(g, mode = "collapse")

  rows[[nm]] <- data.frame(
    network                      = nm,
    nodes                        = vcount(g),
    edges                        = ecount(g),
    regulators_input             = length(tf_in),
    regulators_in_network        = length(in_net),
    regulators_with_targets      = length(with_tgt),
    regulators_ranked            = length(pool),
    genes_pct_of_gene_set        = round(100 * vcount(g) / GENE_SET[[nm]], 1),
    edges_per_node               = round(ecount(g) / vcount(g), 2),
    density                      = signif(edge_density(g), 3),
    clustering                   = round(transitivity(u, type = "localaverage", isolates = "zero"), 3),
    transitivity_global          = round(transitivity(u, type = "global"), 4),
    median_regulator_degree      = median(deg[pool]),
    median_regulator_betweenness = round(median(btw[pool]) / max(btw[pool]), 4),
    top10_edge_share_pct         = round(100 * sum(sort(outd[in_net], decreasing = TRUE)[1:10]) /
                                           ecount(g), 1),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, o("table3_network_properties.csv"), row.names = FALSE)
print(res, row.names = FALSE)
