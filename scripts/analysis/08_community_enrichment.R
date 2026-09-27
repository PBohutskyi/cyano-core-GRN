# 08_community_enrichment.R
# COG and KEGG enrichment of the Louvain communities of the core GRN (Figure 3B,
# Results 3.4.1). Same one-sided hypergeometric test as Table 1 (Eq. 1), applied
# to the genes of each community against the PCC 7942 genome.
#
# Universe: all 2,712 protein-coding genes of PCC 7942. Genes without a COG
# category (or without a KEGG pathway) form an "unassigned" category, which is
# tested with the others. Benjamini-Hochberg q-values are computed within each
# community and source, over the categories with at least one community gene.
#
# Inputs : networks/core/network_of_record.graphml   Louvain community of each node
#          reference/ref_protein_coding_7942.csv       COG and KEGG annotation
#          curation/core_regulator_curation.csv        regulator names
# Outputs: outputs/tables/community_enrichment.csv     all tested categories
#          outputs/tables/community_enrichment_significant.csv   q < 0.05

source("scripts/common.R")
Q_CUTOFF <- 0.05

ann <- read.csv("reference/ref_protein_coding_7942.csv", stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
UNIVERSE <- unique(ann$locus_tag)

## term -> gene tables; genes without a category go to "-" (unassigned)
split_terms <- function(field, sep, keep) {
  out <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
    v <- ann[[field]][i]
    t <- if (is.na(v) || trimws(v) %in% c("", "-")) character(0) else keep(v, sep)
    if (!length(t)) t <- "-"
    data.frame(term = t, gene = ann$locus_tag[i], stringsAsFactors = FALSE)
  }))
  unique(out)
}
cog  <- split_terms("COG_category", "", function(v, s) setdiff(strsplit(v, s)[[1]], c("-", " ")))
kegg <- split_terms("KEGG_Pathway", ",", function(v, s) {
  p <- trimws(strsplit(v, s)[[1]]); p[grepl("^map[0-9]+$", p)] })

enrich <- function(genes, t2g, source_label) {
  genes <- intersect(genes, UNIVERSE)
  N <- length(UNIVERSE); n <- length(genes)
  res <- do.call(rbind, lapply(unique(t2g$term), function(tm) {
    inTerm <- t2g$gene[t2g$term == tm]
    k <- length(intersect(genes, inTerm)); K <- length(inTerm)
    if (k == 0) return(NULL)
    data.frame(source = source_label, term = tm, gene_count = k, set_size = n,
               term_size = K, universe = N, fold = (k / n) / (K / N),
               p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  }))
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res
}

## communities from the network of record
g  <- read_graph(rep_file(NETWORKS$core), format = "graphml")
nm <- node_names(g)
cl <- as.integer(vertex_attr(g, "louvain_cluster"))
cur <- core_curation()
reg_of <- setNames(cur$TF_name, cur$locus_tag)

rows <- list()
for (c0 in sort(unique(cl[!is.na(cl) & cl > 0]))) {
  genes <- nm[cl == c0]
  regs  <- intersect(genes, names(reg_of))
  lab   <- if (length(regs)) paste(reg_of[regs], collapse = ", ") else ""
  r <- rbind(enrich(genes, cog, "COG"), enrich(genes, kegg, "KEGG"))
  r <- r[r$term != "-", ]
  if (!nrow(r)) next
  rows[[length(rows) + 1]] <- cbind(community = c0, regulators = lab,
                                    community_size = length(genes), r)
}
all <- do.call(rbind, rows)
all <- all[order(all$community, all$source, all$p_value), ]
all$gene_ratio <- paste0(all$gene_count, "/", all$set_size)
write.csv(all, o("community_enrichment.csv"), row.names = FALSE)
sig <- all[all$q_value < Q_CUTOFF, ]
write.csv(sig, o("community_enrichment_significant.csv"), row.names = FALSE)

cat("communities tested:", length(rows), "| significant categories (q < 0.05):", nrow(sig), "\n")
print(sig[, c("community", "regulators", "source", "term", "gene_ratio", "fold", "q_value")],
      row.names = FALSE, digits = 3)
