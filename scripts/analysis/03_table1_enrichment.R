# 03_table1_enrichment.R
# Table 1: COG and KEGG enrichment of the tri-homologous core genome and of the
# two duo-homolog gene sets that include PCC 7942, against the PCC 7942 genome.
# One-sided hypergeometric test; Benjamini-Hochberg q-values alongside.
#
# Inputs : reference/ref_protein_coding_7942.csv            annotation and background
#          datasets/Dataset_S3_Core_genome_homologs.xlsx    gene sets
# Outputs: outputs/tables/table1_enrichment.csv             all terms
#          outputs/tables/table1_significant.csv            p < 0.05, as in Table 1

source("scripts/common.R")
P_CUTOFF <- 0.05

ann <- read.csv("reference/ref_protein_coding_7942.csv", stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
S3  <- "datasets/Dataset_S3_Core_genome_homologs.xlsx"
rd  <- function(sh) as.data.frame(read_excel(S3, sheet = sh, .name_repair = "minimal"))
tags <- function(sh) {
  d <- rd(sh); v <- d[["PCC7942_locus_tag"]]
  unique(v[!is.na(v) & trimws(v) != ""])
}
SETS <- list(
  "PCC 7942 - PCC 6803 - PCC 7002" = tags("S3-2. Tri-homologous core"),
  "PCC 7942 - PCC 7002"            = tags("S3-3. Duo-homologs 7942-7002"),
  "PCC 7942 - PCC 6803"            = tags("S3-4. Duo-homologs 7942-6803"))
cat("gene sets:", paste(names(SETS), sapply(SETS, length), collapse = " | "), "\n")

## ---- term -> gene tables ---------------------------------------------------
## COG: one letter per gene may be several, e.g. "EG". Each letter is a term.
cog <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
  v <- ann$COG_category[i]
  if (is.na(v) || v %in% c("", "-")) return(NULL)
  data.frame(term = strsplit(v, "")[[1]], gene = ann$locus_tag[i], stringsAsFactors = FALSE)
}))
cog <- unique(cog[!cog$term %in% c("-", " "), ])

## KEGG: comma-separated pathway ids. Only ko identifiers are kept, since map
## and ko ids are duplicates of the same pathway.
kegg <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
  v <- ann$KEGG_Pathway[i]
  if (is.na(v) || v == "") return(NULL)
  p <- trimws(strsplit(v, ",")[[1]])
  p <- p[grepl("^ko[0-9]+$", p)]
  if (!length(p)) return(NULL)
  data.frame(term = p, gene = ann$locus_tag[i], stringsAsFactors = FALSE)
}))
kegg <- unique(kegg)
cat("annotation: genome", nrow(ann), "genes |", length(unique(cog$term)), "COG categories |",
    length(unique(kegg$term)), "KEGG pathways\n")


## ---- term descriptions ------------------------------------------------------
COG_DESC <- c(
  J = "Translation, ribosomal structure and biogenesis",
  A = "RNA processing and modification",
  K = "Transcription",
  L = "Replication, recombination and repair",
  B = "Chromatin structure and dynamics",
  D = "Cell cycle control, cell division, chromosome partitioning",
  Y = "Nuclear structure",
  V = "Defense mechanisms",
  T = "Signal transduction mechanisms",
  M = "Cell wall, membrane and envelope biogenesis",
  N = "Cell motility",
  Z = "Cytoskeleton",
  W = "Extracellular structures",
  U = "Intracellular trafficking, secretion and vesicular transport",
  O = "Post-translational modification, protein turnover and chaperones",
  X = "Mobilome: prophages and transposons",
  C = "Energy production and conversion",
  G = "Carbohydrate transport and metabolism",
  E = "Amino acid transport and metabolism",
  F = "Nucleotide transport and metabolism",
  H = "Coenzyme transport and metabolism",
  I = "Lipid transport and metabolism",
  P = "Inorganic ion transport and metabolism",
  Q = "Secondary metabolite biosynthesis, transport and catabolism",
  R = "General function prediction only",
  S = "Function unknown")

KEGG_DESC <- c(
  ko00010 = "Glycolysis / gluconeogenesis",
  ko00030 = "Pentose phosphate pathway",
  ko00051 = "Fructose and mannose metabolism",
  ko00061 = "Fatty acid biosynthesis",
  ko00130 = "Ubiquinone and other terpenoid-quinone biosynthesis",
  ko00190 = "Oxidative phosphorylation",
  ko00195 = "Photosynthesis",
  ko00196 = "Photosynthesis - antenna proteins",
  ko00220 = "Arginine biosynthesis",
  ko00230 = "Purine metabolism",
  ko00240 = "Pyrimidine metabolism",
  ko00250 = "Alanine, aspartate and glutamate metabolism",
  ko00260 = "Glycine, serine and threonine metabolism",
  ko00261 = "Monobactam biosynthesis",
  ko00270 = "Cysteine and methionine metabolism",
  ko00300 = "Lysine biosynthesis",
  ko00330 = "Arginine and proline metabolism",
  ko00400 = "Phenylalanine, tyrosine and tryptophan biosynthesis",
  ko00430 = "Taurine and hypotaurine metabolism",
  ko00500 = "Starch and sucrose metabolism",
  ko00520 = "Amino sugar and nucleotide sugar metabolism",
  ko00521 = "Streptomycin biosynthesis",
  ko00523 = "Polyketide sugar unit biosynthesis",
  ko00525 = "Acarbose and validamycin biosynthesis",
  ko00550 = "Peptidoglycan biosynthesis",
  ko00564 = "Glycerophospholipid metabolism",
  ko00620 = "Pyruvate metabolism",
  ko00630 = "Glyoxylate and dicarboxylate metabolism",
  ko00660 = "C5-branched dibasic acid metabolism",
  ko00670 = "One carbon pool by folate",
  ko00680 = "Methane metabolism",
  ko00710 = "Carbon fixation by Calvin cycle",
  ko00760 = "Nicotinate and nicotinamide metabolism",
  ko00860 = "Porphyrin metabolism",
  ko00900 = "Terpenoid backbone biosynthesis",
  ko00970 = "Aminoacyl-tRNA biosynthesis",
  ko01055 = "Biosynthesis of vancomycin group antibiotics",
  ko01100 = "Metabolic pathways",
  ko01110 = "Biosynthesis of secondary metabolites",
  ko01120 = "Microbial metabolism in diverse environments",
  ko01130 = "Biosynthesis of antibiotics",
  ko01200 = "Carbon metabolism",
  ko01210 = "2-Oxocarboxylic acid metabolism",
  ko01212 = "Fatty acid metabolism",
  ko01230 = "Biosynthesis of amino acids",
  ko01501 = "beta-Lactam resistance",
  ko01502 = "Vancomycin resistance",
  ko02010 = "ABC transporters",
  ko02024 = "Quorum sensing",
  ko03010 = "Ribosome",
  ko03018 = "RNA degradation",
  ko03030 = "DNA replication",
  ko03060 = "Protein export",
  ko03070 = "Bacterial secretion system",
  ko03420 = "Nucleotide excision repair",
  ko03440 = "Homologous recombination",
  ko04212 = "Longevity regulating pathway - worm")

## Umbrella KEGG pathways that aggregate most of metabolism. Flagged in the
## output and not reported in Table 1.
UMBRELLA <- c("ko01100", "ko01110", "ko01120", "ko01130", "ko01200")

## ---- one-sided Fisher (hypergeometric) enrichment --------------------------
## The universe is every protein-coding gene in the PCC 7942 genome, not only
## the annotated subset, so gene ratios read against the full gene set (x/1362)
## and match the ratios printed in Table 1.
enrich <- function(genes, t2g, source_label, set_label, universe) {
  N <- length(universe); n <- length(genes)
  terms <- unique(t2g$term)
  res <- do.call(rbind, lapply(terms, function(tm) {
    inTerm <- t2g$gene[t2g$term == tm]
    k <- length(intersect(genes, inTerm)); K <- length(inTerm)
    if (k == 0) return(NULL)
    data.frame(homology = set_label, source = source_label, term = tm,
               gene_count = k, set_size = n, term_size = K, universe = N,
               fold = (k / n) / (K / N),
               p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  }))
  if (is.null(res)) return(NULL)
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res[order(res$p_value), ]
}

UNIVERSE <- unique(ann$locus_tag)
cat("universe:", length(UNIVERSE), "protein-coding genes\n")
all <- do.call(rbind, lapply(names(SETS), function(s)
  rbind(enrich(SETS[[s]], cog,  "COG",  s, UNIVERSE),
        enrich(SETS[[s]], kegg, "KEGG", s, UNIVERSE))))
all$gene_ratio  <- paste0(all$gene_count, "/", all$set_size)
all$description <- ifelse(all$source == "COG", COG_DESC[all$term], KEGG_DESC[all$term])
all$description[is.na(all$description)] <- ""
all$umbrella    <- ifelse(all$term %in% UMBRELLA, "yes", "")
all <- all[, c("homology", "source", "term", "description", "gene_ratio", "gene_count",
               "set_size", "term_size", "universe", "fold", "p_value", "q_value", "umbrella")]
write.csv(all, o("table1_enrichment.csv"), row.names = FALSE)

sig <- all[all$p_value < P_CUTOFF, ]
write.csv(sig, o("table1_significant.csv"), row.names = FALSE)
cat("\nsignificant terms at p <", P_CUTOFF, " (q shown for reference):\n")
for (s in names(SETS)) {
  x <- sig[sig$homology == s, ]
  cat(" ", s, ":", nrow(x), "terms\n")
  if (nrow(x)) print(x[, c("source", "term", "description", "gene_ratio",
                           "fold", "p_value", "q_value", "umbrella")], row.names = FALSE)
}
