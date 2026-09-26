# 01_regulator_counts.R
# Ortholog groups across the three genomes, and the regulator counts shown in
# Figure 2: candidates per prediction source (2A), curated regulators per
# source (2C), and the three-way overlaps of both sets (2B, 2D).
#
# Inputs : datasets/Dataset_S3_Core_genome_homologs.xlsx   ortholog groups
#          datasets/Dataset_S4_Initial_putative_TFs.xlsx   candidates and sources
#          curation/S5_curated_regulators.xlsx             curated regulator lists
# Outputs: outputs/tables/ortholog_groups.csv
#          outputs/tables/figure2_pipeline_counts.csv      Figure 2A and 2C
#          outputs/tables/figure2_tf_venn.csv              Figure 2B and 2D

source("scripts/common.R")
S3 <- "datasets/Dataset_S3_Core_genome_homologs.xlsx"
S4 <- "datasets/Dataset_S4_Initial_putative_TFs.xlsx"
STRAINS <- c("PCC7942", "PCC7002", "PCC6803")
blank <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA", "#N/A", "None")
rd <- function(f, sheet) as.data.frame(read_excel(f, sheet = sheet, .name_repair = "minimal"))

## ---------------------------------------------------------------- ortholog groups
# Each row of the tri-homolog and duo-homolog sheets is one orthologous group.
grp <- list()
take <- function(sheet, sps) {
  d <- rd(S3, sheet)
  for (i in seq_len(nrow(d))) {
    g <- list()
    for (sp in sps) {
      col <- paste0(sp, "_locus_tag")
      if (col %in% names(d) && !blank(d[i, col])) g[[sp]] <- trimws(as.character(d[i, col]))
    }
    if (length(g)) grp[[length(grp) + 1]] <<- g
  }
}
take("S3-2. Tri-homologous core",    STRAINS)
take("S3-3. Duo-homologs 7942-7002", c("PCC7942", "PCC7002"))
take("S3-4. Duo-homologs 7942-6803", c("PCC7942", "PCC6803"))
take("S3-5. Duo-homologs 6803-7002", c("PCC6803", "PCC7002"))
cat("ortholog groups:", length(grp), "\n")
write.csv(do.call(rbind, lapply(grp, function(g)
  data.frame(PCC7942 = ifelse(is.null(g$PCC7942), NA, g$PCC7942),
             PCC7002 = ifelse(is.null(g$PCC7002), NA, g$PCC7002),
             PCC6803 = ifelse(is.null(g$PCC6803), NA, g$PCC6803)))),
  o("ortholog_groups.csv"), row.names = FALSE)

## ---------------------------------------------------------------- Figure 2A and 2C
s4 <- list(PCC7942 = rd(S4, "S4-1. Putative TFs PCC 7942"),
           PCC7002 = rd(S4, "S4-2. Putative TFs PCC 7002"),
           PCC6803 = rd(S4, "S4-3. Putative TFs PCC 6803"))
s5 <- list(PCC7942 = rd(CUR_S5, "S5-1. Curated TFs PCC 7942"),
           PCC7002 = rd(CUR_S5, "S5-2. Curated TFs PCC 7002"),
           PCC6803 = rd(CUR_S5, "S5-3. Curated TFs PCC 6803"))
METH <- c("DeepTF", "ENTRAF", "P2TF", "NCBI")
bars <- do.call(rbind, lapply(STRAINS, function(sp) {
  d <- s4[[sp]]; keep <- d$locus_tag %in% s5[[sp]]$locus_tag
  rbind(
    data.frame(strain = sp, stage = "initial", source = "Total",
               count = nrow(d), stringsAsFactors = FALSE),
    data.frame(strain = sp, stage = "initial", source = METH,
               count = sapply(METH, function(m) sum(d[[m]] == "Y"))),
    data.frame(strain = sp, stage = "curated", source = "Total",
               count = nrow(s5[[sp]])),
    data.frame(strain = sp, stage = "curated", source = METH,
               count = sapply(METH, function(m) sum(d[[m]][keep] == "Y"))))
}))
write.csv(bars, o("figure2_pipeline_counts.csv"), row.names = FALSE)
print(bars)

## ---------------------------------------------------------------- Figure 2B and 2D
# Counted by orthologous group, so circle totals equal the set sizes.
venn <- function(sets, label) {
  used <- setNames(vector("list", 3), STRAINS)
  cnt  <- list()
  for (g in grp) {
    hit <- STRAINS[sapply(STRAINS, function(sp) !is.null(g[[sp]]) && g[[sp]] %in% sets[[sp]])]
    if (length(hit) >= 2) {
      k <- paste(hit, collapse = "&")
      cnt[[k]] <- (if (is.null(cnt[[k]])) 0 else cnt[[k]]) + 1
      for (sp in hit) used[[sp]] <- c(used[[sp]], g[[sp]])
    }
  }
  for (sp in STRAINS) cnt[[sp]] <- length(setdiff(sets[[sp]], used[[sp]]))
  regions <- c(STRAINS, "PCC7942&PCC7002", "PCC7942&PCC6803", "PCC7002&PCC6803",
               "PCC7942&PCC7002&PCC6803")
  out <- data.frame(panel = label, region = regions,
                    count = sapply(regions, function(r) if (is.null(cnt[[r]])) 0 else cnt[[r]]))
  tot <- sapply(STRAINS, function(sp) sum(out$count[grepl(sp, out$region, fixed = TRUE)]))
  stopifnot(all(tot == sapply(sets, length)))
  out
}
f2b <- venn(lapply(s4, function(d) d$locus_tag), "2B initial")
f2d <- venn(lapply(s5, function(d) d$locus_tag), "2D curated")
write.csv(rbind(f2b, f2d), o("figure2_tf_venn.csv"), row.names = FALSE)
print(rbind(f2b, f2d))
