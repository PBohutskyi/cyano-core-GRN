# 05_centralities.R
# Five centrality measures and Integrated Centrality (IC, Equation 2) for the
# network of record of each of the four networks, with the mean and standard
# deviation of every value across the ten seeded runs.
#
# Inputs : networks/<network>/network_of_record.graphml
#          networks/<network>/seeded/top<N>_seed01..10.csv
#          curation/core_regulator_curation.csv, curation/S5_curated_regulators.xlsx
# Outputs: outputs/tables/centralities_{core,s7942,s6803,s7002}.csv  (Dataset S7)

source("scripts/common.R")

for (nm in names(NETWORKS)) {
  net <- NETWORKS[[nm]]
  g_rep <- read_network_of_record(net)
  say(nm, ": network of record ", vcount(g_rep), " nodes, ", ecount(g_rep), " edges")
  rep_df <- ranked_measures(g_rep, net)
  rep_df <- rep_df[, c("locus_tag", MEASURES, paste0(MEASURES, "_norm"), "IC")]

  seed_tabs <- lapply(SEEDS, function(s) ranked_measures(read_seeded_run(net, s), net))

  ## seed mean and sd per regulator, on the normalized scale; a regulator
  ## absent from a run counts as zero in that run
  cols <- c(paste0(MEASURES, "_norm"), "IC")
  agg <- data.frame(locus_tag = rep_df$locus_tag, stringsAsFactors = FALSE)
  for (cl in cols) {
    M <- sapply(seed_tabs, function(d) d[[cl]][match(agg$locus_tag, d$locus_tag)])
    M <- matrix(as.numeric(M), nrow = nrow(agg))
    M[is.na(M)] <- 0
    agg[[paste0(cl, "_seedmean")]] <- rowMeans(M)
    agg[[paste0(cl, "_seedsd")]]   <- apply(M, 1, sd)
  }

  out <- merge(rep_df, agg, by = "locus_tag", all.x = TRUE)
  out <- out[order(-out$IC), ]
  out$IC_rank <- seq_len(nrow(out))
  ## distance of the network-of-record value from the seeded runs, in seed sd
  out$IC_z_vs_seeds <- ifelse(out$IC_seedsd > 0,
                              round((out$IC - out$IC_seedmean) / out$IC_seedsd, 2), NA)
  out <- annotate_regulators(out)
  write.csv(out, o(paste0("centralities_", nm, ".csv")), row.names = FALSE)
  say(nm, ": ", nrow(out), " regulators ranked")
}
