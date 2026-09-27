# 14_network_recurrence.R
# How much of each representative network recurs in its ten seeded runs:
#   edge recurrence  - for each edge of the representative network, the number
#                      of seeded runs that contain the same regulator-target edge
#   node overlap     - genes shared between the representative and each seeded
#                      network, as a percentage of the mean node count
#   regulator overlap- the same for regulators (nodes with outgoing edges)
# Seeded runs are pruned exactly as the network of record (common.R).
#
# Inputs : networks/<network>/network_of_record.graphml, networks/<network>/seeded/
# Outputs: outputs/tables/edge_recurrence.csv            one row per representative edge
#          outputs/tables/network_recurrence_summary.csv one row per network (Dataset S8-6)

source("scripts/common.R")

edges_long <- list(); summ <- list()
for (key in names(NETWORKS)) {
  net <- NETWORKS[[key]]
  rep <- read_network_of_record(net)
  e_rep <- edge_set(rep); n_rep <- node_names(rep); r_rep <- regulators(rep)
  runs <- lapply(SEEDS, function(s) read_seeded_run(net, s))
  e_seed <- lapply(runs, edge_set)

  hits <- Reduce(`+`, lapply(e_seed, function(x) as.integer(e_rep %in% x)))
  edges_long[[key]] <- data.frame(network = net$label, edge = e_rep, seeded_runs_with_edge = hits,
                                  stringsAsFactors = FALSE)
  pct <- function(a, b) 200 * length(intersect(a, b)) / (length(a) + length(b))
  node_ov <- sapply(runs, function(g) pct(n_rep, node_names(g)))
  reg_ov  <- sapply(runs, function(g) pct(r_rep, regulators(g)))

  summ[[key]] <- data.frame(
    network = net$label,
    representative_edges = length(e_rep),
    edges_in_at_least_1_run_percent  = round(100 * mean(hits >= 1), 1),
    edges_in_at_least_5_runs_percent = round(100 * mean(hits >= 5), 1),
    edges_in_all_10_runs_percent     = round(100 * mean(hits == length(SEEDS)), 1),
    node_overlap_median_percent      = round(median(node_ov), 1),
    regulator_overlap_median_percent = round(median(reg_ov), 1),
    stringsAsFactors = FALSE)
  say(net$label, ": ", summ[[key]]$edges_in_at_least_1_run_percent, "% of edges recur in >= 1 run, ",
      summ[[key]]$edges_in_at_least_5_runs_percent, "% in >= 5 runs; regulator overlap ",
      summ[[key]]$regulator_overlap_median_percent, "%")
}
write.csv(do.call(rbind, edges_long), o("edge_recurrence.csv"), row.names = FALSE)
S <- do.call(rbind, summ); rownames(S) <- NULL
write.csv(S, o("network_recurrence_summary.csv"), row.names = FALSE)
print(S, row.names = FALSE)
