# common.R
# Shared settings and functions for the R scripts. Sourced by the analysis and
# dataset scripts; not run on its own. All paths are relative to the
# repository root, which run_all.sh sets as the working directory.

suppressPackageStartupMessages({library(igraph); library(readxl)})

if (!file.exists("scripts/common.R"))
  stop("Run from the repository root (the folder holding run_all.sh). Now in: ", getwd())

## Curated inputs, read only.
CUR_CORE <- "curation/core_regulator_curation.csv"   # 38 core regulators: names, classification
CUR_S5   <- "curation/S5_curated_regulators.xlsx"    # species regulator lists, S5-4 ortholog map
CUR_S6   <- "curation/S6_source_annotation.xlsx"     # gene and edge annotation for Dataset S6
CUR_FIG6 <- "curation/figure6_categories.csv"        # names and categories, Figure 6 regulators

TAB <- "outputs/tables"
FIG <- "outputs/figures"
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
o <- function(f) file.path(TAB, f)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

SEEDS <- 1:10

## The four networks. The network of record of each is the unseeded GENIE3
## network, stored as GraphML (directed, locus-tag node ids, node attribute TF).
## The ten seeded runs of each are GENIE3 link lists, pruned here to the same
## number of edges. cutoff = edges kept; lcc = keep the largest weakly
## connected component after pruning (core GRN only: 1,111 edges give 1,100).
NETWORKS <- list(
  core  = list(dir = "networks/core",    cutoff = 1111, lcc = TRUE,
               label = "Core GRN", reg_list = TRUE),
  s7942 = list(dir = "networks/PCC7942", cutoff = 3500, lcc = FALSE,
               label = "S. elongatus PCC 7942", reg_list = FALSE),
  s6803 = list(dir = "networks/PCC6803", cutoff = 3500, lcc = FALSE,
               label = "Synechocystis sp. PCC 6803", reg_list = FALSE),
  s7002 = list(dir = "networks/PCC7002", cutoff = 3500, lcc = FALSE,
               label = "Picosynechococcus sp. PCC 7002", reg_list = FALSE)
)
rep_file  <- function(net) file.path(net$dir, "network_of_record.graphml")
seed_file <- function(net, s) file.path(net$dir, "seeded", sprintf("top%d_seed%02d.csv", net$cutoff, s))
tf_input_file <- function(net) file.path(net$dir, "regulators_GENIE3_input.csv")

## Largest weakly connected component, where the network definition asks for it.
prune <- function(g, lcc) {
  if (!lcc) return(g)
  cp <- components(g, mode = "weak")
  induced_subgraph(g, which(cp$membership == which.max(cp$csize)))
}

## Network of record, read from GraphML.
read_network_of_record <- function(net) {
  f <- rep_file(net)
  if (!file.exists(f)) stop("missing ", f)
  prune(read_graph(f, format = "graphml"), net$lcc)
}

## A seeded run: the GENIE3 link list (regulatoryGene, targetGene, weight),
## top `cutoff` edges by weight, then the component filter. GENIE3 weights rank
## edges for pruning only; every measure is computed on the unweighted network.
## igraph applies an edge attribute named "weight" automatically, so the weight
## is kept under another name.
read_seeded_run <- function(net, s) {
  f <- seed_file(net, s)
  if (!file.exists(f)) stop("missing ", f)
  ll <- read.csv(f, stringsAsFactors = FALSE)
  names(ll)[1:3] <- c("from", "to", "weight")
  ll <- head(ll[order(-ll$weight), ], net$cutoff)
  g <- graph_from_data_frame(ll[, c("from", "to", "weight")], directed = TRUE)
  E(g)$genie3_weight <- E(g)$weight
  g <- delete_edge_attr(g, "weight")
  prune(g, net$lcc)
}

## Node identifiers, whichever attribute the GraphML writer used.
node_names <- function(g) {
  for (a in c("name", "id", "label")) {
    if (a %in% vertex_attr_names(g)) {
      v <- as.character(vertex_attr(g, a))
      if (length(v) == vcount(g) && !all(is.na(v))) return(v)
    }
  }
  stop("no node-name attribute; found: ", paste(vertex_attr_names(g), collapse = ", "))
}

## Edge set as "from to" strings, for overlap between runs.
edge_set <- function(g) {
  el <- as_edgelist(g, names = FALSE); n <- node_names(g)
  paste(n[el[, 1]], n[el[, 2]])
}
## Shared edges as a percentage of the mean size of the two edge sets.
edge_overlap <- function(a, b) 200 * length(intersect(a, b)) / (length(a) + length(b))

## Stress centrality: number of shortest directed paths through each node
## (Brandes-style accumulation).
stress_centrality <- function(g) {
  n <- vcount(g)
  stress <- numeric(n)
  adj <- adjacent_vertices(g, V(g), mode = "out")
  for (s in seq_len(n)) {
    sigma <- numeric(n); sigma[s] <- 1
    dist  <- rep(-1, n);  dist[s]  <- 0
    order_seen <- integer(0)
    preds <- vector("list", n)
    queue <- s
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      order_seen <- c(order_seen, v)
      for (w in as.integer(adj[[v]])) {
        if (dist[w] < 0) { dist[w] <- dist[v] + 1; queue <- c(queue, w) }
        if (dist[w] == dist[v] + 1) {
          sigma[w] <- sigma[w] + sigma[v]
          preds[[w]] <- c(preds[[w]], v)
        }
      }
    }
    zeta <- numeric(n)
    for (w in rev(order_seen)) {
      for (v in preds[[w]]) zeta[v] <- zeta[v] + sigma[v] * (1 + zeta[w] / sigma[w])
    }
    zeta[s] <- 0
    stress <- stress + zeta
  }
  stress
}

## Eigenvector centrality on the undirected network, with a fixed start for the
## iterative solver so repeated runs give identical values. The caller's random
## number state is restored.
eigen_fixed_start <- function(u, seed = 1) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  eigen_centrality(u, directed = FALSE)$vector
}

## Values are rounded to 12 significant digits, so regulators with equal
## network positions tie exactly on every platform and ranks do not depend on
## floating-point noise in the last digits.
SIG <- 12

## Node measures used in the paper. Eigenvector centrality is computed on the
## undirected network: in a directed regulator-to-target network almost every
## edge ends at a node with no outgoing edges, and the directed eigenvector
## gives every regulator zero.
node_measures <- function(g) {
  data.frame(
    locus_tag   = node_names(g),
    degree      = degree(g, mode = "all"),
    outdegree   = degree(g, mode = "out"),
    kcore       = coreness(g, mode = "all"),
    betweenness = signif(betweenness(g, directed = TRUE, weights = NA), SIG),
    eigenvector = signif(eigen_fixed_start(as_undirected(g, mode = "collapse")), SIG),
    stress      = signif(stress_centrality(g), SIG),
    stringsAsFactors = FALSE
  )
}

## Integrated Centrality (Equation 2): five measures, each divided by its
## maximum over the ranked regulators, summed once.
##   local            degree, k-core
##   global           betweenness, stress
##   community-aware  eigenvector
MEASURES <- c("degree", "kcore", "betweenness", "eigenvector", "stress")

normalize_and_IC <- function(df) {
  for (m in MEASURES) {
    mx <- max(df[[m]], na.rm = TRUE)
    df[[paste0(m, "_norm")]] <- if (is.finite(mx) && mx > 0) signif(df[[m]] / mx, SIG) else 0
  }
  df$IC <- signif(rowSums(df[, paste0(MEASURES, "_norm")], na.rm = TRUE), SIG)
  df
}

## Regulators = nodes with at least one outgoing edge.
regulators <- function(g) node_names(g)[degree(g, mode = "out") > 0]

## Core regulator curation. A regulator is stress-coupled when Stress_category
## is set.
core_curation <- function(path = CUR_CORE) {
  if (!file.exists(path)) stop("missing ", path)
  d <- read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  stopifnot(nrow(d) == 38, !any(duplicated(d$locus_tag)), !any(is.na(d$TF_name)))
  sc <- trimws(ifelse(is.na(d$Stress_category), "", d$Stress_category))
  d$Stress_category <- ifelse(sc == "", NA, sc)
  d$Stress_coupled  <- ifelse(sc == "", "N", "Y")
  d
}

## The ranked pool of a network: every node with outgoing edges and, for the
## core GRN, every curated regulator present in the network. Two core
## regulators (SigF1, Crp) keep only an incoming edge after pruning; they are
## part of the 38 and are ranked with the rest.
ranked_pool <- function(g, net) {
  p <- regulators(g)
  if (net$reg_list) p <- union(p, intersect(core_curation()$locus_tag, node_names(g)))
  p
}

## Measures and IC over the ranked pool of one network, sorted by IC.
## Ties keep network order (R's order() is stable).
ranked_measures <- function(g, net) {
  d <- node_measures(g)
  d <- normalize_and_IC(d[d$locus_tag %in% ranked_pool(g, net), ])
  d <- d[order(-d$IC), ]
  d$IC_rank <- seq_len(nrow(d))
  d
}

## Regulator names and stress classification. Core regulators from the core
## curation; species regulators from the curated S5 lists, where a name exists.
## Stress classification is made for the core regulators only.
annotate_regulators <- function(df) {
  d  <- core_curation()
  nm <- setNames(as.character(d$TF_name), d$locus_tag)
  sc <- as.character(d$Stress_category); sc[is.na(sc)] <- ""
  st <- setNames(ifelse(sc == "", "no", "yes"), d$locus_tag)
  stc <- setNames(sc, d$locus_tag)
  for (sh in c("S5-1. Curated TFs PCC 7942", "S5-2. Curated TFs PCC 7002",
               "S5-3. Curated TFs PCC 6803")) {
    s5 <- as.data.frame(read_excel(CUR_S5, sheet = sh, .name_repair = "minimal"))
    add <- setNames(as.character(s5$TF_name), s5$locus_tag)
    add <- add[!is.na(add) & !(names(add) %in% names(nm))]
    nm <- c(nm, add)
  }
  df$TF_name         <- unname(nm[df$locus_tag])
  df$stress_related  <- unname(st[df$locus_tag])
  df$stress_category <- unname(stc[df$locus_tag])
  front <- c("locus_tag", "TF_name", "stress_related", "stress_category")
  df[, c(front, setdiff(names(df), front))]
}

## Core-to-species ortholog map from S5-4, keyed by the PCC 7942 locus tag.
ortholog_map <- function() {
  s54 <- as.data.frame(read_excel(CUR_S5, sheet = "S5-4. Conserved regulators",
                                  .name_repair = "minimal"))
  list("7942" = setNames(s54$PCC7942_locus_tag, s54$PCC7942_locus_tag),
       "7002" = setNames(s54$PCC7002_locus_tag, s54$PCC7942_locus_tag),
       "6803" = setNames(s54$PCC6803_locus_tag, s54$PCC7942_locus_tag))
}

## Sort in C-locale order, so row order does not depend on the system locale.
c_order <- function(...) order(..., method = "radix")
