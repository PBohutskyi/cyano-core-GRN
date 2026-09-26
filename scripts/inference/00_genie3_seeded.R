# 00_genie3_seeded.R
# Step 00: seeded GENIE3 inference for one network, checked against the seeded
# link lists shipped in networks/<network>/seeded/. Run for all four networks
# by `bash run_all.sh --genie3` (about 30 hours on one core), or one network
# at a time:
#
#   Rscript scripts/inference/00_genie3_seeded.R <network> [seeds]
#     <network>  core, PCC7942, PCC6803 or PCC7002
#     [seeds]    optional, e.g. 1 or 1:10 (default 1:10)
#
# Inputs : expression/expression_<network>.csv         (from Zenodo; genes in rows)
#          networks/<network>/regulators_GENIE3_input.csv
# Outputs: outputs/genie3_rerun/<network>/top<N>_seedNN.csv
#          outputs/genie3_rerun/<network>/linklist_seedNN.csv.gz   (all edges)
#          outputs/genie3_rerun/<network>/run_manifest.csv, sessionInfo.txt
#
# Settings: random forests, K = sqrt, 1,000 trees, one core, set.seed(seed)
# before each run. Identical results need the same GENIE3 version (1.34.0);
# the number of cores must stay at 1.

suppressPackageStartupMessages(library(GENIE3))
if (!file.exists("scripts/common.R")) stop("Run from the repository root.")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("core", "PCC7942", "PCC6803", "PCC7002"))
  stop("usage: Rscript scripts/inference/00_genie3_seeded.R <core|PCC7942|PCC6803|PCC7002> [seeds]")
network <- args[1]
seeds   <- if (length(args) > 1) eval(parse(text = args[2])) else 1:10
top_n   <- if (network == "core") 1111 else 3500

expr_file <- file.path("expression", paste0("expression_", network, ".csv"))
tf_file   <- file.path("networks", network, "regulators_GENIE3_input.csv")
ref_dir   <- file.path("networks", network, "seeded")
out_dir   <- file.path("outputs", "genie3_rerun", network)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(expr_file)) stop("missing ", expr_file, " (download from Zenodo, see README)")

say <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), network, " ", ..., "\n", sep = "")

expr <- as.matrix(read.table(expr_file, header = TRUE, row.names = 1, sep = ",",
                             quote = "", stringsAsFactors = FALSE))
storage.mode(expr) <- "numeric"
tfs <- unique(as.character(read.csv(tf_file, stringsAsFactors = FALSE)$locus_tag))
tfs <- tfs[!is.na(tfs) & tfs != ""]
stopifnot(all(tfs %in% rownames(expr)))
say("genes ", nrow(expr), ", conditions ", ncol(expr), ", regulators ", length(tfs),
    ", GENIE3 ", as.character(packageVersion("GENIE3")))

check <- character(0)
for (s in seeds) {
  t0 <- Sys.time()
  set.seed(s)
  wm <- GENIE3(expr, regulators = tfs, treeMethod = "RF", K = "sqrt",
               nTrees = 1000, nCores = 1, verbose = FALSE)
  ll <- getLinkList(wm)
  write.csv(ll, gzfile(file.path(out_dir, sprintf("linklist_seed%02d.csv.gz", s))),
            row.names = FALSE, quote = FALSE)
  f_top <- file.path(out_dir, sprintf("top%d_seed%02d.csv", top_n, s))
  write.csv(head(ll, top_n), f_top, row.names = FALSE, quote = FALSE)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  ## compare with the shipped seeded run: same edges in the same order
  ref <- read.csv(file.path(ref_dir, basename(f_top)), stringsAsFactors = FALSE)
  new <- read.csv(f_top, stringsAsFactors = FALSE)
  same_edges <- identical(paste(ref[[1]], ref[[2]]), paste(new[[1]], new[[2]]))
  max_w <- max(abs(ref[[3]] - new[[3]]))
  res <- sprintf("seed %d: %s, %d edges in full list, top %d edges %s, max weight difference %.2g, %.1f min",
                 s, network, nrow(ll), top_n, if (same_edges) "IDENTICAL" else "DIFFERENT", max_w, mins)
  say(res); check <- c(check, res)
}
writeLines(check, file.path(out_dir, "comparison_with_shipped_runs.txt"))
if (any(grepl("DIFFERENT", check))) stop("a rerun differs from the shipped seeded run")
write.csv(data.frame(
  item  = c("network", "expr_file", "tf_file", "genes", "conditions", "regulators_used",
            "seeds", "nTrees", "nCores", "treeMethod", "K", "top_n", "GENIE3_version", "R_version"),
  value = c(network, expr_file, tf_file, nrow(expr), ncol(expr), length(tfs),
            paste(range(seeds), collapse = "-"), 1000, 1, "RF", "sqrt", top_n,
            as.character(packageVersion("GENIE3")), R.version.string)),
  file.path(out_dir, "run_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
