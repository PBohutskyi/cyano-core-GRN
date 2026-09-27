# 13_figureS2_species_stability.R
# Supplementary Figure S2: reproducibility of the species GRN results across
# GENIE3 runs, the species counterpart of Figure 4. Each species network of
# record is compared with ten seeded runs of the same pipeline on the same
# input, pruned to the same 3,500 edges.
#   A  shared edges: representative vs seeded, and seeded vs seeded
#   B  regulator out-degree, representative vs the seeded mean
#   C  rank agreement per centrality measure, representative vs each seeded run
#   D  how many seeded runs place a regulator in its species top 15
#   E  recovery of the representative top N in the seeded runs, N = 5 to 20
#   F  overlap between the species top 15 and the core top 15. The core network
#      is held at its representative form and only the species network varies,
#      so the variation shown has one source.
#
# Inputs : networks/<species>/network_of_record.graphml
#          networks/<species>/seeded/top3500_seed01..10.csv
#          outputs/tables/centralities_{core,s7942,s6803,s7002}.csv  (05_centralities.R)
#          curation/S5_curated_regulators.xlsx                       (S5-4 ortholog map)
# Outputs: outputs/figures/figS2_panelA..F, figS2_panelD_<species>, figS2_key,
#          figS2_keyA..F, figS2_assembled                             (PDF and PNG)
#          outputs/tables/figS2_*.csv                                 (plotted values)
#          outputs/tables/seed_centralities_species.csv   (per seed)
#          outputs/tables/seed_outdegree_species.csv      (per seed, Dataset S8-2)
#          outputs/tables/figS2_numbers_for_SI.txt

source("scripts/common.R")
TOPN   <- 15
NRANGE <- 5:20
SPECIES <- c("s7942", "s6803", "s7002")
SP_KEY  <- c(s7942 = "7942", s6803 = "6803", s7002 = "7002")
## Panel A y-axis, fixed so the scale does not exaggerate species differences
A_YMIN <- 40; A_YMAX <- 95

## ------------------------------------------------ style (shared by all figure scripts)
FONT_PDF <- "Times"
FONT_PNG <- "serif"
BASE_CEX <- 1.0
LWD      <- 1.1
PCH_CEX  <- 1.05
TCL      <- -0.3
SRT_MEAS <- 20
ERR_COL  <- "#000000"
ERR_LWD  <- 0.6
## species keep their Figure 2 colours throughout the paper
STRAIN_COL <- c(PCC7942 = "#D55E00", PCC7002 = "#0072B2", PCC6803 = "#009E73")
SP_COL <- c(s7942 = unname(STRAIN_COL["PCC7942"]),
            s6803 = unname(STRAIN_COL["PCC6803"]),
            s7002 = unname(STRAIN_COL["PCC7002"]))
SP_LAB  <- c(s7942 = "PCC 7942", s6803 = "PCC 6803", s7002 = "PCC 7002")
SP_LAB2 <- c(s7942 = "PCC\n7942", s6803 = "PCC\n6803", s7002 = "PCC\n7002")
REP_PCH  <- 24   # filled triangle: representative vs seeded
SEED_PCH <- 21   # filled circle:   seeded vs seeded
MEAS_LAB <- c(degree = "Degree", kcore = "k-core", betweenness = "Betweenness",
              eigenvector = "Eigenvector", stress = "Stress")
W_P <- 3.6; H_P <- 3.2

## transparent = TRUE for the key files, so they can be dropped onto a panel in
## PowerPoint without a white box hiding what is underneath
emit <- function(stem, w, h, drawfun, transparent = FALSE) {
  f  <- file.path(FIG, stem)
  bg <- if (transparent) NA else "white"
  pdf(paste0(f, ".pdf"), width = w, height = h, family = FONT_PDF,
      useDingbats = FALSE, bg = bg)
  par(family = FONT_PDF, bg = if (transparent) "transparent" else "white")
  drawfun(); dev.off()
  png(paste0(f, ".png"), width = w, height = h, units = "in", res = 600,
      family = FONT_PNG, bg = if (transparent) "transparent" else "white")
  par(family = FONT_PNG); drawfun(); dev.off()
  cat("wrote ", stem, ".pdf/.png\n", sep = "")
}
axbox <- function() box(lwd = LWD)
xlab2 <- function(s, line = 1.9) mtext(s, side = 1, line = line, font = 2)
ylab2 <- function(s, line = 2.4) mtext(s, side = 2, line = line, font = 2, las = 0)
catlab <- function(at, labels, cex = 0.9, line = 0.75, srt = 0) {
  u <- par("usr"); y <- u[3] - (u[4] - u[3]) * line * 0.055
  if (srt == 0) text(at, y, labels, xpd = NA, adj = c(0.5, 1), cex = cex)
  else          text(at, y, labels, xpd = NA, adj = c(1, 0.6), cex = cex, srt = srt)
}
mline <- function(i, v, col = ERR_COL, w = 0.20) {
  v <- v[is.finite(v)]
  if (!length(v)) return(rep(NA_real_, 3))
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  segments(i, q[1], i, q[3], lwd = ERR_LWD, col = col)
  segments(i - w, q[2], i + w, q[2], lwd = ERR_LWD * 2, col = col)
  segments(i - w / 2, q[1], i + w / 2, q[1], lwd = ERR_LWD, col = col)
  segments(i - w / 2, q[3], i + w / 2, q[3], lwd = ERR_LWD, col = col)
  q
}
jit <- function(i, v, col, pch = SEED_PCH, seed = 1, cex = PCH_CEX, amt = 0.11) {
  v <- v[is.finite(v)]
  if (!length(v)) return(invisible(NULL))
  set.seed(seed)
  points(jitter(rep(i, length(v)), amount = amt), v, pch = pch, cex = cex,
         col = col, bg = adjustcolor(col, alpha.f = 0.35), lwd = 0.8)
}
darker <- function(col, f = 0.62) {
  v <- col2rgb(col) * f
  rgb(v[1, ], v[2, ], v[3, ], maxColorValue = 255)   # vectorised over col
}

## =========================================================================
## 1. Representative and seeded networks, per species
## =========================================================================
g_rep <- list(); e_rep <- list(); seed_tab <- list(); e_seed <- list(); od_rows <- list()
for (nm in SPECIES) {
  n <- NETWORKS[[nm]]
  g_rep[[nm]] <- read_network_of_record(n)
  e_rep[[nm]] <- edge_set(g_rep[[nm]])
  say(SP_LAB[nm], ": representative network ", vcount(g_rep[[nm]]), " nodes, ",
      length(e_rep[[nm]]), " edges")
  st <- list(); es <- list()
  for (s in SEEDS) {
    g <- read_seeded_run(n, s)
    es[[as.character(s)]] <- edge_set(g)
    x <- degree(g, mode = "out")
    od_rows[[paste(nm, s)]] <- data.frame(network = nm, seed = s,
      regulator = node_names(g)[x > 0], out_degree = x[x > 0])
    d <- ranked_measures(g, n); d$seed <- s; d$network <- nm
    st[[as.character(s)]] <- d[, c("network", "seed", "locus_tag", MEASURES,
                                   paste0(MEASURES, "_norm"), "IC", "IC_rank")]
  }
  seed_tab[[nm]] <- st; e_seed[[nm]] <- es
  say(SP_LAB[nm], ": ", length(SEEDS), " seeded runs done (",
      min(sapply(st, nrow)), " to ", max(sapply(st, nrow)), " regulators ranked)")
}
seed_long <- do.call(rbind, lapply(SPECIES, function(nm) do.call(rbind, seed_tab[[nm]])))
write.csv(seed_long, o("seed_centralities_species.csv"), row.names = FALSE)
od_all <- do.call(rbind, od_rows)
write.csv(od_all, o("seed_outdegree_species.csv"), row.names = FALSE)

rep_tab <- lapply(SPECIES, function(nm)
  read.csv(o(sprintf("centralities_%s.csv", nm)), stringsAsFactors = FALSE))
names(rep_tab) <- SPECIES
for (nm in SPECIES) stopifnot(nrow(rep_tab[[nm]]) > 0)

## =========================================================================
## 2. Panel data
## =========================================================================
## ---- A: edge overlap, one definition for both comparisons
A <- do.call(rbind, lapply(SPECIES, function(nm) {
  pr <- combn(length(SEEDS), 2)
  rbind(
    data.frame(network = nm, group = "rep_vs_seed", run_a = "representative",
               run_b = sprintf("seed%d", SEEDS),
               overlap = sapply(e_seed[[nm]], function(x) edge_overlap(e_rep[[nm]], x)),
               stringsAsFactors = FALSE),
    data.frame(network = nm, group = "seed_vs_seed", run_a = sprintf("seed%d", pr[1, ]),
               run_b = sprintf("seed%d", pr[2, ]),
               overlap = apply(pr, 2, function(k)
                 edge_overlap(e_seed[[nm]][[k[1]]], e_seed[[nm]][[k[2]]])),
               stringsAsFactors = FALSE))
}))
write.csv(A, o("figS2_A_edge_overlap.csv"), row.names = FALSE)

## ---- B: regulator out-degree, representative against the seeded mean
B <- do.call(rbind, lapply(SPECIES, function(nm) {
  d  <- od_all[od_all$network == nm, ]
  m  <- tapply(d$out_degree, d$regulator, mean)
  od <- setNames(degree(g_rep[[nm]], mode = "out"), node_names(g_rep[[nm]]))
  x  <- data.frame(network = nm, locus_tag = names(od), rep = as.numeric(od),
                   seed_mean = as.numeric(m[names(od)]), stringsAsFactors = FALSE)
  x[x$rep > 0 & !is.na(x$seed_mean), ]
}))
rho_od <- sapply(SPECIES, function(nm) {
  x <- B[B$network == nm, ]; cor(x$rep, x$seed_mean, method = "spearman") })
write.csv(B, o("figS2_B_outdegree.csv"), row.names = FALSE)

## ---- C: ranking agreement per measure, representative against each seeded run
C <- do.call(rbind, lapply(SPECIES, function(nm) {
  r <- rep_tab[[nm]]
  do.call(rbind, lapply(c(MEASURES, "IC"), function(mm) {
    v <- sapply(SEEDS, function(s) {
      d <- seed_tab[[nm]][[as.character(s)]]
      i <- intersect(d$locus_tag, r$locus_tag)
      a <- r[[mm]][match(i, r$locus_tag)]; b <- d[[mm]][match(i, d$locus_tag)]
      ## a measure with no variation in a run gives no rank correlation
      if (length(i) < 3 || sd(a, na.rm = TRUE) == 0 || sd(b, na.rm = TRUE) == 0) NA_real_
      else cor(a, b, method = "spearman", use = "complete.obs")
    })
    data.frame(network = nm, measure = mm, seed = SEEDS, rho = v,
               stringsAsFactors = FALSE)
  }))
}))
write.csv(C, o("figS2_C_rank_agreement.csv"), row.names = FALSE)

## ---- D: recurrence in the species top TOPN.
## The species networks carry several dozen regulators each, too many for one
## bar each. Exact pool sizes are in outputs/tables/centralities_s*.csv.
## The distribution is shown instead: how many regulators reach the top TOPN in
## exactly k of the ten seeded runs. The per-regulator table is written out and
## the per-species bar charts are exported separately.
topn_seed <- function(nm, n) lapply(seed_tab[[nm]], function(d)
  d$locus_tag[order(d$IC_rank)][seq_len(min(n, nrow(d)))])
rec_tab <- do.call(rbind, lapply(SPECIES, function(nm) {
  r  <- rep_tab[[nm]]
  tb <- table(factor(unlist(topn_seed(nm, TOPN)), levels = r$locus_tag))
  rng <- t(sapply(r$locus_tag, function(gg) {
    v <- sapply(seed_tab[[nm]], function(d) {
      i <- match(gg, d$locus_tag); if (is.na(i)) NA else d$IC_rank[i] })
    c(min = suppressWarnings(min(v, na.rm = TRUE)), med = median(v, na.rm = TRUE),
      max = suppressWarnings(max(v, na.rm = TRUE))) }))
  data.frame(network = nm, locus_tag = r$locus_tag,
             TF_name = if ("TF_name" %in% names(r)) r$TF_name else NA,
             rep_rank = r$IC_rank, seeds_in_top = as.integer(tb[r$locus_tag]),
             rank_min = rng[, "min"], rank_median = rng[, "med"],
             rank_max = rng[, "max"], stringsAsFactors = FALSE)
}))
rec_tab <- rec_tab[order(match(rec_tab$network, SPECIES),
                         -rec_tab$seeds_in_top, rec_tab$rep_rank), ]
write.csv(rec_tab, o("figS2_D_top15_recurrence.csv"), row.names = FALSE)
D <- sapply(SPECIES, function(nm)
  as.integer(table(factor(rec_tab$seeds_in_top[rec_tab$network == nm],
                          levels = 1:length(SEEDS)))))
rownames(D) <- 1:length(SEEDS)

## ---- E: recovery of the representative top N in the seeded runs
Fd <- do.call(rbind, lapply(SPECIES, function(nm) {
  r <- rep_tab[[nm]]
  do.call(rbind, lapply(NRANGE, function(n) {
    rn <- r$locus_tag[order(r$IC_rank)][seq_len(n)]
    f  <- sapply(topn_seed(nm, n), function(x) length(intersect(x, rn)) / n)
    data.frame(network = nm, n = n, mean = mean(f), lo = min(f), hi = max(f),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(Fd, o("figS2_E_topn_recovery.csv"), row.names = FALSE)

## ---- F: overlap between the species top TOPN and the core top TOPN.
## The cross-validation (Figure 5) identifies regulators shared between the
## networks of record. This panel asks whether the sharing survives resampling
## of the species network. The ortholog map is the one used by
## 15_crossvalidation.R.
ortho <- ortholog_map()
core_rep <- read.csv(o("centralities_core.csv"), stringsAsFactors = FALSE)
## core top TOPN, translated into each species' locus tags
core_top <- function(sk) {
  lt <- core_rep$locus_tag[order(core_rep$IC_rank)][seq_len(TOPN)]
  unname(ortho[[sk]][lt])
}
sp_top <- function(nm, sd = NA) {
  if (is.na(sd)) { r <- rep_tab[[nm]]; r$locus_tag[order(r$IC_rank)][seq_len(TOPN)] }
  else { d <- seed_tab[[nm]][[as.character(sd)]]
         d$locus_tag[order(d$IC_rank)][seq_len(TOPN)] }
}
G <- do.call(rbind, lapply(SPECIES, function(nm) {
  sk <- SP_KEY[[nm]]
  rbind(
    data.frame(network = nm, run = "representative",
               shared = length(intersect(sp_top(nm), core_top(sk))),
               stringsAsFactors = FALSE),
    do.call(rbind, lapply(SEEDS, function(sd)
      data.frame(network = nm, run = sprintf("seed%d", sd),
                 shared = length(intersect(sp_top(nm, sd), core_top(sk))),
                 stringsAsFactors = FALSE))))
}))
## Expected overlap if the species top TOPN were drawn at random from all
## regulators of that species network. The draw is from the whole network, not
## from the conserved subset: most species regulators have no core ortholog, and
## using the conserved subset as the denominator would inflate the null.
## Hits available = core top TOPN whose ortholog is present in that network.
G_NULL <- t(sapply(SPECIES, function(nm) {
  sk   <- SP_KEY[[nm]]
  ct   <- core_top(sk)          # already translated into this species' tags
  hits <- sum(!is.na(ct) & ct %in% rep_tab[[nm]]$locus_tag)
  pool <- nrow(rep_tab[[nm]])
  c(hits = hits, pool = pool, expected = TOPN * hits / pool)
}))
rownames(G_NULL) <- SPECIES
G$expected <- G_NULL[G$network, "expected"]
G$p <- mapply(function(nm, k) phyper(k - 1, G_NULL[nm, "hits"],
                G_NULL[nm, "pool"] - G_NULL[nm, "hits"], TOPN, lower.tail = FALSE),
              G$network, G$shared)
write.csv(G, o("figS2_F_core_overlap.csv"), row.names = FALSE)
write.csv(data.frame(network = SPECIES, G_NULL), o("figS2_F_null.csv"), row.names = FALSE)

## =========================================================================
## 3. Panels
## =========================================================================
## A: two jittered columns per species, representative first
panelA <- function() {
  par(mar = c(4.2, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  if (max(A$overlap) > A_YMAX || min(A$overlap) < A_YMIN)
    warning("panel A: values outside [", A_YMIN, ", ", A_YMAX, "] are clipped")
  yl <- c(A_YMIN, A_YMAX)
  plot(NA, xlim = c(0.4, length(SPECIES) + 0.6), ylim = yl, axes = FALSE,
       xlab = "", ylab = "")
  for (k in seq_along(SPECIES)) {
    nm <- SPECIES[k]; col <- SP_COL[nm]
    for (j in 1:2) {
      gp <- c("rep_vs_seed", "seed_vs_seed")[j]
      v  <- A$overlap[A$network == nm & A$group == gp]
      xj <- k + c(-0.18, 0.18)[j]
      jit(xj, v, col, pch = c(REP_PCH, SEED_PCH)[j], seed = k * 10 + j,
          cex = PCH_CEX - 0.15, amt = 0.09)
      mline(xj, v, w = 0.15)
    }
    text(k, yl[2] - 0.03 * (yl[2] - yl[1]), sprintf("%.1f", median(A$overlap[
         A$network == nm & A$group == "seed_vs_seed"])),
         cex = 0.82, col = darker(col), font = 2)
  }
  axis(2, at = seq(A_YMIN, A_YMAX, 10), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = seq_along(SPECIES), labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(seq_along(SPECIES), SP_LAB2[SPECIES], cex = 0.88)
  ylab2("Shared edges (%)")
  xlab2("Species GRN", line = 2.9)
  axbox()
}

## B: out-degree, three species in one panel
panelB <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  lim <- range(c(B$rep, B$seed_mean)); lim[1] <- max(1, lim[1] * 0.8)
  lim[2] <- lim[2] * 1.35
  plot(NA, xlim = lim, ylim = lim, log = "xy", axes = FALSE, xlab = "", ylab = "")
  abline(0, 1, lty = 2, lwd = LWD, col = "#7F7F7F")
  for (nm in SPECIES) {
    x <- B[B$network == nm, ]
    points(x$rep, x$seed_mean, pch = 21, cex = PCH_CEX - 0.2, lwd = 0.7,
           col = SP_COL[nm], bg = adjustcolor(SP_COL[nm], alpha.f = 0.4))
  }
  axis(1, lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, lwd = LWD, cex.axis = 0.95, tcl = TCL)
  ylab2("Out-degree, mean of 10 seeded runs", line = 2.2)
  xlab2("Out-degree, representative network")
  ## legend() does not apply text.font to plotmath entries, so the three rho
  ## labels are drawn with text() instead, which honours font = 2
  u  <- par("usr")
  lx <- 10^(u[1] + 0.03 * (u[2] - u[1]))
  for (k in seq_along(SPECIES)) {
    nm <- SPECIES[k]
    ly <- 10^(u[4] - (0.05 + 0.075 * (k - 1)) * (u[4] - u[3]))
    text(lx, ly, adj = c(0, 0.5), cex = 0.78, font = 2,
         col = darker(SP_COL[nm]),
         labels = bquote(.(SP_LAB[[nm]]) * ":" ~ rho ==
                         .(sprintf("%.3f", rho_od[[nm]]))))
  }
  axbox()
}

## C: rank agreement per measure, species offset within each measure
panelC <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  ms <- c(MEASURES, "IC")
  plot(NA, xlim = c(0.5, length(ms) + 0.5), ylim = c(0, 1.10), axes = FALSE,
       xlab = "", ylab = "")
  abline(h = seq(0, 1, 0.25), col = "#E8E8E8", lwd = 0.7)
  off <- c(-0.24, 0, 0.24)
  for (i in seq_along(ms)) for (k in seq_along(SPECIES)) {
    nm <- SPECIES[k]
    v  <- C$rho[C$network == nm & C$measure == ms[i]]
    jit(i + off[k], v, SP_COL[nm], seed = i * 10 + k, cex = PCH_CEX - 0.3, amt = 0.07)
    mline(i + off[k], v, w = 0.10)
  }
  axis(2, at = seq(0, 1, 0.25), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = seq_along(ms), labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(seq_along(ms), c(unname(MEAS_LAB[MEASURES]), "IC"), cex = 0.78, srt = SRT_MEAS)
  ylab2(expression(bold(paste("Agreement with representative (", rho, ")"))))
  xlab2("Centrality measure", line = 2.5)
  axbox()
}

## D: distribution of recurrence, grouped bars
panelD <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  bp <- barplot(t(D), beside = TRUE, col = SP_COL[SPECIES], border = "black",
                lwd = 0.5, axes = FALSE, names.arg = rep("", nrow(D)),
                ylim = c(0, max(D) * 1.12), space = c(0, 0.55))
  axis(2, lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = colMeans(bp), labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(colMeans(bp), rownames(D), cex = 0.85)
  ylab2("Number of regulators")
  xlab2(paste0("Seeded runs placing a regulator in the top ", TOPN))
  axbox()
}

## per-species version of D, exported separately, not part of the assembly
panelD_one <- function(nm) function() {
  d <- rec_tab[rec_tab$network == nm & rec_tab$seeds_in_top > 0, ]
  d <- d[order(d$seeds_in_top, -d$rep_rank), ]
  lab <- if (all(is.na(d$TF_name))) d$locus_tag else d$TF_name
  lab_cex <- if (nrow(d) > 30) 0.42 else 0.6
  mw <- max(strwidth(lab, units = "inches", cex = lab_cex, font = 2))
  par(mar = c(3.6, (mw + 0.55) / 0.2, 0.6, 0.6), mgp = c(2.4, 0.6, 0),
      cex = BASE_CEX, las = 1)
  bp <- barplot(d$seeds_in_top, horiz = TRUE, col = SP_COL[nm], border = "black",
                lwd = 0.5, xlim = c(0, 10.4), axes = FALSE,
                names.arg = rep("", nrow(d)), space = 0.35)
  axis(1, at = seq(0, 10, 2), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, at = bp, labels = lab, tick = FALSE, line = -0.15,
       cex.axis = lab_cex, font = 2)
  ylab2("Regulator", line = mw / 0.2 + 0.9)
  xlab2(paste0("Seeded runs placing it in the top ", TOPN, ", ", SP_LAB[nm]))
  abline(v = length(SEEDS), lty = 3, lwd = LWD, col = "#7F7F7F")
  axbox()
}

## E: recovery curves, one per species
panelE <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  plot(NA, xlim = range(NRANGE), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  abline(h = seq(0, 1, 0.25), col = "#E8E8E8", lwd = 0.7)
  for (nm in SPECIES) {
    x <- Fd[Fd$network == nm, ]
    polygon(c(x$n, rev(x$n)), c(x$lo, rev(x$hi)),
            col = adjustcolor(SP_COL[nm], alpha.f = 0.16),
            border = adjustcolor(SP_COL[nm], alpha.f = 0.55), lwd = 0.4)
  }
  for (nm in SPECIES) {
    x <- Fd[Fd$network == nm, ]
    lines(x$n, x$mean, lwd = 1.8, col = SP_COL[nm])
    points(x$n, x$mean, pch = 21, cex = PCH_CEX - 0.35, col = SP_COL[nm],
           bg = adjustcolor(SP_COL[nm], alpha.f = 0.45), lwd = 0.7)
  }
  abline(v = TOPN, lty = 3, lwd = LWD, col = "#7F7F7F")
  axis(1, at = seq(min(NRANGE), max(NRANGE), 5), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, at = seq(0, 1, 0.25), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  ylab2("Representative top N recovered")
  xlab2("Cutoff N")
  axbox()
}

## F: overlap with the core top 15. Box and whisker rather than a jittered
## column, so the panel reads differently from Figure 4 and stays compact with
## three species. Whiskers span the full range of the seeded runs, not a Tukey
## rule, because there are only ten values per species and nothing is hidden.
panelF <- function() {
  par(mar = c(4.2, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  vs <- lapply(SPECIES, function(nm) G$shared[G$network == nm & G$run != "representative"])
  vr <- sapply(SPECIES, function(nm) G$shared[G$network == nm & G$run == "representative"])
  yl <- c(0, min(8, max(c(unlist(vs), vr, G_NULL[, "expected"])) + 1.6))
  plot(NA, xlim = c(0.4, length(SPECIES) + 0.6), ylim = yl, axes = FALSE,
       xlab = "", ylab = "")
  for (k in seq_along(SPECIES)) {
    nm <- SPECIES[k]; col <- SP_COL[nm]; v <- vs[[k]]
    q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
    segments(k, min(v), k, max(v), lwd = ERR_LWD, col = ERR_COL)      # full range
    segments(c(k - 0.09, k - 0.09), c(min(v), max(v)),
             c(k + 0.09, k + 0.09), c(min(v), max(v)), lwd = ERR_LWD, col = ERR_COL)
    rect(k - 0.22, q[1], k + 0.22, q[3], col = adjustcolor(col, alpha.f = 0.45),
         border = ERR_COL, lwd = ERR_LWD)
    segments(k - 0.22, q[2], k + 0.22, q[2], lwd = ERR_LWD * 2.4, col = ERR_COL)
    points(k, vr[k], pch = REP_PCH, cex = PCH_CEX, col = ERR_COL,
           bg = darker(col), lwd = 0.7)                                # representative
    ## null for this species: the denominator differs between networks, so the
    ## expected value is drawn per species rather than as one line
    segments(k - 0.30, G_NULL[k, "expected"], k + 0.30, G_NULL[k, "expected"],
             lty = 2, lwd = LWD * 2, col = "#7F7F7F")
  }
  axis(2, at = seq(0, floor(yl[2]), 2), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  if (max(c(unlist(vs), vr)) > yl[2])
    warning("panel F: values above ", yl[2], " are clipped")
  axis(1, at = seq_along(SPECIES), labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(seq_along(SPECIES), SP_LAB2[SPECIES], cex = 0.88)
  ylab2(paste0("Shared with the core top ", TOPN, " (IC)"))
  xlab2("Species GRN", line = 2.9)
  axbox()
}

## Per-panel keys, each covering only what its own panel uses, exported as
## separate transparent files. The combined key below is also written.
g_ <- "#8C8C8C"
## pch follows the panel: circles where the panel draws points, squares where
## it draws bars or boxes
key_species <- function(where = "topleft", title = "Species network", pch = 22) {
  legend(where, bty = "n", cex = 0.9, y.intersp = 1.3,
         pt.cex = if (pch == 21) PCH_CEX else 1.5,
         title = title, title.adj = 0, title.font = 2,
         legend = unname(SP_LAB[SPECIES]), pch = pch,
         col = if (pch == 21) SP_COL[SPECIES] else "black",
         pt.bg = if (pch == 21) adjustcolor(SP_COL[SPECIES], alpha.f = 0.4)
                 else SP_COL[SPECIES], pt.lwd = 0.8)
}
key_A <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species(pch = 21)
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         legend = c("Representative vs seeded run", "Seeded run vs seeded run",
                    "Median and interquartile range"),
         pch = c(REP_PCH, SEED_PCH, NA), lty = c(NA, NA, 1),
         lwd = c(NA, NA, ERR_LWD * 2.4), pt.cex = c(PCH_CEX, PCH_CEX, NA),
         col = c(g_, g_, ERR_COL),
         pt.bg = c(adjustcolor(g_, alpha.f = 0.4), adjustcolor(g_, alpha.f = 0.4), NA))
}
key_B <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species(pch = 21)
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         legend = c("Regulator", "1:1 line"), pch = c(21, NA), lty = c(NA, 2),
         lwd = c(NA, LWD), pt.cex = c(PCH_CEX, NA), col = c(g_, g_),
         pt.bg = c(adjustcolor(g_, alpha.f = 0.4), NA))
}
key_C <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species(pch = 21)
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         legend = c("One seeded run", "Median and interquartile range"),
         pch = c(SEED_PCH, NA), lty = c(NA, 1), lwd = c(NA, ERR_LWD * 2.4),
         pt.cex = c(PCH_CEX - 0.3, NA), col = c(g_, ERR_COL),
         pt.bg = c(adjustcolor(g_, alpha.f = 0.4), NA))
}
key_D <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species()
}
key_E <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species(pch = 21)
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         legend = c("Mean across seeded runs", "Range across seeded runs",
                    paste0("Cutoff reported in the main text (N = ", TOPN, ")")),
         pch = c(21, 22, NA), lty = c(1, NA, 3), lwd = c(1.8, NA, LWD),
         pt.cex = c(PCH_CEX - 0.35, 1.6, NA),
         col = c(g_, adjustcolor(g_, alpha.f = 0.55), g_),
         pt.bg = c(adjustcolor(g_, alpha.f = 0.45),
                   adjustcolor(g_, alpha.f = 0.16), NA))
}
key_F <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); key_species()
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         legend = c("Representative species network",
                    "Median and interquartile range across seeded runs",
                    "Full range across seeded runs",
                    "Expected overlap at random"),
         pch = c(REP_PCH, 22, NA, NA), lty = c(NA, NA, 1, 2),
         lwd = c(NA, NA, ERR_LWD, LWD * 2),
         pt.cex = c(PCH_CEX, 1.6, NA, NA),
         col = c(ERR_COL, ERR_COL, ERR_COL, "#7F7F7F"),
         pt.bg = c(g_, adjustcolor(g_, alpha.f = 0.45), NA, NA))
}

## combined key, added to the assembled figure in PowerPoint
key_panel <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new()
  g <- "#8C8C8C"
  legend("topleft", bty = "n", cex = 0.92, y.intersp = 1.3, pt.cex = 1.5,
         title = "Species network", title.adj = 0, title.font = 2,
         legend = unname(SP_LAB[SPECIES]), pch = 22,
         col = "black", pt.bg = SP_COL[SPECIES], pt.lwd = 0.6)
  legend("topright", bty = "n", cex = 0.86, y.intersp = 1.3,
         title = "Symbols and lines", title.adj = 0, title.font = 2,
         legend = c("Representative network",
                    "Seeded run",
                    "Median and interquartile range",
                    "Full range across seeded runs",
                    "Mean across seeded runs (panel E)",
                    "Range across seeded runs, shaded (panel E)",
                    "1:1 line (panel B)",
                    paste0("Cutoff reported in the main text (N = ", TOPN, ")"),
                    "Expected overlap at random, per species (panel F)"),
         pch    = c(REP_PCH, SEED_PCH, NA, NA, 21, 22, NA, NA, NA),
         lty    = c(NA, NA, 1, 1, 1, NA, 2, 3, 2),
         lwd    = c(NA, NA, ERR_LWD * 2.4, ERR_LWD, 1.8, NA, LWD, LWD, LWD * 2),
         pt.cex = c(PCH_CEX, PCH_CEX - 0.15, NA, NA, PCH_CEX - 0.35, 1.6, NA, NA, NA),
         col    = c(ERR_COL, g, ERR_COL, ERR_COL, g,
                    adjustcolor(g, alpha.f = 0.55), g, g, g),
         pt.bg  = c(g, adjustcolor(g, alpha.f = 0.4), NA, NA,
                    adjustcolor(g, alpha.f = 0.45),
                    adjustcolor(g, alpha.f = 0.16), NA, NA, NA))
}

## =========================================================================
## 4. Write panels, key and the assembled figure
## =========================================================================
emit("figS2_panelA", W_P, H_P, panelA)
emit("figS2_panelB", W_P, H_P, panelB)
emit("figS2_panelC", W_P, H_P, panelC)
emit("figS2_panelD", W_P, H_P, panelD)
emit("figS2_panelE", W_P, H_P, panelE)
emit("figS2_panelF", W_P, H_P, panelF)
emit("figS2_key",    6.4,  2.9, key_panel, transparent = TRUE)
for (kk in list(c("A", "key_A"), c("B", "key_B"), c("C", "key_C"),
                c("D", "key_D"), c("E", "key_E"), c("F", "key_F")))
  emit(paste0("figS2_key", kk[1]), 4.6, 1.7, get(kk[2]), transparent = TRUE)
for (nm in SPECIES)
  emit(paste0("figS2_panelD_", sub("^s", "", nm)), W_P * 1.2, H_P * 2.2, panelD_one(nm))
emit("figS2_assembled", 3 * W_P, 2 * H_P, function() {
  layout(matrix(1:6, 2, 3, byrow = TRUE))
  panelA(); panelB(); panelC()
  panelD(); panelE(); panelF()
})

## =========================================================================
## 5. Numbers the SI needs
## =========================================================================
con <- file(o("figS2_numbers_for_SI.txt"), "w")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
qs <- function(v) sprintf("%.1f%% (IQR %.1f to %.1f)", median(v),
                          quantile(v, 0.25), quantile(v, 0.75))
w("Supplementary Figure S2 source numbers (13_figureS2_species_stability.R)")
for (nm in SPECIES)
  w("  ", SP_LAB[nm], ": representative network ", vcount(g_rep[[nm]]), " nodes, ",
    length(e_rep[[nm]]), " edges, ", nrow(rep_tab[[nm]]), " regulators ranked")
w("")
w("A. Shared edges between runs:")
for (nm in SPECIES) {
  w("  ", SP_LAB[nm], " seeded vs seeded: ",
    qs(A$overlap[A$network == nm & A$group == "seed_vs_seed"]))
  w("  ", SP_LAB[nm], " representative vs seeded: ",
    qs(A$overlap[A$network == nm & A$group == "rep_vs_seed"]))
}
w("")
w("B. Regulator out-degree, representative vs the seeded mean, Spearman rho:")
for (nm in SPECIES) w("  ", SP_LAB[nm], ": ", sprintf("%.3f", rho_od[[nm]]),
                      " over ", sum(B$network == nm), " regulators")
w("")
w("C. Rank correlation with the representative network, median over ",
  length(SEEDS), " seeded runs:")
for (nm in SPECIES) {
  w("  ", SP_LAB[nm], ":")
  for (mm in c(MEASURES, "IC"))
    w("    ", ifelse(mm == "IC", "IC", MEAS_LAB[mm]), ": ",
      sprintf("%.3f", median(C$rho[C$network == nm & C$measure == mm], na.rm = TRUE)),
      " (min ", sprintf("%.3f", min(C$rho[C$network == nm & C$measure == mm], na.rm = TRUE)),
      ")", if (any(is.na(C$rho[C$network == nm & C$measure == mm])))
             sprintf(" [%d of %d runs undefined]",
                     sum(is.na(C$rho[C$network == nm & C$measure == mm])), length(SEEDS)) else "")
}
w("")
w("D. Regulators reaching the top ", TOPN, " by IC:")
for (nm in SPECIES) {
  r <- rec_tab[rec_tab$network == nm, ]
  w("  ", SP_LAB[nm], ": ", sum(r$seeds_in_top == length(SEEDS)), " in all ",
    length(SEEDS), " seeded runs, ", sum(r$seeds_in_top >= 8), " in at least 8, ",
    sum(r$seeds_in_top > 0), " in at least one")
}
w("")
w("E. Recovery of the representative top N in the seeded runs:")
for (nm in SPECIES) {
  x <- Fd[Fd$network == nm, ]
  w("  ", SP_LAB[nm], ": N = 10 ", sprintf("%.0f%%", 100 * x$mean[x$n == 10]),
    ", N = ", TOPN, " ", sprintf("%.0f%%", 100 * x$mean[x$n == TOPN]),
    ", N = 20 ", sprintf("%.0f%%", 100 * x$mean[x$n == 20]))
}
w("")
w("F. Regulators shared between the species top ", TOPN, " and the core top ", TOPN,
  " (core network held fixed at its representative form):")
for (nm in SPECIES) {
  v <- G$shared[G$network == nm & G$run != "representative"]
  w("  ", SP_LAB[nm], ": representative ",
    G$shared[G$network == nm & G$run == "representative"], " of ", TOPN,
    "; seeded runs median ", median(v), ", range ", min(v), " to ", max(v))
}
for (nm in SPECIES) {
  v <- G$shared[G$network == nm & G$run != "representative"]
  w("  ", SP_LAB[nm], " expected at random ",
    sprintf("%.1f", G_NULL[nm, "expected"]), " of ", TOPN,
    " (", G_NULL[nm, "hits"], " of the core top ", TOPN,
    " have an ortholog among the ", G_NULL[nm, "pool"], " regulators of that network)",
    "; fold over random ", sprintf("%.1f", median(v) / G_NULL[nm, "expected"]),
    ", hypergeometric p at the median ",
    signif(phyper(median(v) - 1, G_NULL[nm, "hits"],
                  G_NULL[nm, "pool"] - G_NULL[nm, "hits"], TOPN,
                  lower.tail = FALSE), 3))
}
close(con)

