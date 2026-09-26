# 11_figure4_core_stability.R
# Figure 4: reproducibility of the core GRN result across GENIE3 runs. The core
# network of record is compared with ten seeded runs of the same pipeline on
# the same input. All seeded networks are pruned exactly as the network of
# record (top 1,111 edges, largest connected component).
#   A  shared edges: network of record vs each seeded run, and seeded vs seeded
#   B  regulator out-degree, network of record vs the seeded mean
#   C  rank agreement per centrality measure, network of record vs each run
#   D  how many seeded runs place each regulator in the top 15 by IC
#   E  stress-coupled regulators in each seeded top 15, against the null
#   F  recovery of the network-of-record top N in the seeded runs, N = 5 to 20
#
# Inputs : networks/core/network_of_record.graphml
#          networks/core/seeded/top1111_seed01..10.csv
#          outputs/tables/centralities_core.csv   (05_centralities.R)
# Outputs: outputs/figures/fig4_panelA..F, fig4_panelB_labelled, fig4_key,
#          fig4_keyB, fig4_keyF, fig4_assembled          (PDF and PNG)
#          outputs/tables/fig4_*.csv                     (plotted values)
#          outputs/tables/seed_centralities_core.csv     (per seed, Dataset S8-3)
#          outputs/tables/seed_outdegree_core.csv        (per seed, Dataset S8-2)
#          outputs/tables/fig4_numbers_for_SI.txt

source("scripts/common.R")
TOPN   <- 15                # cutoff used for the headline enrichment test
NRANGE <- 5:20              # cutoffs scanned in panel F

## ------------------------------------------------ style (shared by all figure scripts)
FONT_PDF <- "Times"; FONT_PNG <- "serif"; BASE_CEX <- 1.0; LWD <- 1.1
## one accent for the representative network, one for the seeded runs, and the
## stress-coupled pair; all drawn from the Figure 2 colourblind-safe palette
REP_COL    <- "#000000"
SEED_COL   <- "#0072B2"
STRESS_COL <- c(yes = "#CC79A7", no = "#BFBFBF")
MEAS_COL   <- c(degree = "#56B4E9", kcore = "#009E73", betweenness = "#E69F00",
                eigenvector = "#CC79A7", stress = "#7F7F7F")
MEAS_LAB   <- c(degree = "Degree", kcore = "k-core", betweenness = "Betweenness",
                eigenvector = "Eigenvector", stress = "Stress")
PCH_CEX  <- 1.05            # symbol size for the jittered distributions
TCL      <- -0.3            # tick length, same on every axis of every panel
SRT_MEAS <- 20              # rotation of the centrality labels in panel C (degrees)
A_YMIN   <- 50; A_YMAX <- 85   # panel A y-axis range
ERR_COL  <- "#000000"       # median and interquartile bars
ERR_LWD  <- 0.6             # thin, so the points stay readable
W_P <- 3.6; H_P <- 3.2      # one panel of the 3 x 2 grid

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
## axis titles in bold
xlab2 <- function(s, line = 1.9) mtext(s, side = 1, line = line, font = 2)
ylab2 <- function(s, line = 2.4) mtext(s, side = 2, line = line, font = 2, las = 0)
## base R silently drops categorical labels that would overlap, so they are
## drawn with text() instead and every category stays labelled
catlab <- function(at, labels, cex = 0.9, line = 0.75, srt = 0) {
  u <- par("usr"); y <- u[3] - (u[4] - u[3]) * line * 0.055
  if (srt == 0) text(at, y, labels, xpd = NA, adj = c(0.5, 1), cex = cex)
  else          text(at, y, labels, xpd = NA, adj = c(1, 0.6), cex = cex, srt = srt)
}
## median and interquartile range, drawn over a jittered column
mline <- function(i, v, col = ERR_COL, w = 0.28) {
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  segments(i, q[1], i, q[3], lwd = ERR_LWD, col = col)
  segments(i - w, q[2], i + w, q[2], lwd = ERR_LWD * 2, col = col)
  segments(i - w / 2, q[1], i + w / 2, q[1], lwd = ERR_LWD, col = col)
  segments(i - w / 2, q[3], i + w / 2, q[3], lwd = ERR_LWD, col = col)
  q
}
jit <- function(i, v, col, seed = 1, cex = PCH_CEX) {
  set.seed(seed)
  points(jitter(rep(i, length(v)), amount = 0.15), v, pch = 21, cex = cex,
         col = col, bg = adjustcolor(col, alpha.f = 0.35), lwd = 0.8)
}
## a darker version of a fill colour, for text that must stay legible
darker <- function(col, f = 0.62) {
  v <- col2rgb(col) * f
  rgb(v[1], v[2], v[3], maxColorValue = 255)
}
fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

## =========================================================================
## 1. Network of record and seeded core networks
## =========================================================================
net <- NETWORKS$core
g_rep <- read_network_of_record(net)
nm_rep <- node_names(g_rep)
e_rep <- edge_set(g_rep)
say("core network of record: ", vcount(g_rep), " nodes, ", length(e_rep), " edges")

g_seed <- lapply(SEEDS, function(s) read_seeded_run(net, s))
seed_tab <- lapply(SEEDS, function(s) {
  d <- ranked_measures(g_seed[[s]], net); d$seed <- s
  d[, c("seed", "locus_tag", MEASURES, paste0(MEASURES, "_norm"), "IC", "IC_rank")]
})
names(g_seed) <- names(seed_tab) <- as.character(SEEDS)
seed_long <- annotate_regulators(do.call(rbind, seed_tab))
write.csv(seed_long, o("seed_centralities_core.csv"), row.names = FALSE)

rep_core <- read.csv(o("centralities_core.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(rep_core) == 38)
is_stress <- setNames(rep_core$stress_related == "yes", rep_core$locus_tag)
tfname    <- setNames(rep_core$TF_name, rep_core$locus_tag)
N <- nrow(rep_core); K <- sum(is_stress)
say("pool ", N, " regulators, ", K, " stress-coupled")

## =========================================================================
## 2. Panel data
## =========================================================================
## ---- A: edge overlap. Both comparisons use one definition, the shared
## fraction of the two pruned edge sets, so the two columns are on one scale.
e_seed <- lapply(g_seed, edge_set)
rep_vs_seed <- sapply(e_seed, function(x) edge_overlap(e_rep, x))
pr <- combn(length(SEEDS), 2)
seed_vs_seed <- apply(pr, 2, function(k) edge_overlap(e_seed[[k[1]]], e_seed[[k[2]]]))
A <- rbind(data.frame(group = "rep_vs_seed", run_a = "representative",
                      run_b = sprintf("seed%d", SEEDS), overlap = rep_vs_seed),
           data.frame(group = "seed_vs_seed", run_a = sprintf("seed%d", pr[1, ]),
                      run_b = sprintf("seed%d", pr[2, ]), overlap = seed_vs_seed))
wA <- suppressWarnings(wilcox.test(rep_vs_seed, seed_vs_seed))
write.csv(A, o("fig4_A_edge_overlap.csv"), row.names = FALSE)

## ---- B: regulator out-degree, representative against the seeded mean
## out-degree of every regulator in every pruned seeded run
od <- do.call(rbind, lapply(SEEDS, function(s) {
  g <- g_seed[[s]]; x <- degree(g, mode = "out")
  data.frame(seed = s, regulator = node_names(g)[x > 0], out_degree = x[x > 0])
}))
write.csv(od, o("seed_outdegree_core.csv"), row.names = FALSE)
od_mean <- tapply(od$out_degree, od$regulator, mean)
od_rep  <- setNames(degree(g_rep, mode = "out"), nm_rep)
B <- data.frame(locus_tag = names(od_rep), rep = as.numeric(od_rep),
                seed_mean = as.numeric(od_mean[names(od_rep)]), stringsAsFactors = FALSE)
B <- B[B$rep > 0 & !is.na(B$seed_mean), ]
B$TF_name <- unname(tfname[B$locus_tag])
B$stress  <- ifelse(is_stress[B$locus_tag] %in% TRUE, "yes", "no")
rho_od <- cor(B$rep, B$seed_mean, method = "spearman")
ct_od  <- suppressWarnings(cor.test(B$rep, B$seed_mean, method = "spearman"))
sl_od  <- unname(coef(lm(log10(B$seed_mean) ~ log10(B$rep)))[2])   # slope on log-log
write.csv(B, o("fig4_B_outdegree.csv"), row.names = FALSE)

## ---- C: ranking agreement per measure, representative against each seeded run
rho_vs <- function(col) sapply(SEEDS, function(s) {
  d <- seed_tab[[as.character(s)]]
  i <- intersect(d$locus_tag, rep_core$locus_tag)
  cor(rep_core[[col]][match(i, rep_core$locus_tag)], d[[col]][match(i, d$locus_tag)],
      method = "spearman")
})
C <- do.call(rbind, lapply(c(MEASURES, "IC"), function(m)
  data.frame(measure = m, seed = SEEDS, rho = rho_vs(m), stringsAsFactors = FALSE)))
write.csv(C, o("fig4_C_rank_agreement.csv"), row.names = FALSE)

## ---- D: how many seeded runs place each regulator in the top TOPN by IC
topn_seed <- function(n) lapply(seed_tab, function(d)
  d$locus_tag[order(d$IC_rank)][seq_len(min(n, nrow(d)))])
top_by_seed <- topn_seed(TOPN)
rec <- table(factor(unlist(top_by_seed), levels = rep_core$locus_tag))
rng <- t(sapply(rep_core$locus_tag, function(g) {
  r <- sapply(seed_tab, function(d) { i <- match(g, d$locus_tag)
                                      if (is.na(i)) NA else d$IC_rank[i] })
  c(min = suppressWarnings(min(r, na.rm = TRUE)), med = median(r, na.rm = TRUE),
    max = suppressWarnings(max(r, na.rm = TRUE))) }))
D <- data.frame(locus_tag = rep_core$locus_tag,
                TF_name = unname(tfname[rep_core$locus_tag]),
                rep_rank = rep_core$IC_rank,
                seeds_in_top = as.integer(rec[rep_core$locus_tag]),
                rank_min = rng[, "min"], rank_median = rng[, "med"], rank_max = rng[, "max"],
                stress = ifelse(is_stress[rep_core$locus_tag], "yes", "no"),
                stringsAsFactors = FALSE)
D <- D[order(-D$seeds_in_top, D$rep_rank), ]
write.csv(D, o("fig4_D_top15_recurrence.csv"), row.names = FALSE)

## ---- E: stress-coupled regulators in each seeded top TOPN, against the null
E <- data.frame(seed = SEEDS,
                k_stress = sapply(top_by_seed, function(x) sum(is_stress[x], na.rm = TRUE)),
                stringsAsFactors = FALSE)
E$p <- phyper(E$k_stress - 1, K, N - K, TOPN, lower.tail = FALSE)
EXP   <- TOPN * K / N
k_rep <- sum(is_stress[rep_core$locus_tag[rep_core$IC_rank <= TOPN]])
p_rep <- phyper(k_rep - 1, K, N - K, TOPN, lower.tail = FALSE)
write.csv(E, o("fig4_E_stress_by_seed.csv"), row.names = FALSE)

## ---- F: recovery of the representative top N in the seeded runs, over a range
## of N. Shows the agreement is not specific to the cutoff chosen.
Fd <- do.call(rbind, lapply(NRANGE, function(n) {
  rep_n <- rep_core$locus_tag[order(rep_core$IC_rank)][seq_len(n)]
  f <- sapply(topn_seed(n), function(x) length(intersect(x, rep_n)) / n)
  data.frame(n = n, mean = mean(f), lo = min(f), hi = max(f), stringsAsFactors = FALSE)
}))
write.csv(Fd, o("fig4_F_topn_recovery.csv"), row.names = FALSE)

## =========================================================================
## 3. Panels
## =========================================================================
GRP_LAB <- c("Representative\nvs 10 seeded", "10 seeded\nvs 10 seeded")

panelA <- function() {
  par(mar = c(4.2, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  if (max(A$overlap) > A_YMAX)
    warning("panel A: values above A_YMAX (", A_YMAX, ") are clipped; raise it")
  plot(NA, xlim = c(0.5, 2.5), ylim = c(A_YMIN, A_YMAX), axes = FALSE,
       xlab = "", ylab = "")
  cols <- c(REP_COL, SEED_COL)
  for (i in 1:2) {
    v <- A$overlap[A$group == c("rep_vs_seed", "seed_vs_seed")[i]]
    jit(i, v, cols[i], seed = i)
    q <- mline(i, v)
    text(i, A_YMAX - 0.05 * (A_YMAX - A_YMIN), sprintf("%.1f", q[2]), cex = 0.88, col = darker(cols[i]), font = 2)
  }
  axis(2, at = seq(A_YMIN, A_YMAX, 5), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = 1:2, labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(1:2, GRP_LAB, cex = 0.85)
  ylab2("Shared edges (%)")
  xlab2("Core GRN comparison", line = 2.7)
  text(1.5, A_YMIN + 0.04 * (A_YMAX - A_YMIN), fmt_p(wA$p.value), cex = 0.8, col = "#4D4D4D")
  axbox()
}

## panel B, optionally with regulator names. The labelled version is a reading
## aid: labels are placed by a simple rule and may overlap, so the unlabelled
## version is the one to assemble from.
panelB <- function(labels = FALSE) {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  lim <- range(c(B$rep, B$seed_mean)); lim[1] <- max(1, lim[1] * 0.8); lim[2] <- lim[2] * 1.35
  plot(B$rep, B$seed_mean, log = "xy", xlim = lim, ylim = lim, axes = FALSE,
       xlab = "", ylab = "", pch = 21, cex = PCH_CEX + 0.1, lwd = 0.9,
       col = STRESS_COL[B$stress], bg = adjustcolor(STRESS_COL[B$stress], alpha.f = 0.45))
  abline(0, 1, lty = 2, lwd = LWD, col = "#7F7F7F")
  axis(1, lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, lwd = LWD, cex.axis = 0.95, tcl = TCL)
  ylab2("Out-degree, mean of 10 seeded runs", line = 2.2)
  xlab2("Out-degree, representative network")
  legend("topleft", bty = "n", cex = 0.8, adj = 0, y.intersp = 1.15,
         legend = c(as.expression(bquote(rho == .(sprintf("%.3f", rho_od)) *
                                         "," ~ .(fmt_p(ct_od$p.value)))),
                    as.expression(bquote("slope" == .(sprintf("%.2f", sl_od))))))
  if (labels) {
    side <- ifelse(seq_len(nrow(B)) %% 2 == 0, 1, -1)
    lx <- B$rep * (1 + 0.30 * side); ly <- B$seed_mean * (1 - 0.17 * side)
    segments(B$rep, B$seed_mean, lx, ly, lwd = 0.5, col = "#9A9A9A")
    text(lx, ly, B$TF_name, cex = 0.5, adj = ifelse(side > 0, 0, 1))
  }
  axbox()
}

panelC <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  ms <- c(MEASURES, "IC"); cols <- c(MEAS_COL[MEASURES], IC = REP_COL)
  plot(NA, xlim = c(0.5, length(ms) + 0.5), ylim = c(0, 1.10), axes = FALSE,
       xlab = "", ylab = "")
  abline(h = seq(0, 1, 0.25), col = "#E8E8E8", lwd = 0.7)
  for (i in seq_along(ms)) {
    v <- C$rho[C$measure == ms[i]]
    jit(i, v, cols[i], seed = i, cex = PCH_CEX - 0.1)
    q <- mline(i, v, w = 0.26)
    text(i, 1.07, sprintf("%.2f", q[2]), cex = 0.72, col = darker(cols[i]), font = 2)
  }
  axis(2, at = seq(0, 1, 0.25), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = seq_along(ms), labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(seq_along(ms), c(unname(MEAS_LAB[MEASURES]), "IC"), cex = 0.78, srt = SRT_MEAS)
  ylab2(expression(bold(paste("Agreement with representative (", rho, ")"))))
  xlab2("Centrality measure", line = 2.5)
  axbox()
}

panelD <- function() {
  d <- D[D$seeds_in_top > 0, ]
  d <- d[order(d$seeds_in_top, -d$rep_rank), ]
  ## left margin sized to the longest regulator name, so no space is wasted
  lab_cex <- 0.66
  mw <- max(strwidth(d$TF_name, units = "inches", cex = lab_cex, font = 2))
  par(mar = c(3.6, (mw + 0.55) / 0.2, 0.6, 0.6), mgp = c(2.4, 0.6, 0),
      cex = BASE_CEX, las = 1)
  bp <- barplot(d$seeds_in_top, horiz = TRUE, col = STRESS_COL[d$stress],
                border = "black", lwd = 0.7, xlim = c(0, 10.4), axes = FALSE,
                names.arg = rep("", nrow(d)), space = 0.35)
  axis(1, at = seq(0, 10, 2), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, at = bp, labels = d$TF_name, tick = FALSE, line = -0.15,
       cex.axis = lab_cex, font = 2)
  ylab2("Regulator", line = mw / 0.2 + 0.9)
  xlab2(paste0("Seeded runs placing it in the top ", TOPN))
  abline(v = length(SEEDS), lty = 3, lwd = LWD, col = "#7F7F7F")
  axbox()
}

panelE <- function(annot = TRUE) {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  lo <- max(0, min(floor(EXP) - 1, min(E$k_stress) - 1))
  hi <- min(TOPN, max(ceiling(EXP) + 1, max(E$k_stress) + 1))
  x <- lo:hi; y <- as.integer(table(factor(E$k_stress, levels = x)))
  bp <- barplot(y, col = adjustcolor(STRESS_COL["yes"], alpha.f = 0.55),
                border = "black", lwd = 0.7, ylim = c(0, max(y) + 2.1), axes = FALSE,
                names.arg = rep("", length(x)), space = 0.3)
  axis(2, at = 0:max(y + 1), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(1, at = bp, labels = FALSE, lwd = LWD, tcl = TCL)
  catlab(bp, x, cex = 0.88)
  ylab2("Number of seeded runs")
  xlab2(paste0("Stress-coupled regulators in the top ", TOPN))
  xat <- function(v) approx(x, bp, v, rule = 2)$y
  abline(v = xat(EXP),   lty = 2, lwd = 1.6, col = "#7F7F7F")
  abline(v = xat(k_rep), lty = 1, lwd = 1.6, col = REP_COL)
  if (annot) {
    text(xat(EXP), max(y) + 1.75, sprintf("expected %.1f", EXP), cex = 0.78,
         col = "#7F7F7F", adj = c(0, 0.5))
    text(xat(k_rep), max(y) + 0.85, "representative", cex = 0.78, adj = c(1.08, 0.5))
  }
  axbox()
}

panelF <- function() {
  par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.5, 0.7, 0), cex = BASE_CEX, las = 1)
  plot(NA, xlim = range(NRANGE), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  abline(h = seq(0, 1, 0.25), col = "#E8E8E8", lwd = 0.7)
  polygon(c(Fd$n, rev(Fd$n)), c(Fd$lo, rev(Fd$hi)),
          col = adjustcolor(SEED_COL, alpha.f = 0.20),
          border = adjustcolor(SEED_COL, alpha.f = 0.75), lwd = 0.5)
  lines(Fd$n, Fd$mean, lwd = 2.0, col = SEED_COL)
  points(Fd$n, Fd$mean, pch = 21, cex = PCH_CEX - 0.15, col = SEED_COL,
         bg = adjustcolor(SEED_COL, alpha.f = 0.45), lwd = 0.8)
  abline(v = TOPN, lty = 3, lwd = LWD, col = "#7F7F7F")
  axis(1, at = seq(min(NRANGE), max(NRANGE), 5), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  axis(2, at = seq(0, 1, 0.25), lwd = LWD, cex.axis = 0.95, tcl = TCL)
  ylab2("Representative top N recovered")
  xlab2("Cutoff N")
  axbox()
}

## per-panel keys, exported separately for assembly in PowerPoint
key_B <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new()
  legend("center", bty = "n", cex = 0.9, y.intersp = 1.3, pt.cex = PCH_CEX + 0.1,
         legend = c("Stress-coupled regulator", "Not stress-coupled", "1:1 line"),
         pch = c(21, 21, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, LWD),
         col = c(STRESS_COL["yes"], STRESS_COL["no"], "#7F7F7F"),
         pt.bg = c(adjustcolor(STRESS_COL["yes"], alpha.f = 0.45),
                   adjustcolor(STRESS_COL["no"],  alpha.f = 0.45), NA))
}
key_F <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new()
  legend("center", bty = "n", cex = 0.9, y.intersp = 1.3,
         legend = c("Mean of 10 seeded runs", "Range across the 10 seeded runs",
                    paste0("Cutoff reported in the main text (N = ", TOPN, ")")),
         pch = c(21, 22, NA), lty = c(1, NA, 3), lwd = c(2.0, NA, LWD),
         pt.cex = c(PCH_CEX - 0.15, 1.8, NA),
         col = c(SEED_COL, adjustcolor(SEED_COL, alpha.f = 0.75), "#7F7F7F"),
         pt.bg = c(adjustcolor(SEED_COL, alpha.f = 0.45),
                   adjustcolor(SEED_COL, alpha.f = 0.20), NA))
}

key_panel <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new()
  legend("center",
         legend = c("Stress-coupled regulator", "Not stress-coupled",
                    "Representative network", "Seeded runs",
                    "Expected under a random draw", "Median and interquartile range"),
         fill   = c(STRESS_COL["yes"], STRESS_COL["no"], NA, NA, NA, NA),
         border = c("black", "black", NA, NA, NA, NA),
         lty = c(NA, NA, 1, 1, 2, 1), lwd = c(NA, NA, 1.6, 2.0, 1.6, 2.3),
         col = c(NA, NA, REP_COL, SEED_COL, "#7F7F7F", "#4D4D4D"),
         bty = "n", cex = 0.9, y.intersp = 1.25)
}

## =========================================================================
## 4. Write panels, key and the assembled figure
## =========================================================================
## Standalone panels omit the in-plot annotations of panel E, which are added
## in PowerPoint. The assembled version carries them.
emit("fig4_panelA", W_P, H_P, panelA)
emit("fig4_panelB", W_P, H_P, function() panelB(labels = FALSE))
emit("fig4_panelB_labelled", W_P * 1.5, H_P * 1.5, function() panelB(labels = TRUE))
emit("fig4_panelC", W_P, H_P, panelC)
emit("fig4_panelD", W_P, H_P, panelD)
emit("fig4_panelE", W_P, H_P, function() panelE(annot = FALSE))
emit("fig4_panelF", W_P, H_P, panelF)
emit("fig4_key",    3.0,  1.9, key_panel, transparent = TRUE)
emit("fig4_keyB",   2.6,  1.2, key_B, transparent = TRUE)
emit("fig4_keyF",   3.2,  1.2, key_F, transparent = TRUE)
emit("fig4_assembled", 3 * W_P, 2 * H_P, function() {
  layout(matrix(1:6, 2, 3, byrow = TRUE))
  panelA(); panelB(labels = FALSE); panelC()
  panelD(); panelE(annot = TRUE);   panelF()
})

## =========================================================================
## 5. Numbers the SI and the figure legend need
## =========================================================================
con <- file(o("fig4_numbers_for_SI.txt"), "w")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
qs <- function(v) sprintf("%.1f%% (IQR %.1f to %.1f)", median(v),
                          quantile(v, 0.25), quantile(v, 0.75))
w("Figure 4 source numbers, core GRN (11_figure4_core_stability.R)")
w("Representative network ", vcount(g_rep), " nodes, ", length(e_rep), " edges; ",
  length(SEEDS), " seeded runs of the same pipeline")
w("")
w("A. Shared edges, representative vs each seeded run: ", qs(rep_vs_seed))
w("   seeded vs seeded (45 pairs): ", qs(seed_vs_seed))
w("   Wilcoxon rank sum, ", fmt_p(wA$p.value))
w("")
w("B. Regulator out-degree, representative vs the seeded mean: Spearman rho = ",
  sprintf("%.3f", rho_od), ", ", fmt_p(ct_od$p.value), ", n = ", nrow(B),
  "; log-log slope ", sprintf("%.2f", sl_od))
w("")
w("C. Rank correlation with the representative network, median over ",
  length(SEEDS), " seeded runs:")
for (m in c(MEASURES, "IC"))
  w("   ", ifelse(m == "IC", "IC", MEAS_LAB[m]), ": ",
    sprintf("%.3f", median(C$rho[C$measure == m])),
    " (min ", sprintf("%.3f", min(C$rho[C$measure == m])), ")")
w("")
w("D. Top-", TOPN, " recurrence: ", sum(D$seeds_in_top == length(SEEDS)),
  " regulators appear in all ", length(SEEDS), " seeded runs, ",
  sum(D$seeds_in_top >= 8), " in at least 8, ", sum(D$seeds_in_top > 0), " in at least one")
w("   of those appearing in at least 8, ",
  sum(D$seeds_in_top >= 8 & D$stress == "yes"), " are stress-coupled")
w("")
w("E. Stress-coupled regulators in the top ", TOPN, ":")
w("   representative ", k_rep, " of ", TOPN, ", p = ", signif(p_rep, 3))
w("   seeded runs median ", median(E$k_stress), ", range ", min(E$k_stress),
  " to ", max(E$k_stress), "; ", sum(E$p < 0.05), " of ", length(SEEDS),
  " significant at 0.05")
w("   expected under a random draw ", sprintf("%.1f", EXP),
  " (pool ", N, ", ", K, " stress-coupled)")
w("")
w("F. Recovery of the representative top N in the seeded runs:")
for (n in c(10, TOPN, 20))
  w("   N = ", n, ": mean ", sprintf("%.0f%%", 100 * Fd$mean[Fd$n == n]),
    " (range ", sprintf("%.0f", 100 * Fd$lo[Fd$n == n]), " to ",
    sprintf("%.0f", 100 * Fd$hi[Fd$n == n]), ")")
close(con)

