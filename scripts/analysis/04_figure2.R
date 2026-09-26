# 04_figure2.R
# Figure 2: candidate and curated regulators per prediction source (panels A
# and C) and their overlap across the three strains (panels B and D). Each
# panel is also written on its own, with a separate key.
#
# Inputs : outputs/tables/figure2_pipeline_counts.csv   (01_regulator_counts.R)
#          outputs/tables/figure2_tf_venn.csv           (01_regulator_counts.R)
# Outputs: outputs/figures/fig2_panelA to fig2_panelD, fig2_key, fig2_assembled

source("scripts/common.R")

## ---- style (shared by all figure scripts) -----------------------------------
FONT_PDF <- "Times"; FONT_PNG <- "serif"; BASE_CEX <- 1.0; LWD <- 1.1
COLS <- c(Total = "#7F7F7F", DeepTF = "#56B4E9", ENTRAF = "#009E73",
          P2TF  = "#E69F00", NCBI   = "#CC79A7")
STRAIN_COL <- c(PCC7942 = "#D55E00", PCC7002 = "#0072B2", PCC6803 = "#009E73")
W_BAR <- 4.4; H_BAR <- 2.3; W_VENN <- 3.9; H_VENN <- 3.9

emit <- function(stem, w, h, drawfun) {
  f <- file.path(FIG, stem)
  pdf(paste0(f, ".pdf"), width = w, height = h, family = FONT_PDF, useDingbats = FALSE)
  par(family = FONT_PDF); drawfun(); dev.off()
  png(paste0(f, ".png"), width = w, height = h, units = "in", res = 600, family = FONT_PNG)
  par(family = FONT_PNG); drawfun(); dev.off()
  cat("wrote ", stem, ".pdf/.png\n", sep = "")
}

bars <- read.csv(o("figure2_pipeline_counts.csv"), stringsAsFactors = FALSE)
venn <- read.csv(o("figure2_tf_venn.csv"), stringsAsFactors = FALSE)
STR <- c("PCC7942", "PCC7002", "PCC6803")
LAB <- c("PCC 7942", "PCC 7002", "PCC 6803")
SRC <- c("Total", "DeepTF", "ENTRAF", "P2TF", "NCBI")
KEY <- c("Total TFs", "DeepTF", "ENTRAF", "P2TF", "NCBI")

## ------------------------------------------------ panels A and C: grouped bars
bar_panel <- function(stage, ymax, ystep) {
  m <- sapply(STR, function(s) sapply(SRC, function(k)
    bars$count[bars$strain == s & bars$stage == stage & bars$source == k]))
  par(mar = c(2.9, 3.0, 0.5, 0.5), mgp = c(1.8, 0.45, 0), cex = BASE_CEX, las = 1, tcl = -0.3)
  bp <- barplot(m, beside = TRUE, col = COLS[SRC], border = "black", lwd = LWD,
                ylim = c(0, ymax), axes = FALSE, names.arg = rep("", 3), space = c(0, 1.1))
  axis(2, at = seq(0, ymax, by = ystep), lwd = LWD, cex.axis = 0.95)
  mtext("Count", side = 2, line = 1.9, las = 0)
  axis(1, at = colMeans(bp), labels = LAB, tick = FALSE, line = 0.0, cex.axis = 1.0)
  mtext("Strain", side = 1, line = 1.8)
  box(lwd = LWD)
}

## ------------------------------------------------ standalone key
key_panel <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  legend("center", legend = KEY, fill = COLS[SRC], border = "black",
         bty = "n", cex = 1.1, y.intersp = 1.25)
}

## ------------------------------------------------ panels B and D: area-scaled Venn
## Circle radius is proportional to the square root of the set size, so circle
## AREA is proportional to the number of regulators. Region labels are placed at
## the centroid of each region, found by sampling, so they follow the geometry.
venn_panel <- function(panel) {
  g <- function(r) venn$count[venn$panel == panel & venn$region == r]
  n <- c(g("PCC7942") + g("PCC7942&PCC7002") + g("PCC7942&PCC6803") + g("PCC7942&PCC7002&PCC6803"),
         g("PCC7002") + g("PCC7942&PCC7002") + g("PCC7002&PCC6803") + g("PCC7942&PCC7002&PCC6803"),
         g("PCC6803") + g("PCC7942&PCC6803") + g("PCC7002&PCC6803") + g("PCC7942&PCC7002&PCC6803"))
  r  <- 0.95 * sqrt(n / mean(n))
  D  <- 0.60 * mean(r)                       # centre offset from origin
  ang <- c(140, 40, 270) * pi / 180
  cx <- D * cos(ang) * 1.55; cy <- D * sin(ang) * 1.15

  par(mar = c(0.2, 0.2, 0.2, 0.2))
  lim <- max(abs(c(cx, cy)) + r) + 0.30
  plot(NA, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1, axes = FALSE, xlab = "", ylab = "")
  th <- seq(0, 2 * pi, length.out = 512)
  for (i in 1:3)
    polygon(cx[i] + r[i] * cos(th), cy[i] + r[i] * sin(th),
            col = adjustcolor(STRAIN_COL[i], alpha.f = 0.33), border = "black", lwd = LWD)

  # centroid of each region, by sampling
  gx <- seq(-lim, lim, length.out = 700); gy <- gx
  G  <- expand.grid(x = gx, y = gy)
  inC <- sapply(1:3, function(i) (G$x - cx[i])^2 + (G$y - cy[i])^2 <= r[i]^2)
  code <- apply(inC, 1, function(z) paste(as.integer(z), collapse = ""))
  place <- function(bits, value) {
    k <- code == bits
    if (!any(k)) return(invisible())
    text(mean(G$x[k]), mean(G$y[k]), value, cex = 1.2, font = 2)
  }
  place("100", g("PCC7942")); place("010", g("PCC7002")); place("001", g("PCC6803"))
  place("110", g("PCC7942&PCC7002")); place("101", g("PCC7942&PCC6803"))
  place("011", g("PCC7002&PCC6803")); place("111", g("PCC7942&PCC7002&PCC6803"))

  for (i in 1:3) {
    lx <- cx[i] * 1.42; ly <- cy[i] + ifelse(i == 3, -r[i] - 0.22, r[i] + 0.22)
    text(lx, ly, gsub("PCC", "PCC ", STR[i]), cex = 1.05, font = 2, col = STRAIN_COL[i])
  }
}

## ------------------------------------------------ write everything
emit("fig2_panelA", W_BAR,  H_BAR,  function() bar_panel("initial", 225, 50))
emit("fig2_panelB", W_VENN, H_VENN, function() venn_panel("2B initial"))
emit("fig2_panelC", W_BAR,  H_BAR,  function() bar_panel("curated",  90, 20))
emit("fig2_panelD", W_VENN, H_VENN, function() venn_panel("2D curated"))
emit("fig2_key",    1.7,    1.7,    key_panel)
emit("fig2_assembled", W_BAR + W_VENN, 2 * H_BAR + 0.4, function() {
  layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.15, 1))
  bar_panel("initial", 225, 50); venn_panel("2B initial")
  bar_panel("curated",  90, 20); venn_panel("2D curated")
})
