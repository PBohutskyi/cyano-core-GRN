# 02_figure1A.R
# Figure 1A: gene-level three-way Venn diagram of the three genomes.
# Shared regions come from Dataset S3; genome totals are the numbers of
# protein-coding genes in reference/genome_gene_lists.csv.
#
# Inputs : datasets/Dataset_S3_Core_genome_homologs.xlsx
#          reference/genome_gene_lists.csv
# Outputs: outputs/tables/figure1A_regions.csv
#          outputs/figures/fig1_panelA.{pdf,png}

source("scripts/common.R")
S3 <- "datasets/Dataset_S3_Core_genome_homologs.xlsx"

## ---- style (shared by all figure scripts) -----------------------------------
FONT_PDF <- "Times"; FONT_PNG <- "serif"; LWD <- 1.1
STRAIN_COL <- c(PCC7942 = "#D55E00", PCC7002 = "#0072B2", PCC6803 = "#009E73")
W_VENN <- 3.9; H_VENN <- 3.9
emit <- function(stem, w, h, drawfun) {
  f <- file.path(FIG, stem)
  pdf(paste0(f, ".pdf"), width = w, height = h, family = FONT_PDF, useDingbats = FALSE)
  par(family = FONT_PDF); drawfun(); dev.off()
  png(paste0(f, ".png"), width = w, height = h, units = "in", res = 600, family = FONT_PNG)
  par(family = FONT_PNG); drawfun(); dev.off()
  cat("wrote ", stem, ".pdf/.png\n", sep = "")
}

## ---- counts ----------------------------------------------------------------
nrows <- function(sh) nrow(as.data.frame(read_excel(S3, sheet = sh, .name_repair = "minimal")))
tri  <- nrows("S3-2. Tri-homologous core")
d_ab <- nrows("S3-3. Duo-homologs 7942-7002")
d_ac <- nrows("S3-4. Duo-homologs 7942-6803")
d_bc <- nrows("S3-5. Duo-homologs 6803-7002")

lists <- read.csv("reference/genome_gene_lists.csv", stringsAsFactors = FALSE)
tot <- sapply(c("PCC7942", "PCC6803", "PCC7002"),
              function(s) sum(!is.na(lists[[s]]) & trimws(lists[[s]]) != ""))

only_a <- tot[["PCC7942"]] - (tri + d_ab + d_ac)
only_b <- tot[["PCC7002"]] - (tri + d_ab + d_bc)
only_c <- tot[["PCC6803"]] - (tri + d_ac + d_bc)
stopifnot(only_a >= 0, only_b >= 0, only_c >= 0)

reg <- data.frame(
  region = c("PCC7942", "PCC7002", "PCC6803", "PCC7942&PCC7002",
             "PCC7942&PCC6803", "PCC7002&PCC6803", "PCC7942&PCC7002&PCC6803"),
  count  = c(only_a, only_b, only_c, d_ab, d_ac, d_bc, tri))
write.csv(reg, o("figure1A_regions.csv"), row.names = FALSE)
print(reg)
cat("\ngenome totals:", paste(names(tot), tot, collapse = "  "), "\n")

## ---- Venn, circle area proportional to genome size --------------------------
venn_panel <- function() {
  g <- function(r) reg$count[reg$region == r]
  n <- c(tot[["PCC7942"]], tot[["PCC7002"]], tot[["PCC6803"]])
  r <- 0.95 * sqrt(n / mean(n)); D <- 0.60 * mean(r)
  ang <- c(140, 40, 270) * pi / 180
  cx <- D * cos(ang) * 1.55; cy <- D * sin(ang) * 1.15
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  lim <- max(abs(c(cx, cy)) + r) + 0.30
  plot(NA, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1, axes = FALSE, xlab = "", ylab = "")
  th <- seq(0, 2 * pi, length.out = 512)
  for (i in 1:3)
    polygon(cx[i] + r[i] * cos(th), cy[i] + r[i] * sin(th),
            col = adjustcolor(STRAIN_COL[i], alpha.f = 0.33), border = "black", lwd = LWD)
  gx <- seq(-lim, lim, length.out = 700)
  G  <- expand.grid(x = gx, y = gx)
  inC <- sapply(1:3, function(i) (G$x - cx[i])^2 + (G$y - cy[i])^2 <= r[i]^2)
  code <- apply(inC, 1, function(z) paste(as.integer(z), collapse = ""))
  place <- function(bits, v) { k <- code == bits
    if (any(k)) text(mean(G$x[k]), mean(G$y[k]), format(v, big.mark = ","), cex = 1.15, font = 2) }
  place("100", g("PCC7942")); place("010", g("PCC7002")); place("001", g("PCC6803"))
  place("110", g("PCC7942&PCC7002")); place("101", g("PCC7942&PCC6803"))
  place("011", g("PCC7002&PCC6803")); place("111", g("PCC7942&PCC7002&PCC6803"))
  nm <- c("PCC 7942", "PCC 7002", "PCC 6803")
  for (i in 1:3)
    text(cx[i] * 1.42, cy[i] + ifelse(i == 3, -r[i] - 0.22, r[i] + 0.22),
         nm[i], cex = 1.05, font = 2, col = STRAIN_COL[i])
}
emit("fig1_panelA", W_VENN, H_VENN, venn_panel)
