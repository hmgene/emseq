#outlier Y1, E2
library(data.table)
tt=fread("../../bigdata/2025-08-01/filtered_3rep_3x_per_group.csv.gz") 
thr <- 10  # minimum total counts per group

## remove outliers
tt[, grep("Y1|E2", names(tt),value=T) := NULL ]
for (g in c("E", "W", "Y")) {
  tt[, paste0("t", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\..?CpG"), names(tt), value = TRUE)]
  tt[, paste0("c", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\.CpG"), names(tt), value = TRUE)]
  tt[, paste0("p", g) := get(paste0("c", g)) / get(paste0("t", g))]
}

res <- tt[abs(pW - pE) > 0.1 | abs(pY-pW), {
  # Initialize results as numeric
  pvEW <- NA_real_; pvWY <- NA_real_; l2EW <- NA_real_; l2WY <- NA_real_;
  # EW comparison
  if (tE > 0 & tW > 0) {
    pvEW <- prop.test(x = c(cE, cW), n = c(tE, tW))$p.value
    l2EW <- log2((cW / tW) / (cE / tE))  # log2 fold change
  }
  # WY comparison
  if (tW > 0 & tY > 0) {
    pvWY <- prop.test(x = c(cW, cY), n = c(tW, tY))$p.value
    l2WY <- log2((cY / tY) / (cW / tW))  # log2 fold change
  }
  .(pvEW = pvEW, l2EW = l2EW, pvWY = pvWY, l2WY = l2WY)
}, by = .(chrom, start)]


library(data.table)
library(VennDiagram)

# Assume res is your data.table with columns: pvEW, l2EW, pvWY, l2WY

# Define significant positive/negative for each comparison
sig_EW_pos <- res[pvEW < 0.01 & l2EW > 0, paste(chrom, start)]
sig_EW_neg <- res[pvEW < 0.01 & l2EW < 0, paste(chrom, start)]

sig_WY_pos <- res[pvWY < 0.01 & l2WY > 0, paste(chrom, start)]
sig_WY_neg <- res[pvWY < 0.01 & l2WY < 0, paste(chrom, start)]

# Example: Venn diagram for positive changes
venn.plot <- venn.diagram(
  x = list(EW = sig_EW_pos, WY = sig_WY_pos),
  category.names = c("EW up", "WY up"),
  filename = NULL,  # NULL to plot to R device
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Significant Positive Changes"
)

grid.draw(venn.plot)

# Example: Venn diagram for negative changes
venn.plot_neg <- venn.diagram(
  x = list(EW = sig_EW_neg, WY = sig_WY_neg),
  category.names = c("EW down", "WY down"),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Significant Negative Changes"
)

grid.draw(venn.plot_neg)

