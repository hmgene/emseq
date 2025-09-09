
binom=function(tt,delta=0.1){
    #outlier Y1, E2
    library(data.table)
    #tt=fread("../../bigdata/2025-08-01/filtered_3rep_3x_per_group.csv.gz") 

    ## remove outliers
    tt[, grep("Y1|E2", names(tt),value=T) := NULL ]
    for (g in c("E", "W", "Y")) {
      tt[, paste0("t", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\.uCpG"), names(tt), value = TRUE)]
      tt[, paste0("c", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\.CpG"), names(tt), value = TRUE)]
      tt[, paste0("p", g) := get(paste0("c", g)) / (get(paste0("t", g))+get(paste0("c",g)))]
    }

    res <- tt[abs(pW - pE) > delta | abs(pY-pW) > delta, {
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
    return( res)
}


pl_venn=function(res){
    library(data.table)
    library(VennDiagram)
    library(grid)

    # Define the four sets based on p-value < 0.01 and log2FC direction
    EW_up  <- res[pvEW < 0.01 & l2EW > 0, paste(chrom, start)]
    EW_down <- res[pvEW < 0.01 & l2EW < 0, paste(chrom, start)]
    WY_up  <- res[pvWY < 0.01 & l2WY > 0, paste(chrom, start)]
    WY_down <- res[pvWY < 0.01 & l2WY < 0, paste(chrom, start)]

    # Combine all four sets in one Venn diagram
    venn.plot <- venn.diagram(
      x = list("EW up" = EW_up, "EW down" = EW_down,
               "WY up" = WY_up, "WY down" = WY_down),
      filename = NULL,           # NULL to plot in R device
      fill = c("red", "pink", "blue", "lightblue"),
      alpha = 0.5,
      cex = 1.2,
      cat.cex = 1.2,
      margin = 0.1,
      main = "Significant Changes: EW vs WY"
    )

    # Draw
    grid.newpage()
    grid.draw(venn.plot)
}

