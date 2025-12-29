# demoPlotMatchTables.R
source("R/mls_tools.R")

# --- Meta Data Definitions ---
# Helper class replacement (List)
TableHeaderLabels <- function(id, x, y, z, r, t, scaleR, scan) {
  list(lblID_or_list=id, lblX=x, lblY=y, lblZorNone=z, lblRadiusOrNone=r,
       lblTime=t, scaleRadius=scaleR, lblScanOrNone=scan)
}

# 1. Load Data
tbForest <- tableRead("stemListForest.txt")
metaForest <- TableHeaderLabels("ID", "X", "Y", NULL, "H", NULL, 0.15/20.0, NULL)

tbPlot <- tableRead("stemListPlot.txt")
metaPlot <- TableHeaderLabels(c("PLOT", "ID"), "X", "Y", NULL, "STEMDIAM", NULL, 0.005, NULL)

# 2. Match
searchDist <- 5.0
estX0 <- 1545813.11
estY0 <- 6935928.56

# High-performance call
result <- tableMatchAndGetLinkTableAndTrafoMatLocal(
  tbForest, tbPlot, metaForest, metaPlot,
  searchDist, estX0, estY0
)

# 3. Stats & Save
stats <- tableCalculateLinkStats(result$tbLink, "X_trafo_p", "Y_trafo_p", "X_w", "Y_w")
cat(sprintf("Linked Trees: Avg Diff %.4f, SD %.4f, Range %.4f - %.4f\n",
            stats[[1]], stats[[2]], stats[[3]], stats[[4]]))

tableSave(result$tbLink, "OUT_stemListPlot_trafo.txt")

# 4. Transform other file
if(file.exists("stemListOther.txt")) {
  tbOther <- tableRead("stemListOther.txt")

  # Vectorized Transform
  pts <- as.matrix(tbOther[, c("X", "Y")])
  pts_new <- homogeneousTrafoMatrix2D(result$M, pts)
  tbOther$X <- pts_new[,1]
  tbOther$Y <- pts_new[,2]

  tableSave(tbOther, "OUT_stemListOther_trafo.txt")
}
