#' @title Helper class replacement (List)
#' @description Meta Data Definitions. Helper class replacement (List) to define table header labels.
#'
#' @param id Column name for ID, or list of IDs.
#' @param x Column name for X.
#' @param y Column name for Y.
#' @param z Column name for Z, or NULL.
#' @param r Column name for Radius, or NULL.
#' @param t Column name for Time.
#' @param scaleR Scaling factor for radius.
#' @param scan Column name for Scan, or NULL.
#'
#' @return A list containing the mapped metadata column labels and scales.
#' @export
TableHeaderLabels <- function(id, x, y, z, r, t, scaleR, scan) {
  list(lblID_or_list=id, lblX=x, lblY=y, lblZorNone=z, lblRadiusOrNone=r,
       lblTime=t, scaleRadius=scaleR, lblScanOrNone=scan)
}

#' @title Run Plot Match Demo
#' @description Runs a demonstration of the plot matching process using sample data (Air_Positions.txt and Ground_Positions.txt).
#'
#' @return None. Outputs statistics to the console and saves transformed tables to a temporary directory.
#' @export
runDemoPlotMatch <- function() {
  # 1. Load Data safely using system.file
  pathForest <- system.file("extdata", "Air_Positions.txt", package = "TreeMatcheR")
  if (pathForest == "") stop("Sample file Air_Positions.txt not found. Ensure it is in inst/extdata/")
  tbForest <- tableRead(pathForest)
  metaForest <- TableHeaderLabels("ID", "X", "Y", NULL, "H", NULL, 0.15/20.0, NULL)

  pathPlot <- system.file("extdata", "Ground_Positions.txt", package = "TreeMatcheR")
  if (pathPlot == "") stop("Sample file Ground_Positions.txt not found. Ensure it is in inst/extdata/")
  tbPlot <- tableRead(pathPlot)
  metaPlot <- TableHeaderLabels(c("PLOT", "ID"), "X", "Y", NULL, "STEMDIAM", NULL, 0.005, NULL)

  # 2. Match
  searchDist <- 5.0
  estX0 <- 311.99
  estY0 <- 212.44

  # High-performance call
  result <- tableMatchAndGetLinkTableAndTrafoMatLocal(
    tbForest, tbPlot, metaForest, metaPlot,
    searchDist, estX0, estY0
  )

  # 3. Stats & Save to a temporary directory
  stats <- tableCalculateLinkStats(result$tbLink, "X_trafo_p", "Y_trafo_p", "X_w", "Y_w")
  cat(sprintf("Linked Trees: Avg Diff %.4f, SD %.4f, Range %.4f - %.4f\n",
              stats[[1]], stats[[2]], stats[[3]], stats[[4]]))

  out_path_1 <- file.path(tempdir(), "OUT_Ground_Positions_trafo.txt")
  tableSave(result$tbLink, out_path_1)
  cat(sprintf("Saved transformed plot data to temporary directory: %s\n", out_path_1))

  # 4. Transform other file
  pathOther <- system.file("extdata", "Ground_Positions.txt", package = "TreeMatcheR")
  if(file.exists(pathOther) && pathOther != "") {
    tbOther <- tableRead(pathOther)

    # Vectorized Transform
    pts <- as.matrix(tbOther[, c("X", "Y")])
    pts_new <- homogeneousTrafoMatrix2D(result$M, pts)
    tbOther$X <- pts_new[,1]
    tbOther$Y <- pts_new[,2]

    out_path_2 <- file.path(tempdir(), "OUT_GroundPosOther_trafo.txt")
    tableSave(tbOther, out_path_2)
    cat(sprintf("Saved other transformed data to temporary directory: %s\n", out_path_2))
  }
}
