# QuadTree.R
# Implements a Spatial Grid Index for fast neighborhood search

#' @title Create a Grid Index from a Matrix of coordinates
#' @description Create a Grid Index from a Matrix of coordinates. Implements a Spatial Grid Index for fast neighborhood search.
#'
#' @param pts Matrix (N x 2) with X, Y.
#' @param cellSize Size of grid bins.
#'
#' @return An object of class "GridIndex" containing the index, cellSize, and pts.
#' @export
gridIndexCreate <- function(pts, cellSize) {
  # Calculate integer cell coordinates
  # We use floor() to bucket points
  ix <- floor(pts[,1] / cellSize)
  iy <- floor(pts[,2] / cellSize)

  # Create a unique key string for each cell "ix_iy"
  keys <- paste(ix, iy, sep="_")

  # Split indices by key (creates a hash map equivalent)
  # returns list( "10_5" = c(1, 4, 9), ... )
  index <- split(seq_len(nrow(pts)), keys)

  structure(list(
    index = index,
    cellSize = cellSize,
    pts = pts
  ), class = "GridIndex")
}

#' @title Query the Grid for candidates near multiple query points
#' @description Query the Grid for candidates near multiple query points.
#'
#' @param grid The GridIndex object.
#' @param queryPts Matrix (M x 2).
#' @param radius Search radius.
#'
#' @return List of integer vectors (indices in the original pts matrix).
#' @export
gridIndexQueryBatch <- function(grid, queryPts, radius) {
  # Range of cells to check (radius determines how many neighbor cells)
  delta <- ceiling(radius / grid$cellSize)

  # Query cell coords
  qx <- floor(queryPts[,1] / grid$cellSize)
  qy <- floor(queryPts[,2] / grid$cellSize)

  # For each query point, we need to check (2*delta+1)^2 cells
  # This part is tricky to fully vectorize without exploding memory,
  # but in R, lookups into the 'index' list are fast.

  # We will iterate over the neighbor offsets only (usually 3x3 = 9 iterations)
  # and collect candidates.

  n_query <- nrow(queryPts)
  candidates <- vector("list", n_query)

  # Offsets
  offsets <- expand.grid(dx = -delta:delta, dy = -delta:delta)

  for(i in 1:nrow(offsets)) {
    # Calc neighbor keys
    nx <- qx + offsets$dx[i]
    ny <- qy + offsets$dy[i]
    keys <- paste(nx, ny, sep="_")

    # Batch lookup
    # 'grid$index' is a list. Subsetting it with a vector of keys returns a sublist.
    matches <- grid$index[keys]

    # Matches is a list of length n_query. NULL where no key found.
    # We append these indices to our candidates
    # This loop logic is the only slight bottleneck, but much faster than tree traversal
    for(j in which(!vapply(matches, is.null, logical(1)))) {
      candidates[[j]] <- c(candidates[[j]], matches[[j]])
    }
  }

  return(candidates)
}

#' @title Helper to verify exact distance after grid pruning
#' @description Helper to verify exact distance after grid pruning. Returns indices of World points within radius of Query points.
#'
#' @param queryPts Matrix of query points.
#' @param worldPts Matrix of world points.
#' @param candidatesList List of candidate indices.
#' @param radius Search radius.
#' @param queryRadii Radii for query points.
#' @param worldRadii Radii for world points.
#' @param useWeighted Logical flag to use weighted distances.
#'
#' @return List of vectors (list of lists containing valid indices and distances).
#' @export
filterByDistance <- function(queryPts, worldPts, candidatesList, radius, queryRadii, worldRadii, useWeighted) {

  res <- vector("list", nrow(queryPts))
  r2 <- radius^2

  # This loop iterates over LOCAL points (usually fewer than world).
  # Inside, we perform vectorized distance checks against candidates.
  for(i in seq_along(candidatesList)) {
    idxs <- candidatesList[[i]]
    if(length(idxs) == 0) next

    # Get world coords for these candidates
    w_subset <- worldPts[idxs, , drop=FALSE]
    q_pt <- queryPts[i, ]

    # Vectorized Dist calculation
    dx <- w_subset[,1] - q_pt[1]
    dy <- w_subset[,2] - q_pt[2]
    d2 <- dx*dx + dy*dy

    # Basic filter
    keep <- d2 < r2

    if(any(keep)) {
      # Refine for weighted if necessary
      valid_indices <- idxs[keep]
      valid_dists <- sqrt(d2[keep])

      # Just return raw data for the scorer to handle, or handle scoring here?
      # To keep it generic, we return list(indices, distances)
      res[[i]] <- list(indices = valid_indices, dists = valid_dists)
    }
  }
  return(res)
}
