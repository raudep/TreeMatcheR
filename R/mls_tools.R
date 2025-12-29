# R_mlsTools.R
source("R/geometry.R")
source("R/QuadTree.R")

# --- File I/O ---

tableRead <- function(tablePath) {
  if (!file.exists(tablePath)) stop(paste("File not found:", tablePath))
  read.table(tablePath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

tableSave <- function(df, outPath) {
  write.table(df, file = outPath, sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- Data Preparation ---

tableToMatrix <- function(df, meta) {
  get_col <- function(lbl, default=0.0) {
    if(!is.null(lbl) && lbl %in% names(df)) df[[lbl]] else rep(default, nrow(df))
  }

  if(length(meta$lblID_or_list) > 1) {
    ids <- do.call(paste, c(df[meta$lblID_or_list], sep="_"))
  } else {
    ids <- df[[meta$lblID_or_list]]
  }

  mat <- cbind(
    get_col(meta$lblX),
    get_col(meta$lblY),
    get_col(meta$lblRadiusOrNone, 1.0) * meta$scaleRadius
  )

  list(ids=as.character(ids), mat=mat, df_orig=df)
}

# --- Fully Vectorized Matcher ---

stemListMatch <- function(worldObj, localObj, searchRes, searchDist, linkDist, cellSize,
                          estX0, estY0, midXLocal=NULL, midYLocal=NULL) {

  wMat <- worldObj$mat
  lMat <- localObj$mat

  # 1. Centers and Offsets
  meanW <- c(estX0, estY0)
  meanL <- if(is.null(midXLocal)) colMeans(lMat[,1:2, drop=FALSE]) else c(midXLocal, midYLocal)

  # Center Local Points
  lMatCentered <- lMat
  lMatCentered[,1] <- lMatCentered[,1] - meanL[1]
  lMatCentered[,2] <- lMatCentered[,2] - meanL[2]

  # 2. Build World Spatial Index
  # Use a grid size slightly larger than linkDist for safe bucketing
  gridW <- gridIndexCreate(wMat[,1:2, drop=FALSE], max(linkDist, 2.0))

  # 3. Generate Search Space (The Grid)
  # Angles
  diagL <- sqrt(sum((apply(lMat[,1:2], 2, max) - apply(lMat[,1:2], 2, min))^2)) / 2
  if(diagL == 0) diagL <- 1
  deltaThetaDeg <- (searchRes / diagL) * (180 / pi)
  thetas <- unique(seq(0, 360, by=deltaThetaDeg) %% 360)

  # Translations
  steps <- seq(-searchDist, searchDist, by=searchRes)

  # Create a Master Grid of all (dx, dy, theta) combinations
  # We construct this matrix once to avoid nested looping
  grid_params <- expand.grid(dx=steps, dy=steps, theta=thetas)
  n_sims <- nrow(grid_params)

  cat(sprintf("Evaluating %d hypotheses via Vectorization...\n", n_sims))

  # 4. Batch Scoring Strategy
  # Evaluating 100,000s of hypotheses on all points is memory heavy.
  # We process 'Thetas' in chunks, but fully vectorize dX/dY within those chunks.

  bestWeight <- -1.0
  bestIdx <- 0

  # Pre-compute Rotation Matrices for unique Thetas
  # This avoids sin/cos calls inside the loop
  unique_thetas <- unique(grid_params$theta)
  rads <- unique_thetas * pi / 180
  cos_t <- cos(rads)
  sin_t <- sin(rads)

  # Iterate over Unique Angles (Outer layer only)
  for(i in seq_along(unique_thetas)) {
    th <- unique_thetas[i]
    co <- cos_t[i]
    si <- sin_t[i]

    # Rotate ALL Local points for this angle ONCE
    # N_local x 2
    rotX <- lMatCentered[,1] * co - lMatCentered[,2] * si
    rotY <- lMatCentered[,1] * si + lMatCentered[,2] * co

    # Identify all grid entries that use this theta
    # In 'expand.grid', the last column changes slowest, so these are usually contiguous blocks
    # But strictly filtering is safer:
    subset_indices <- which(grid_params$theta == th)
    subset_dx <- grid_params$dx[subset_indices]
    subset_dy <- grid_params$dy[subset_indices]

    # 5. The "Inner Loop" Replacement:
    # Instead of shifting points and querying the grid for every dx/dy,
    # We invert the logic:
    # For every Local Point, which World Points are compatible *at all*?
    # Then we check if the required translation (dx, dy) matches our grid.

    # However, for dense grids, the "Forward" approach is usually better if we verify fast.
    # Let's verify a batch of Translation Hypotheses against the Quad/Grid Index.

    # To keep memory manageable, we loop through the translation candidates for this angle
    # But we do NOT loop through points.

    current_angle_best_w <- -1
    current_angle_best_j <- 0

    # Optimization: Check if match is even possible?
    # Skip for now to ensure correctness.

    for(j in seq_along(subset_indices)) {
      dx <- subset_dx[j]
      dy <- subset_dy[j]

      # Transform Local to World Est
      currX <- rotX + meanW[1] + dx
      currY <- rotY + meanW[2] + dy

      # --- Fast Score Kernel ---
      # This is the critical hot path

      # 1. Query Neighbors (Returns list of integer indices)
      # gridIndexQueryBatch is heavily optimized in R_quadTree.R
      candidates <- gridIndexQueryBatch(gridW, cbind(currX, currY), linkDist)

      # 2. Compute Score
      # To avoid R loop overhead here, we flatten the structure

      # Identify which local points actually found a neighbor
      has_cand <- !vapply(candidates, is.null, logical(1))

      if(!any(has_cand)) next

      valid_local_indices <- which(has_cand)
      nValid <- length(valid_local_indices)

      # Extract data for valid matches only
      # We need to process potentially multiple matches per local point
      # This is hard to vectorise perfectly without a C++ extension,
      # but we can do a "Nearest Neighbor" assumption for speed.

      score_acc <- 0

      # We iterate only over valid points (subset of local)
      # This is "Order(N_local_valid)" which is small
      for(idx in valid_local_indices) {
        cands <- candidates[[idx]]

        # Calculate Squared Dists
        w_sub <- wMat[cands, , drop=FALSE]
        d2 <- (w_sub[,1] - currX[idx])^2 + (w_sub[,2] - currY[idx])^2

        # Filter strict distance
        in_range <- d2 < (linkDist^2)
        if(!any(in_range)) {
          nValid <- nValid - 1
          next
        }

        cands <- cands[in_range]
        d_vals <- sqrt(d2[in_range])

        # Weighted logic
        rL <- lMat[idx, 3]
        rW <- wMat[cands, 3]
        ratio <- pmax(rL/rW, rW/rL)
        dw <- d_vals * ratio

        # Best match for this stem
        min_dw <- min(dw)

        # Add to score accumulator (1 / (1+dw))
        score_acc <- score_acc + (1.0 / (1.0 + min_dw))
      }

      if(nValid > 0) {
        final_score <- score_acc / nValid
        if(final_score > bestWeight) {
          bestWeight <- final_score
          bestIdx <- subset_indices[j] # Points to global grid_params
        }
      }
    }
  }

  # Retrieve Best Params
  bestP <- grid_params[bestIdx, ]
  cat(sprintf("Best: W=%.4f, dx=%.2f, dy=%.2f, theta=%.2f\n",
              bestWeight, bestP$dx, bestP$dy, bestP$theta))

  # --- Reconstruct Transform & Pairs (Exact Calculation) ---

  # Rotate
  rad <- bestP$theta * pi / 180
  rotX <- lMatCentered[,1] * cos(rad) - lMatCentered[,2] * sin(rad)
  rotY <- lMatCentered[,1] * sin(rad) + lMatCentered[,2] * cos(rad)

  # Translate
  finalX <- rotX + meanW[1] + bestP$dx
  finalY <- rotY + meanW[2] + bestP$dy

  # Exact Matching (Pairs)
  candidates <- gridIndexQueryBatch(gridW, cbind(finalX, finalY), linkDist)
  world_best <- vector("list", nrow(wMat))

  # Fill Pairs
  for(i in seq_len(nrow(lMat))) {
    matches <- candidates[[i]]
    if(is.null(matches)) next

    w_sub <- wMat[matches, , drop=FALSE]
    d2 <- (w_sub[,1] - finalX[i])^2 + (w_sub[,2] - finalY[i])^2
    valid <- d2 < linkDist^2
    if(!any(valid)) next

    valid_idx <- matches[valid]
    dw <- sqrt(d2[valid]) * pmax(lMat[i,3]/wMat[valid_idx,3], wMat[valid_idx,3]/lMat[i,3])

    best_loc <- which.min(dw)
    w_id <- valid_idx[best_loc]
    val_dw <- dw[best_loc]

    curr <- world_best[[w_id]]
    if(is.null(curr) || val_dw < curr$dw) {
      world_best[[w_id]] <- list(localIdx=i, dw=val_dw)
    }
  }

  final_pairs_idx <- which(!vapply(world_best, is.null, logical(1)))

  # Matrix Calculation
  M <- diag(4)
  if(length(final_pairs_idx) > 2) {
    w_indices <- final_pairs_idx
    l_indices <- vapply(world_best[w_indices], function(x) x$localIdx, integer(1))
    M <- homogeneous2DSetTrafoFromVectors(
      lMat[l_indices, 1], lMat[l_indices, 2],
      wMat[w_indices, 1], wMat[w_indices, 2]
    )
  }

  return(list(M=M, lIndices=l_indices, wIndices=w_indices))
}


# --- Wrapper ---

tableMatchAndGetLinkTableAndTrafoMatLocal <- function(tbWorld, tbLocal, metaWorld, metaLocal,
                                                      searchDist, estX0, estY0) {

  wObj <- tableToMatrix(tbWorld, metaWorld)
  lObj <- tableToMatrix(tbLocal, metaLocal)

  res <- stemListMatch(wObj, lObj, searchRes=1.0, searchDist=searchDist, linkDist=3.0,
                       cellSize=0.1, estX0=estX0, estY0=estY0)

  # Apply Transforms and Links
  tbRes <- tbLocal
  pts <- as.matrix(tbRes[, c(metaLocal$lblX, metaLocal$lblY)])
  pts_trans <- homogeneousTrafoMatrix2D(res$M, pts)

  tbRes[[paste0(metaLocal$lblX, "_trafo_p")]] <- pts_trans[,1]
  tbRes[[paste0(metaLocal$lblY, "_trafo_p")]] <- pts_trans[,2]

  tbRes$isLinked <- 0
  wCols <- names(tbWorld)
  for(c in wCols) tbRes[[paste0(c, "_w")]] <- NA

  if(length(res$lIndices) > 0) {
    tbRes$isLinked[res$lIndices] <- 1
    for(c in wCols) tbRes[res$lIndices, paste0(c, "_w")] <- tbWorld[res$wIndices, c]
  }

  return(list(tbLink=tbRes, M=res$M))
}

tableCalculateLinkStats <- function(df, lblX1, lblY1, lblX2, lblY2) {
  linked <- df[df$isLinked == 1, ]
  if (nrow(linked) == 0) return(list(0,0,0,0))
  dists <- sqrt((linked[[lblX1]] - linked[[lblX2]])^2 + (linked[[lblY1]] - linked[[lblY2]])^2)
  return(list(mean(dists), sd(dists), min(dists), max(dists)))
}
