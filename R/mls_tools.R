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

# --- Matcher ---

stemListMatch <- function(worldObj, localObj, searchRes, searchDist, linkDist, cellSize,
                          estX0, estY0, midXLocal=NULL, midYLocal=NULL) {

  wMat <- worldObj$mat
  lMat <- localObj$mat
  nWorld <- nrow(wMat)
  nLocal <- nrow(lMat)

  # FIX 1: Python calculates World Center via Bounding Box (Min+Max)/2, not Mean
  wMin <- apply(wMat[,1:2], 2, min)
  wMax <- apply(wMat[,1:2], 2, max)
  wSize <- wMax - wMin

  # Python Logic: estimatedPlotMiddleXWorld IS PASSED IN as estX0
  # If estX0 is used, Python centers the grid around it.
  meanW <- c(estX0, estY0)

  # FIX 2: Local Center is Bounding Box Middle, NOT Centroid
  if(is.null(midXLocal)) {
    lMin <- apply(lMat[,1:2], 2, min)
    lMax <- apply(lMat[,1:2], 2, max)
    lSize <- lMax - lMin
    meanL <- lMin + 0.5 * lSize
  } else {
    meanL <- c(midXLocal, midYLocal)
  }

  # Center Local Points relative to Bounding Box Middle
  lMatCentered <- lMat
  lMatCentered[,1] <- lMatCentered[,1] - meanL[1]
  lMatCentered[,2] <- lMatCentered[,2] - meanL[2]

  # 3. Build World Spatial Index
  gridW <- gridIndexCreate(wMat[,1:2, drop=FALSE], max(linkDist, 2.0))

  # 4. Generate Search Grid (Exact Python Replication)

  # Angle Steps
  # Python: diagLocal = sqrt( (0.5*DX)^2 + (0.5*DY)^2 )
  lRange <- apply(lMat[,1:2], 2, max) - apply(lMat[,1:2], 2, min)
  diagL <- sqrt( (0.5*lRange[1])^2 + (0.5*lRange[2])^2 )
  if(diagL == 0) diagL <- 1

  deltaThetaDeg <- (searchRes / diagL) * (180 / pi)
  # Python: Ntheta = ceil(360 / delta)
  Ntheta <- ceiling(360.0 / deltaThetaDeg)
  # Python loop: k * deltaThetaDeg for k in 0..Ntheta
  # Note: Python range is 0 to Ntheta-1
  thetas <- (0:(Ntheta-1)) * deltaThetaDeg

  # Grid Steps
  # Python: Nx = ceil( (2.0*searchDist) / res )
  Nx <- ceiling( (2.0 * searchDist) / searchRes )
  Ny <- ceiling( (2.0 * searchDist) / searchRes )

  # Python start offset: startXlocal = estX - searchDist
  # The actual dx added to the centered local is:
  # dx = (estX - searchDist - meanL_X) + i*res
  # But we handle meanL subtraction in lMatCentered.
  # So we just need the relative shift from World Center:
  # shift = -searchDist + i*res

  stepsX <- -searchDist + (0:(Nx-1)) * searchRes
  stepsY <- -searchDist + (0:(Ny-1)) * searchRes

  grid_params <- expand.grid(dx=stepsX, dy=stepsY, theta=thetas)

  cat(sprintf("Evaluating %d hypotheses (Method: Exact Python Replication)...\n", nrow(grid_params)))

  bestWeight <- -1.0
  bestIdx <- 0

  # Pre-compute Rotations
  unique_thetas <- unique(grid_params$theta)
  rads <- unique_thetas * pi / 180
  cos_t <- cos(rads)
  sin_t <- sin(rads)

  for(i in seq_along(unique_thetas)) {
    th <- unique_thetas[i]
    co <- cos_t[i]
    si <- sin_t[i]

    # Rotate ALL local points
    rotX <- lMatCentered[,1] * co - lMatCentered[,2] * si
    rotY <- lMatCentered[,1] * si + lMatCentered[,2] * co

    subset_indices <- which(grid_params$theta == th)
    subset_dx <- grid_params$dx[subset_indices]
    subset_dy <- grid_params$dy[subset_indices]

    for(j in seq_along(subset_indices)) {
      dx <- subset_dx[j]
      dy <- subset_dy[j]

      # Transform to World Est
      currX <- rotX + meanW[1] + dx
      currY <- rotY + meanW[2] + dy

      # Query Neighbors
      candidates <- gridIndexQueryBatch(gridW, cbind(currX, currY), linkDist)

      has_cand <- !vapply(candidates, is.null, logical(1))
      if(!any(has_cand)) next

      valid_indices <- which(has_cand)

      # Python One-to-One Logic
      min_dw_per_world <- rep(Inf, nWorld)

      for(idx in valid_indices) {
        cands <- candidates[[idx]]

        w_sub <- wMat[cands, , drop=FALSE]
        d2 <- (w_sub[,1] - currX[idx])^2 + (w_sub[,2] - currY[idx])^2

        in_range <- d2 < (linkDist^2)
        if(!any(in_range)) next

        cands <- cands[in_range]
        d_vals <- sqrt(d2[in_range])

        rL <- lMat[idx, 3]
        rW <- wMat[cands, 3]

        ratio <- rep(1.0, length(cands))
        valid_r <- (rL > 0.001 & rW > 0.001)
        if(any(valid_r)) ratio[valid_r] <- pmax(rL/rW[valid_r], rW[valid_r]/rL)

        dw <- d_vals * ratio

        best_loc_idx <- which.min(dw)
        best_w_id <- cands[best_loc_idx]
        best_dw <- dw[best_loc_idx]

        if(best_dw < min_dw_per_world[best_w_id]) {
          min_dw_per_world[best_w_id] <- best_dw
        }
      }

      valid_matches <- min_dw_per_world[is.finite(min_dw_per_world)]

      if(length(valid_matches) > 0) {
        score_acc <- sum(1.0 / (1.0 + valid_matches))
        final_score <- score_acc / nLocal

        if(final_score > bestWeight) {
          bestWeight <- final_score
          bestIdx <- subset_indices[j]
        }
      }
    }
  }

  bestP <- grid_params[bestIdx, ]
  cat(sprintf("Best: W=%.4f, dx=%.2f, dy=%.2f, theta=%.2f\n",
              bestWeight, bestP$dx, bestP$dy, bestP$theta))

  # 5. Reconstruct Matches
  rad <- bestP$theta * pi / 180
  rotX <- lMatCentered[,1] * cos(rad) - lMatCentered[,2] * sin(rad)
  rotY <- lMatCentered[,1] * sin(rad) + lMatCentered[,2] * cos(rad)
  finalX <- rotX + meanW[1] + bestP$dx
  finalY <- rotY + meanW[2] + bestP$dy

  candidates <- gridIndexQueryBatch(gridW, cbind(finalX, finalY), linkDist)
  world_best <- vector("list", nWorld)

  for(i in seq_len(nLocal)) {
    matches <- candidates[[i]]
    if(is.null(matches)) next

    w_sub <- wMat[matches, , drop=FALSE]
    d2 <- (w_sub[,1] - finalX[i])^2 + (w_sub[,2] - finalY[i])^2
    valid <- d2 < linkDist^2
    if(!any(valid)) next

    valid_idx <- matches[valid]
    d_vals <- sqrt(d2[valid])

    rL <- lMat[i, 3]; rW <- wMat[valid_idx, 3]
    ratio <- rep(1.0, length(valid_idx))
    valid_r <- (rL > 0.001 & rW > 0.001)
    if(any(valid_r)) ratio[valid_r] <- pmax(rL/rW[valid_r], rW[valid_r]/rL)

    dw <- d_vals * ratio

    best_loc <- which.min(dw)
    w_id <- valid_idx[best_loc]
    val_dw <- dw[best_loc]

    curr <- world_best[[w_id]]
    if(is.null(curr) || val_dw < curr$dw) {
      world_best[[w_id]] <- list(localIdx=i, dw=val_dw)
    }
  }

  final_pairs_idx <- which(!vapply(world_best, is.null, logical(1)))

  M <- diag(4)
  l_indices <- integer(0)
  w_indices <- integer(0)

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
    # Vectorized assignment
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
