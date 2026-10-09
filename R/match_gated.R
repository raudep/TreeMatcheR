# match_gated.R
# Gated variant of stemListMatch(): size acts as a hard gate (not a soft weight),
# rotation is bounded, and final links are found by a global one-to-one assignment
# after iterative rigid refinement. Filter the field table to the trees you expect
# ALS to detect BEFORE calling; every field tree passed in takes part.

#' @title Gated stem list matching
#' @description Aligns a local (field) stem list to a world (ALS) list and links trees
#'   one-to-one. Unlike \code{stemListMatch()}, a pair is only admissible if its size
#'   ratio is <= \code{maxRatio}, so a small field tree can never be linked to a large
#'   ALS tree however close it is. All field trees passed in are used for alignment
#'   and linking, so filter \code{localObj} beforehand (e.g. to dominant trees).
#'
#' @param worldObj,localObj Objects from \code{tableToMatrix()} (cols: x, y, size).
#' @param estX0,estY0 Estimated plot centre in world coordinates.
#' @param searchRes Translation grid step (m) for the coarse search.
#' @param searchDist Half-width (m) of the translation search window.
#' @param linkDist Max horizontal distance (m) for an admissible pair.
#' @param maxRatio Max size ratio max(rL/rW, rW/rL) for an admissible pair (hard gate).
#' @param maxRotDeg Rotation search half-range in degrees around \code{rotCentreDeg};
#'   NULL searches the full 360 degrees as the original function does.
#' @param rotCentreDeg Centre of the rotation search (e.g. magnetic declination).
#' @param nRefine Number of refine iterations (re-link, refit rigid transform).
#' @param sizeW Weight (0-1) of size agreement in the alignment score; 0 = distance only.
#' @param midXLocal,midYLocal Optional local rotation centre (default: bbox middle of
#'   the trees passed in).
#' @param verbose Print progress.
#'
#' @return List: M (4x4), lIndices, wIndices (row indices into the full tables),
#'   theta, score, peakRatio (2nd-best distinct / best score; near 1 = ambiguous),
#'   used (local rows used; rows with NA x, y or size are skipped), pairDist,
#'   pairRatio, flag.
#' @export
stemListMatchGated <- function(worldObj, localObj, estX0, estY0,
                               searchRes = 1.0, searchDist = 5.0, linkDist = 3.0,
                               maxRatio = 1.75, maxRotDeg = 15, rotCentreDeg = 0,
                               nRefine = 3, sizeW = 0.5,
                               midXLocal = NULL, midYLocal = NULL, verbose = TRUE) {
  wMat <- worldObj$mat
  lMat <- localObj$mat
  flag <- "ok"

  # 1. Field trees to use: all rows with valid x, y and size -------------------
  #    (filtering to detectable trees is the caller's job)
  ok <- !is.na(lMat[, 1]) & !is.na(lMat[, 2]) & !is.na(lMat[, 3])
  cand <- which(ok)
  if (any(!ok)) warning(sprintf("%d field tree(s) with NA x, y or size skipped.", sum(!ok)))
  if (length(cand) < 3) {
    flag <- "few_trees_rotation_fixed"
    maxRotDeg <- 0
  }
  lC <- lMat[cand, , drop = FALSE]

  # Rotation centre: bbox middle of the field trees passed in (as in the original)
  if (is.null(midXLocal)) {
    meanL <- (apply(lMat[ok, 1:2, drop = FALSE], 2, min) +
              apply(lMat[ok, 1:2, drop = FALSE], 2, max)) / 2
  } else meanL <- c(midXLocal, midYLocal)
  meanW <- c(estX0, estY0)
  cx <- lC[, 1] - meanL[1]; cy <- lC[, 2] - meanL[2]

  gridW <- gridIndexCreate(wMat[, 1:2, drop = FALSE], max(linkDist, 2.0))

  # Size ratio; unknown sizes are not penalised (ratio 1)
  sizeRatio <- function(a, b) {
    r <- rep(1.0, length(b))
    v <- !is.na(b) & b > 0.001 & !is.na(a) & a > 0.001
    r[v] <- pmax(a / b[v], b[v] / a)
    r
  }

  # 2. Hypothesis grid -------------------------------------------------------
  lRange <- apply(lMat[ok, 1:2, drop = FALSE], 2, function(v) diff(range(v)))
  diagL <- max(sqrt(sum((0.5 * lRange)^2)), 1)
  dTh <- (searchRes / diagL) * 180 / pi
  if (is.null(maxRotDeg)) {
    thetas <- (0:(ceiling(360 / dTh) - 1)) * dTh
  } else {
    K <- ceiling(maxRotDeg / dTh)
    thetas <- rotCentreDeg + (-K:K) * dTh
  }
  steps <- -searchDist + (0:ceiling(2 * searchDist / searchRes)) * searchRes
  hyp <- expand.grid(dx = steps, dy = steps, theta = thetas)
  hyp$score <- 0
  if (verbose) cat(sprintf("Evaluating %d hypotheses on %d field trees...\n",
                           nrow(hyp), length(cand)))

  # Score: one-to-one (per world tree keep the best local), gated by distance AND size,
  # each pair weighted by closeness and size agreement
  for (th in thetas) {
    a <- th * pi / 180
    rx <- cx * cos(a) - cy * sin(a)
    ry <- cx * sin(a) + cy * cos(a)
    idxTh <- which(hyp$theta == th)
    for (h in idxTh) {
      qx <- rx + meanW[1] + hyp$dx[h]
      qy <- ry + meanW[2] + hyp$dy[h]
      cl <- gridIndexQueryBatch(gridW, cbind(qx, qy), linkDist)
      best <- c()
      for (i in which(!vapply(cl, is.null, logical(1)))) {
        w <- cl[[i]]
        d <- sqrt((wMat[w, 1] - qx[i])^2 + (wMat[w, 2] - qy[i])^2)
        g <- d < linkDist & sizeRatio(lC[i, 3], wMat[w, 3]) <= maxRatio
        if (!any(g)) next
        # contribution: closeness x size agreement (1 at equal size, 0 at maxRatio)
        sc <- (1 / (1 + d[g])) * (1 - sizeW * log(sizeRatio(lC[i, 3], wMat[w[g], 3])) / log(maxRatio))
        k <- which.max(sc); wid <- as.character(w[g][k])
        if (is.null(best[wid]) || is.na(best[wid]) || sc[k] > best[wid]) best[wid] <- sc[k]
      }
      if (length(best)) hyp$score[h] <- sum(best) / length(cand)
    }
  }

  b <- which.max(hyp$score)
  bp <- hyp[b, ]
  # Ambiguity: best score among hypotheses clearly different from the winner
  far <- sqrt((hyp$dx - bp$dx)^2 + (hyp$dy - bp$dy)^2) > linkDist |
         abs(hyp$theta - bp$theta) > 3 * dTh
  peakRatio <- if (any(far) && bp$score > 0) max(hyp$score[far]) / bp$score else NA
  if (verbose) cat(sprintf("Coarse best: W=%.4f dx=%.2f dy=%.2f theta=%.2f peakRatio=%.2f\n",
                           bp$score, bp$dx, bp$dy, bp$theta, peakRatio))

  # 3. Initial rigid transform from the winning hypothesis -------------------
  a <- bp$theta * pi / 180
  R <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  M <- diag(4)
  M[1:2, 1:2] <- R
  M[1:2, 4] <- meanW + c(bp$dx, bp$dy) - R %*% meanL

  # 4. Global one-to-one gated assignment + rigid refit (ICP-style) -----------
  assignPairs <- function(M) {
    p <- homogeneousTrafoMatrix2D(M, lC[, 1:2, drop = FALSE])
    cl <- gridIndexQueryBatch(gridW, p, linkDist)
    rows <- list()
    for (i in which(!vapply(cl, is.null, logical(1)))) {
      w <- unique(cl[[i]])
      d <- sqrt((wMat[w, 1] - p[i, 1])^2 + (wMat[w, 2] - p[i, 2])^2)
      r <- sizeRatio(lC[i, 3], wMat[w, 3])
      g <- d < linkDist & r <= maxRatio
      if (any(g)) rows[[length(rows) + 1]] <- data.frame(l = i, w = w[g], d = d[g], r = r[g])
    }
    if (!length(rows)) return(data.frame(l = integer(0), w = integer(0), d = numeric(0), r = numeric(0)))
    E <- do.call(rbind, rows)
    # Additive cost: distance and size mismatch on comparable [0,1] scales
    E$cost <- E$d / linkDist + log(E$r) / log(maxRatio)
    if (requireNamespace("clue", quietly = TRUE) && nrow(E) > 1) {
      ul <- sort(unique(E$l)); uw <- sort(unique(E$w))
      n <- max(length(ul), length(uw))
      C <- matrix(1e6, n, n)
      C[cbind(match(E$l, ul), match(E$w, uw))] <- E$cost
      s <- as.integer(clue::solve_LSAP(C))
      sel <- cbind(seq_len(n), s)
      sel <- sel[sel[, 1] <= length(ul) & sel[, 2] <= length(uw), , drop = FALSE]
      sel <- sel[C[sel] < 1e6, , drop = FALSE]
      key <- paste(ul[sel[, 1]], uw[sel[, 2]])
      E <- E[paste(E$l, E$w) %in% key, ]
    } else {
      # Fallback: globally sorted greedy (cheapest admissible pair first)
      E <- E[order(E$cost), ]
      keep <- !duplicated(E$l)
      E <- E[keep, ]; E <- E[!duplicated(E$w), ]
    }
    E
  }

  E <- assignPairs(M)
  for (it in seq_len(nRefine)) {
    if (nrow(E) < 3) break
    Mnew <- homogeneous2DSetTrafoFromVectors(lC[E$l, 1], lC[E$l, 2], wMat[E$w, 1], wMat[E$w, 2])
    # Keep the refit only if the rotation stays inside the allowed range
    thNew <- atan2(Mnew[2, 1], Mnew[1, 1]) * 180 / pi
    if (!is.null(maxRotDeg) && abs(((thNew - rotCentreDeg + 180) %% 360) - 180) > maxRotDeg + dTh) break
    Enew <- assignPairs(Mnew)
    if (nrow(Enew) < nrow(E)) break
    M <- Mnew; E <- Enew
  }
  if (nrow(E) < 3 && flag == "ok") flag <- "fewer_than_3_links"
  if (!is.na(peakRatio) && peakRatio > 0.9 && flag == "ok") flag <- "ambiguous_alignment"

  list(M = M, lIndices = cand[E$l], wIndices = E$w,
       theta = atan2(M[2, 1], M[1, 1]) * 180 / pi,
       score = bp$score, peakRatio = peakRatio, used = cand,
       pairDist = E$d, pairRatio = E$r, flag = flag)
}

#' @title Gated matcher wrapper
#' @description Drop-in alternative to \code{tableMatchAndGetLinkTableAndTrafoMatLocal()}
#'   using \code{stemListMatchGated()}. Adds columns linkDist, sizeRatio.
#'   Pass only the field trees you want matched; apply the returned M to the full
#'   field table with \code{homogeneousTrafoMatrix2D()} if you need all trees moved.
#' @param tbWorld,tbLocal Data frames (ALS tops, field trees already filtered).
#' @param metaWorld,metaLocal Metadata from \code{TableHeaderLabels()}.
#' @param estX0,estY0 Estimated plot centre in world coordinates.
#' @param ... Passed to \code{stemListMatchGated()}.
#' @return List: tbLink, M, theta, diag (score, peakRatio, flag, nUsed, nLinked).
#' @export
tableMatchGated <- function(tbWorld, tbLocal, metaWorld, metaLocal, estX0, estY0, ...) {
  wObj <- tableToMatrix(tbWorld, metaWorld)
  lObj <- tableToMatrix(tbLocal, metaLocal)
  res <- stemListMatchGated(wObj, lObj, estX0 = estX0, estY0 = estY0, ...)

  tbRes <- tbLocal
  pts <- homogeneousTrafoMatrix2D(res$M, as.matrix(tbRes[, c(metaLocal$lblX, metaLocal$lblY)]))
  tbRes[[paste0(metaLocal$lblX, "_trafo_p")]] <- pts[, 1]
  tbRes[[paste0(metaLocal$lblY, "_trafo_p")]] <- pts[, 2]
  tbRes$isLinked <- 0L
  tbRes$linkDist <- NA_real_; tbRes$sizeRatio <- NA_real_
  for (cc in names(tbWorld)) tbRes[[paste0(cc, "_w")]] <- NA
  if (length(res$lIndices)) {
    tbRes$isLinked[res$lIndices] <- 1L
    tbRes$linkDist[res$lIndices] <- res$pairDist
    tbRes$sizeRatio[res$lIndices] <- res$pairRatio
    for (cc in names(tbWorld)) tbRes[res$lIndices, paste0(cc, "_w")] <- tbWorld[res$wIndices, cc]
  }
  list(tbLink = tbRes, M = res$M, theta = res$theta,
       diag = data.frame(score = res$score, peakRatio = res$peakRatio, flag = res$flag,
                         nUsed = length(res$used), nLinked = length(res$lIndices)))
}
