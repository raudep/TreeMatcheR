# geometry.R

# --- Matrix Utilities (Vectorized) ---

# Create a 4x4 Identity Matrix
trafoMatrixCreateAsUnit <- function() {
  diag(4)
}

# Apply 2D transform to a matrix of points (N x 2 or N x 3)
# Returns N x 2 matrix
homogeneousTrafoMatrix2D <- function(M, pts) {
  # pts: Matrix with columns X, Y
  n <- nrow(pts)
  # Create homogeneous coordinates [x, y, 1]
  # t() is used because matrix mult requires 3xN
  v <- rbind(t(pts[, 1:2]), rep(1, n))

  # Apply M (3x3 part is sufficient for 2D, but we accept 4x4)
  # We only need rows 1, 2, 4 and cols 1, 2, 4 of M for 2D X/Y
  # But assuming standard 4x4:
  res <- M[1:2, 1:2] %*% v[1:2, ] + as.vector(M[1:2, 4])

  # Return as N x 2
  return(t(res))
}

# Estimate 2D Rigid Body Transformation (SVD method)
# X1, Y1: Local points (Source)
# X2, Y2: World points (Target)
homogeneous2DSetTrafoFromVectors <- function(X1, Y1, X2, Y2) {
  # Centroids
  mx1 <- mean(X1); my1 <- mean(Y1)
  mx2 <- mean(X2); my2 <- mean(Y2)

  # Center
  p1 <- rbind(X1 - mx1, Y1 - my1)
  p2 <- rbind(X2 - mx2, Y2 - my2)

  # Covariance
  H <- p1 %*% t(p2)

  # SVD
  s <- svd(H)
  R <- s$v %*% t(s$u)

  # Reflection fix
  if (det(R) < 0) {
    s$v[, 2] <- -s$v[, 2]
    R <- s$v %*% t(s$u)
  }

  # Translation
  t_vec <- c(mx2, my2) - R %*% c(mx1, my1)

  # 4x4 Matrix
  M <- diag(4)
  M[1:2, 1:2] <- R
  M[1:2, 4] <- t_vec

  return(M)
}
