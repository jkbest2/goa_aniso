## Round to nearest integer *away* from a central value
round_out <- function(x, center = 0) {
  ifelse(x < center, floor(x), ceiling(x))
}

## Calculate the cross product of two vectors to get the vector normal to the
## plane containing those two vectors.
vec_cross <- function(u1, u2) {
  i1 <- c(2, 3, 1)
  i2 <- c(3, 1, 2)

  u1[i1] * u2[i2] - u1[i2] * u2[i1]
}

## Get vector normal to the plane that contains points `v1`, `v2`, and `v3`
norm_vec <- function(v1, v2, v3) {
  u1 <- v1 - v2
  u2 <- v1 - v3
  vec_cross(u1, u2)
}

## Get slope vector from normal vector
dz_from_normvec <- function(normvec) {
  c(normvec[1] / normvec[3], normvec[2] / normvec[3])
}

## Get the aspect of the slope in radians, with north at zero.
aspect_from_dz <- function(dz) {
  if (dz[2] == 0) {
    az <- 0
  } else {
    az <- atan2(dz[1], dz[2])
  }
  az
}

## Calculate H anisotropy matrix
det1_gamma <- function(beta) {
  1/2 * (sqrt(beta^2 + 4) - beta)
}

# Construct anisotropy matrix
aniso_h <- function(theta, beta, gamma = NULL, det1 = TRUE) {
  if (is.null(gamma) && det1) {
    gamma <- det1_gamma(beta)
  }
  if (is.null(gamma)) {
    stop("Need a value for gamma")
  }
  beta >= 0 || stop("Need beta > 0")
  gamma > 0 || stop("Need gamma > 0")
  v <- c(cos(theta), sin(theta))
  diag(c(gamma, gamma)) + beta * v %*% t(v)
}

aniso_poly <- function(center, theta, beta, rho) {
  if (inherits(center, "sfc_POINT")) center <- c(st_coordinates(center))
  H <- aniso_h(theta, beta)
  a <- t(seq(0, 2 * pi, length.out = 257))
  a[257] <- a[1]
  xy <- rbind(cos(a), sin(a))
  ell <- rho * H %*% xy
  st_polygon(list(t(ell + center)))
}
