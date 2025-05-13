
    w <- matrix(c(52, 74, 73, 63, 98, 89, 72, 73, 75), nrow = 3, byrow = TRUE)    
    dz_ew <- ((w[1,1] + 2 * w[2,1] + w[3,1]) - (w[1,3] + 2 * w[2,3] + w[3,3]))  / (8 * 30)
    dz_ns <- ((w[1,1] + 2 * w[1,2] + w[1,3]) - (w[3,1] + 2 * w[3,2] + w[3,3])) / (8 * 30)
    s <- sqrt((dz_ew ^ 2) + (dz_ns ^ 2))
    s <- atan(s) * (180 / pi)

# Horn slope function (average gradient using Sobel-style kernel)
horn_slope_raster <- function(elev, cellsize)
{
  slope_mat <- matrix(NA, nrow = nrow(elev), ncol = ncol(elev))

  for (i in 2:(nrow(elev) - 1)) {
    for (j in 2:(ncol(elev) - 1)) {
      # Get 3x3 window
      z1 <- elev[i-1, j-1]; z2 <- elev[i-1, j]; z3 <- elev[i-1, j+1]
      z4 <- elev[i, j-1];   z5 <- elev[i, j];   z6 <- elev[i, j+1]
      z7 <- elev[i+1, j-1]; z8 <- elev[i+1, j]; z9 <- elev[i+1, j+1]
      
      dzdx <- ((z3 + 2*z6 + z9) - (z1 + 2*z4 + z7)) / (8 * cellsize)
      dzdy <- ((z7 + 2*z8 + z9) - (z1 + 2*z2 + z3)) / (8 * cellsize)
      
      slope <- atan(sqrt(dzdx^2 + dzdy^2))
      slope_mat[i, j] <- slope * (180 / pi)  # Convert to degrees
    }
  }
  
  return(slope_mat)
}

# Steepest descent function
steepest_slope_raster <- function(elev, cellsize) {
  slope_mat <- matrix(NA, nrow = nrow(elev), ncol = ncol(elev))
  
  directions <- matrix(c(
    -1, -1, 1.4142,   -1, 0, 1.0,   -1, 1, 1.4142,
     0, -1, 1.0,                      0, 1, 1.0,
     1, -1, 1.4142,    1, 0, 1.0,    1, 1, 1.4142
  ), byrow = TRUE, ncol = 3)

  for (i in 2:(nrow(elev) - 1)) {
    for (j in 2:(ncol(elev) - 1)) {
      max_slope <- 0
      for (k in 1:nrow(directions)) {
        di <- directions[k, 1]
        dj <- directions[k, 2]
        dist <- directions[k, 3] * cellsize
        
        neighbor_elev <- elev[i + di, j + dj]
        dz <- elev[i, j] - neighbor_elev
        
        slope_deg <- atan(dz / dist) * (180 / pi)
        if (slope_deg > max_slope) {
          max_slope <- slope_deg
        }
      }
      slope_mat[i, j] <- max_slope
    }
  }
  
  return(slope_mat)
}

# Small example DEM
dem <- matrix(c(52, 74, 73, 63, 98, 89, 72, 73, 75), nrow = 3, byrow = TRUE)

cellsize <- 30  # Assuming square 1x1 grid

# Compute both slopes
slope_horn <- horn_slope_raster(dem, cellsize)
slope_steepest <- steepest_slope_raster(dem, cellsize)

# Print rounded results
cat("Horn Slope (degrees):\n")
print(round(slope_horn, 2))

cat("\nSteepest Descent Slope (degrees):\n")
print(round(slope_steepest, 2))
