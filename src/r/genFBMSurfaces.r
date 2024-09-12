library(NLMR)


ncol <- 256
nrow <- 256
fractal_dim <- 1.2

for (i in 1:10)
{
    fbm <- NLMR::nlm_fbm(ncol, nrow, fractal_dim)
    raster::writeRaster(fbm, filename = paste0("../nlogo/scratch/fbm_raster", i), format = "ascii")

}