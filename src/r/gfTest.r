gaussian_field <- NLMR::nlm_gaussianfield(ncol = 127, 
    nrow = 127, 
    autocorr_range = 5,
    mag_var = 8, 
    nug = 5)
    
gf <- terra::rast(gaussian_field)
gf_df <- terra::as.data.frame(gf, xy = TRUE)
y

library(tidyverse)
library(tidyterra)

ggplot() +
    geom_spatraster(data = terra::rast(fbm_raster), aes(fill = layer))
#####
fbm_raster  <- NLMR::nlm_fbm(ncol = 5, nrow = 5, fract_dim = 1.5)
fbm_df <- as.data.frame(fbm_raster, xy = TRUE)
fbm_df[,1] <- floor(fbm_df[,1])
fbm_df[,2] <- floor(fbm_df[,2])
fbm_df[1:2,]


r <- raster(ncol=3, nrow=3)
values(r) <- sqrt(1:ncell(r))
r[3:5] <- NA
as.data.frame(fbm_raster, xy = TRUE)
s <- stack(r, r*2)
as.data.frame(s)
as.data.frame(s, na.rm=TRUE)
