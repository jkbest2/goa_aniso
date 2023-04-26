library(tidyverse)
library(sf)
## Raster library; includes `terrain` function to get slope, aspect
library(terra)
## Other raster library, interfaces nicely with `sf`
library(stars)
## For downloading bathymetry data from NOAA
library(marmap)
## Shapefiles for coastlines
library(rnaturalearth)
## CMOcean scales for plotting
library(cmocean)

source("10_data_funs.R")

## utm_crs <- st_crs("+proj=utm +zone=6 ellps=WGS84")
crs <- st_crs(3338)

## Buffer distance: used to crop the raster down to the area that matters and to
## enlarge the area of the mesh triangles when extracting slope and aspect
## values, smoothing them somewhat over adjacent triangles.
buffer_dist <- 25e3

## Get the bounding box of the mesh so that we have bathymetry covering all
## triangles in the mesh; use `round_out` to ensure we don't miss anything on
## the edges.
mesh_pt <- read_rds("data/mesh_pt.rds")
mesh_outerbound <- read_rds("data/mesh_outerbound.rds")
## Need to transform to 4326 (lat/long) here for `getNOAA.bathy`
mesh_bbox <- mesh_outerbound |>
  st_buffer(buffer_dist + 10e3) |>
  st_transform(st_crs(4326)) |>
  st_bbox() |>
  round_out(center = c(Inf, Inf, -Inf, -Inf))
## Warping the bathymetry raster below is cropping the southern end of the
## raster too much, losing part of the mesh. Need to go further south to cover
## the mesh after warping. The extra raster will be cropped later (before slope
## and aspect are computed), so this doesn't introduce a lot of extra
## computation.
mesh_bbox[1] <- mesh_bbox[1] - 5
mesh_bbox[2] <- mesh_bbox[2] - 5
mesh_bbox[3] <- mesh_bbox[3] + 5
mesh_bbox[4] <- mesh_bbox[4] + 5

## Download bathymetry from NOAA and save for re-use. NOTE reducing the
## resolution here would be another way to smooth out small-scale variation.
goa_ba <- getNOAA.bathy(lon1 = mesh_bbox[1],
                        lon2 = mesh_bbox[3],
                        lat1 = mesh_bbox[2],
                        lat2 = mesh_bbox[4],
                        resolution = 1,
                        keep = TRUE,
                        antimeridian = FALSE,
                        path = "data")

## Convert to a standard matrix to make it easier to convert to a data frame,
## allowing us to construct a `stars` object
goa_ba2 <- goa_ba
class(goa_ba2) <- NULL

## Pivot to a data frame so we can reconstruct the raster; wasn't able to figure
## out how to cleanly convert directly from a matrix.
goa_df <- as_tibble(goa_ba2) %>%
  mutate(lon = as.numeric(rownames(goa_ba))) %>%
  pivot_longer(cols = !lon,
               names_to = "lat",
               names_transform = as.numeric,
               values_to = "elev")

goa2 <- st_as_stars(goa_df)

## goa_st <- st_as_stars(goa_df) |>
##   st_set_crs(st_crs(4326))

## goa <- st_as_stars(goa_df) %>%
##   st_set_crs(4326) %>%
##   st_transform(crs = st_crs("+proj=utm +zone=6 ellps=WGS84"))

mesh_tri <- read_rds("data/mesh_tri.rds")
mesh_outerbound <- read_rds("data/mesh_outerbound.rds")

## A function to take the mean of directional observations, as suggested in this
## GIS StackExchange answer: https://gis.stackexchange.com/a/147135. Expects
## aspects in radians. `w` argument is optional weights. `w` must be either
## length 1 or the same length as `asp`. `w` is normalized.
mean_aspect <- function(asp, w = 1) {
  if (length(w) != 1 && length(w) != length(asp))
      stop("w must be length 1 or same length as asp")
  ## Normalize the weights
  if (length(w) > 1)
      w <- w / mean(w)
  xs <- w * cos(asp)
  ys <- w * sin(asp)
  atan2(sum(ys), sum(xs))
}

## Using st_transform to go from lat-long to UTM gives a curvilinear raster,
## which the `terrain` function can't handle. So instead resample the raster to
## a regular grid on UTM using st_warp. This results in lots of missing pixels,
## but not in any areas we care about. It also crops out some of the southern
## areas, but a quick plot shows that the raster still covers the mesh
## completely.
goa_bathy <- st_as_stars(goa_df) |>
  st_set_crs(4326) |>
  st_warp(crs = crs) |>
  ## Need to buffer when cropping, or you get missing values in the `terrain`
  ## calculations later on.
  st_crop(st_buffer(mesh_outerbound, buffer_dist + 10e3))
### Plot to double check that mesh is completely covered by the warped raster
## ggplot() +
##   geom_stars(data = goa_bathy) +
##   geom_sf(data = mesh_tri, fill = NA)

write_rds(goa_bathy, "data/goa_bathy.rds")
goa_te <- rast(goa_bathy)

## Convert the mesh triangles to a `terra` vector object so that they can be
## used in `extract` below. Don't need to transform coordinate systems now that
## the raster is `warp`ed
mesh_tri_sv <- mesh_tri |>
  st_buffer(dist = buffer_dist) |>
  vect()

## Calculate the slope and aspect. Use radians so that calculating mean aspect
## (function above) is straightforward. Slope units will be relative anyway,
## though it might be worth considering using grade as a percent rather than
## radians or degrees here, though that is a nonlinear tranformation.
goa_terrain <- terrain(goa_te, v = c("slope", "aspect"), unit = "radians")

tri_df <- extract(goa_terrain, mesh_tri_sv) |>
  group_by(ID) |>
  ## Use percent gradient as slope measure rather than radians
  mutate(slope = tan(slope)) |>
  summarize(slope = mean(slope),
            ## Weight aspect mean calculation by slope. Aspects of higher slopes
            ## should have more effect. Also removes any 90° slopes in areas
            ## with zero slope.
            aspect = mean_aspect(aspect, w = slope),
            .groups = "drop") |>
  mutate(asp90 = aspect + pi / 2,
         ## Make sure all angles are in [-pi, pi]
         aspect = atan2(sin(aspect), cos(aspect)),
         asp90 = atan2(sin(asp90), cos(asp90))) |>
  st_sf(geometry = mesh_tri)
write_rds(tri_df, "data/tri_df.rds")

### Plots for checking slope and aspects
terrplot <- ggplot() +
  geom_sf(aes(fill = aspect, alpha = slope),
          data = tri_df,
          color = NA) +
  scale_fill_cmocean(name = "phase")

## slope_scale <- 5e5
## tri_vecs <- tri_df |>
##   st_centroid() |>
##   mutate(ctr.X = st_coordinates(mesh_tri)[, 1],
##          ctr.Y = st_coordinates(mesh_tri)[, 2],
##          slope = slope * slope_scale,
##          end.X = ctr.X + slope * cos(asp90),
##          end.Y = ctr.Y + slope * sin(asp90)) |>
##   rowwise() |>
##   mutate(lsmat = list(matrix(c(ctr.X, end.X, ctr.Y, end.Y), nrow = 2)),
##          geometry = st_sfc(st_linestring(lsmat))) |>
##   st_sf(sf_column_name = "geometry", crs = crs)

## ggplot() + geom_sf(data = tri_vecs)

## terrplot +
##   geom_sf(data = tri_vecs)

## ggplot() +
##   geom_stars(data = goa_bathy) +
##   geom_sf(data = mesh_tri, fill = NA) +
##   geom_sf(data = tri_vecs) +
##   ggtitle("UTM") +
##   guides(fill = "none")

## tri_df |>
##   mutate(asp_cl = atan2(sin(aspect), cos(aspect)),
##          asp_cut = cut(asp_cl, seq(-pi, pi, length.out = 8)),
##          asp90_cl = atan2(sin(asp90), cos(asp90)),
##          asp90_cut = cut(asp90_cl, seq(-pi, pi, length.out = 6))) |>
##   ggplot() +
##   geom_sf(aes(fill = asp_cut), color = NA, alpha = 0.8) +#, alpha = slope), color = NA) +
##   scale_fill_cmocean(name = "phase", discrete = TRUE) +
##   scale_alpha_continuous() +
##   geom_sf(data = ak_state)

## data <- list(slope = tri_df$slope,
##              theta = tri_df$asp90)

## ell_poly <- function(dfr) {
##   aniso_poly(st_coordinates(dfr$center),
##              dfr$asp90,
##              dfr$slope * 50,
##              rho = 20e3)
## }
## tt <- tri_df |>
##   mutate(center = st_centroid(mesh_tri)) |>
##   st_drop_geometry() |>
##   transmute(center = center,
##             theta = asp90,
##             beta = slope * 50,
##             rho = 20e3)
## tri_ells <- map(seq_len(nrow(tt)),
##                 \(r) aniso_poly(tt$center[r],
##                                 tt$theta[r],
##                                 tt$beta[r],
##                                 tt$rho[r])) |>
##   st_sfc(crs = crs) |>
##   st_sf()

## ak_state <- ne_states("united states of america", returnclass = "sf") |>
##   filter(name == "Alaska")

## tri_bbox <- st_bbox(tri_ells)

## bound_df <- read_rds("data/bound_df.rds")

## tri_ells |>
##   slice_sample(n = 500) |>
## ggplot() +
##   geom_sf(data = ak_state, color = "darkred") +
##   geom_sf(fill = NA) +
##   geom_sf(data = bound_df, fill = NA, color = "darkblue") +
##   coord_sf(xlim = tri_bbox[c(1, 3)],
##            ylim = tri_bbox[c(2, 4)])
