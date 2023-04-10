library(tidyverse)
library(sf)
## Raster library; includes `terrain` function to get slope, aspect
library(terra)
## Other raster library, interfaces nicely with `sf`
library(stars)
library(starsExtra)
## For downloading bathymetry data from NOAA
library(marmap)
## Shapefiles for coastlines
library(rnaturalearth)
## CMOcean scales for plotting
library(cmocean)

source("10_data_funs.R")

utm_crs <- st_crs("+proj=utm +zone=6 ellps=WGS84")

## Get the bounding box of the mesh so that we have bathymetry covering all
## triangles in the mesh; use `round_out` to ensure we don't miss anything on
## the edges.
mesh_pt <- read_rds("data/mesh_pt.rds")
mesh_bbox <- st_bbox(st_transform(mesh_pt, st_crs(4326))) |>
  round_out(center = c(Inf, Inf, -Inf, -Inf))

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

goa_st <- st_as_stars(goa_df) |>
  st_set_crs(st_crs(4326))

## goa <- st_as_stars(goa_df) %>%
##   st_set_crs(4326) %>%
##   st_transform(crs = st_crs("+proj=utm +zone=6 ellps=WGS84"))

## goa_te <- rast(goa_df, crs = "+proj=lonlat")
## goa_te_slope <- terrain(goa_te, "slope", unit = "radians", neighbors = 8)
## goa_te_aspect <- terrain(goa_te, "aspect", unit = "radians", neighbors = 8)

## goa_st <- st_as_stars(goa_te) %>%
##   st_transform(utm_crs)
## goa_st_slope <- st_as_stars(goa_te_slope) %>%
##   st_transform(utm_crs)
## goa_st_aspect <- st_as_stars(goa_te_aspect) %>%
##   st_transform(utm_crs)

## ## Get a polygon of Alaska to delineate coastal areas etc
## ak_coast <- ne_states("united states of america", returnclass = "sf") %>%
##   filter(name == "Alaska") %>%
##   st_transform(utm_crs) %>%
##   st_crop(goa_utm) %>%
##   sf:::select.sf(name)
##

## goa_utm <- st_transform(goa_st, utm_crs)

## st_extract(goa_utm, at = mesh_pt[2345])

mesh <- read_rds("data/mesh.rds")

## Transform to lat/long so that the mesh vertices match the projection of the
## bathymetry.
mesh_ll <- st_transform(mesh_pt, st_crs(4326))

mesh_elev <- st_extract(goa_st, mesh_ll) |>
  st_transform(utm_crs)

mesh_coords <- cbind(st_coordinates(mesh_elev), mesh_elev$elev)

## A function that calculates the vector normal to the triangle formed by the
## coordinates, including elevation.
tri_normvecs <- map(seq_len(nrow(mesh$graph$tv)),
    \(i) {
      tris <- mesh$graph$tv[i, ]
      norm_vec(mesh_coords[tris[1], ],
               mesh_coords[tris[2], ],
               mesh_coords[tris[3], ])
    })


tri_dz <- map(tri_normvecs, dz_from_normvec)
tri_dznorm <- map_dbl(tri_dz, norm, type = "2")

tri_asp <- map_dbl(tri_dz, aspect_from_dz)

mesh_tri <- read_rds("data/mesh_tri.rds")
tri_df <- tibble(geometry = mesh_tri,
       slope = tri_dznorm,
       aspect = tri_asp,
       asp_deg = 180 * aspect / pi,
       asp90 = aspect + pi / 2) |>
  st_sf()
write_rds(tri_df, "data/tri_df.rds")

tri_df |>
  ggplot() +
  geom_sf(aes(fill = asp90, alpha = slope), color = NA) +
  scale_fill_cmocean(name = "phase") +
  scale_alpha_continuous()

data <- list(slope = tri_df$slope,
             theta = tri_df$asp90)

## ell_poly <- function(dfr) {
##   aniso_poly(st_coordinates(dfr$center),
##              dfr$asp90,
##              dfr$slope * 50,
##              rho = 20e3)
## }
tt <- tri_df |>
  mutate(center = st_centroid(geometry)) |>
  st_drop_geometry() |>
  transmute(center = center,
            theta = asp90,
            beta = slope * 50,
            rho = 20e3)
tri_ells <- map(seq_len(nrow(tt)),
                \(r) aniso_poly(tt$center[r],
                                tt$theta[r],
                                tt$beta[r],
                                tt$rho[r])) |>
  st_sfc(crs = utm_crs) |>
  st_sf()

ak_state <- ne_states("united states of america", returnclass = "sf") |>
  filter(name == "Alaska")

tri_bbox <- st_bbox(tri_ells)

tri_ells |>
  ## slice_sample(n = 250) |>
ggplot() +
  geom_sf(data = ak_state, color = "darkred") +
  geom_sf(fill = NA) +
  coord_sf(xlim = tri_bbox[c(1, 3)],
           ylim = tri_bbox[c(2, 4)])
