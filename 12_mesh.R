library(tidyverse)
library(sf)
library(INLA)


locs <- read_rds("data/goa_bts.rds")|>
  select() |>
  distinct() |>
  st_transform(st_crs("+proj=utm +zone=6 ellps=WGS84")) |>
  st_coordinates()

## Create mesh; fine resolution probably important for this case
mesh <- inla.mesh.2d(locs,
                     ## boundary = hull,
                     offset = c(1e3, 250e3),
                     max.edge = c(1e3, 200e3),
                     max.n = 1000,
                     min.angle = c(30, 21),
                     cutoff = 5e3)
                     ## update to `sp` package throws errors here; will deal with
                     ## CRS outside of INLA
                     ## crs = sp::CRS("+proj=utm +zone=6 +datum=WGS84"))
write_rds(mesh, "data/mesh.rds")

## Create `sf` data frame with mesh vertex locations. Useful for plotting etc.
mesh_pt <- map(seq_len(mesh$n), \(i) st_point(mesh$loc[i, 1:2])) |>
  st_sfc(crs = st_crs("+proj=utm +zone=6 +datum=WGS84"))
write_rds(mesh_pt, "data/mesh_pt.rds")

## Construct Voronoi polygons around each mesh vertext. Used for calculating
## integration weights and plotting
mesh_vor <- mesh_pt |>
  st_combine() |>
  st_voronoi() |>
  st_collection_extract() |>
  st_crop(mesh_pt) |>
  st_sfc()
write_rds(mesh_vor, "data/mesh_vor.rds")

  ## vor <- st_voronoi(st_combine(mesh_df))
  ## vor_df <- st_sf(geometry = st_collection_extract(vor), crs = crs) %>%
  ##   ## Join with mesh_df to preserve index. Needs to happen before intersection
  ##   ## with boundary so we don't lose points outside the boundary that still
  ##   ## contribute area within the boundary.
  ##   st_join(mesh_df, join = st_contains)
  ## vor_df2 <- st_intersection(vor_df, bound_df) %>%
  ##   mutate(area = st_area(geometry))

## helper function to construct `sf` polygons from each triangle in the mesh
st_mesh_tri <- function(idx, mesh) {
  tv <- mesh$graph$tv[idx, ]
  tv <- c(tv, tv[1])
  st_polygon(list(mesh$loc[tv, 1:2]))
}
## Construct and save an `sf` data frame with each mesh triangle. Useful later
## for interfacing with bathymetry raster and similar.
mesh_tri <- map(seq_len(nrow(mesh$graph$tv)),
                st_mesh_tri,
                mesh = mesh) |>
  st_sfc(crs = st_crs("+proj=utm +zone=6 ellps=WGS84"))
write_rds(mesh_tri, "data/mesh_tri.rds")
