library(tidyverse)
library(sf)
library(INLA)

## crs <- st_crs("+proj=utm +zone=6 +datum=WGS84")
crs <- st_crs(3338)

locs <- read_rds("data/goa_bts.rds")|>
  select() |>
  distinct() |>
  st_transform(crs) |>
  st_coordinates()

hull <- read_rds("data/goa_hauls_hull.rds") |>
  st_transform(crs)

## Create mesh; resolution will probably be important as it also determines how
## fine the gradients for the anisotropy field are calculated, though it is
## possible to decouple these later on. The primary parameter that determines
## the number of mesh points is the `cutoff`. At 25e3 we get roughly 700
## vertices, at 20e3 we get 950, increasing as the cutoff decreases.
mesh <- inla.mesh.2d(locs,
                     boundary = as(hull, "Spatial"),
                     offset = c(1e3, 250e3),
                     max.edge = c(10e3, 200e3),
                     max.n = 1000,
                     min.angle = c(30, 21),
                     cutoff = 20e3)
                     ## update to `sp` package throws errors here; will deal with
                     ## CRS outside of INLA
                     ## crs = sp::CRS("+proj=utm +zone=6 +datum=WGS84"))
write_rds(mesh, "data/mesh.rds")

## Create `sf` data frame with mesh vertex locations. Useful for plotting etc.
mesh_pt <- map(seq_len(mesh$n), \(i) st_point(mesh$loc[i, 1:2])) |>
  ## st_sfc(crs = st_crs("+proj=utm +zone=6 +datum=WGS84"))
  st_sfc(crs = crs)
write_rds(mesh_pt, "data/mesh_pt.rds")

## Use the concave hull around the bottom trawl survey with a 10km
## buffer for good measure.
bound_df <- st_sf(hull, crs = crs) |>
  st_buffer(dist = 10e3)
write_rds(bound_df, "data/bound_df.rds")

## Create an `sf` dataframe with mesh vertices as features and an index column
## so that voronoi polygons can be associated back to corresponding vertices.
mesh_df <- mesh_pt |>
  st_sf() |>
  rowid_to_column("vertex_id")

## Construct Voronoi polygons around each mesh vertext. Used for calculating
## integration weights and plotting
mesh_vor <- mesh_pt |>
  st_combine() |>
  st_voronoi() |>
  st_collection_extract() |>
  st_crop(mesh_pt) |>
  st_sfc()
write_rds(mesh_vor, "data/mesh_vor.rds")

vor <- st_voronoi(st_combine(mesh_pt))
vor_df <- st_sf(geometry = st_collection_extract(vor), crs = crs) %>%
  ## Join with mesh_df to preserve index. Needs to happen before intersection
  ## with boundary so we don't lose points outside the boundary that still
  ## contribute area within the boundary.
  st_join(st_sf(mesh_pt), join = st_contains)
vor_df2 <- st_intersection(vor_df, bound_df) %>%
  mutate(area = st_area(geometry))

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
  ## st_sfc(crs = st_crs("+proj=utm +zone=6 ellps=WGS84"))
  st_sfc(crs = crs)
write_rds(mesh_tri, "data/mesh_tri.rds")

mesh_outerbound <- mesh_tri |>
  st_union()
write_rds(mesh_outerbound, "data/mesh_outerbound.rds")
