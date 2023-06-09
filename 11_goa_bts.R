library(tidyverse)
library(sf)

crs <- st_crs(3338)

## Download the survey data into a temporary directory
goa_bts_files <- c("goa1984_1987.zip",
                   "goa1990_1999.zip",
                   "goa2001_2005.zip",
                   "goa2007_2013.zip",
                   "goa2015_2019.zip")
goa_bts_urls <- paste0("https://apps-afsc.fisheries.noaa.gov/RACE/groundfish/survey_data/downloads/",
                       goa_bts_files)

walk(seq_along(goa_bts_files),
     \(idx) {
       if (!file.exists(file.path("rawdata", goa_bts_files[idx])))
         download.file(goa_bts_urls[idx],
                       destfile = file.path("rawdata", goa_bts_files[idx]))})

## Read the survey data
bts_cols <- cols(
  LATITUDE = col_double(),
  LONGITUDE = col_double(),
  STATION = col_character(),
  STRATUM = col_integer(),
  YEAR = col_integer(),
  DATETIME = col_datetime(format = "%m/%d/%Y %H:%M"),
  WTCPUE = col_double(),
  NUMCPUE = col_double(),
  COMMON = col_character(),
  SCIENTIFIC = col_character(),
  SID = col_double(),
  BOT_DEPTH = col_double(),
  BOT_TEMP = col_double(),
  SURF_TEMP = col_double(),
  VESSEL = col_factor(),
  CRUISE = col_factor(),
  HAUL = col_integer()
)

## The second and fourth tables have more than 50k rows. Column names are
## duplicated at the 50,000th row causing parsing warnings. These are safe to
## ignore as these rows are filtered out because they won't have valid
## latitudes or longitudes.
goa_bts <- map_df(
  file.path("rawdata", goa_bts_files),
  read_csv,
  col_types = bts_cols,
  na = c("", "NA", "-9999"),
  skip_empty_rows = TRUE) |>
  filter(!is.na(LATITUDE),
         !is.na(LONGITUDE)) |>
  st_as_sf(coords = c("LONGITUDE",
                      "LATITUDE"),
           crs = st_crs(4326)) |>
  st_transform(crs)
write_rds(goa_bts, "data/goa_bts.rds")

goa_hauls <- goa_bts |>
  select(STATION, STRATUM, YEAR, DATETIME, BOT_DEPTH, BOT_TEMP, SURF_TEMP,
         VESSEL, CRUISE, HAUL) |>
  distinct()
write_rds(goa_hauls, "data/goa_hauls.rds")

## Station 154-96 was sampled by two different vessels in 1990 at the same time,
## so we need to join by vessel as well as station and date-time
goa_arrowtooth <- goa_bts |>
  st_drop_geometry() |>
  filter(SCIENTIFIC == "Atheresthes stomias") |>
  select(STATION, DATETIME, VESSEL, WTCPUE, NUMCPUE, COMMON, SCIENTIFIC) |>
  right_join(goa_hauls,
             by = c("STATION", "DATETIME", "VESSEL")) |>
  mutate(WTCPUE = replace_na(WTCPUE, 0)) |>
  st_sf()
write_rds(goa_arrowtooth, "data/goa_arrowtooth.rds")

## Use the new `st_concave_hull` function to get a tighter polygon around the
## catch locations using a concave hull (requires updated version of GEOS and
## sf)
goa_hauls_hull <- st_concave_hull(st_combine(goa_hauls),
                                  ratio = 0.075,
                                  allow_holes = TRUE)
## Consider subtracting the non-marine areas of the map from the concave hull
ak_state <- rnaturalearth::ne_states("united states of america",
                                     returnclass = "sf") |>
  filter(name == "Alaska") |>
  st_transform(crs)
goa_hauls_hull <- st_difference(goa_hauls_hull, ak_state)
write_rds(goa_hauls_hull, "data/goa_hauls_hull.rds")

### Plot of haul locations and the concave hull.
## goa_hauls_hull |>
##   ## st_buffer(dist = 10e3, nQuadSegs = 1e10) |>
##   ggplot() +
##   geom_sf() +
##   geom_sf(data = goa_hauls)
