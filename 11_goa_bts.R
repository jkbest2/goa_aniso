library(tidyverse)
library(sf)

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
           crs = st_crs(4326))
write_rds(goa_bts, "data/goa_bts.rds")

goa_hauls <- goa_bts |>
  select(STATION, STRATUM, YEAR, DATETIME, BOT_DEPTH, BOT_TEMP, SURF_TEMP,
         VESSEL, CRUISE, HAUL) |>
  distinct()
write_rds(goa_hauls, "data/goa_hauls.rds")

goa_arrowtooth <- goa_bts |>
  st_drop_geometry() |>
  filter(SCIENTIFIC == "Atheresthes stomias") |>
  select(STATION, DATETIME, WTCPUE, NUMCPUE, COMMON, SCIENTIFIC) |>
  right_join(goa_hauls,
             by = c("STATION", "DATETIME")) |>
  mutate(WTCPUE = replace_na(WTCPUE, 0)) |>
  st_sf()
write_rds("goa_arrowtooth", "data/goa_arrowtooth.rds")
