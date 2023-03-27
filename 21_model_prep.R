library(tidyverse)
library(sf)
if (interactive()) {
  devtools::load_all("~/dev/vecaniso")
} else {
  library(vecaniso)
}

mesh <- read_rds("data/mesh.rds")
goa_arrowtooth <- read_rds("data/goa_arrowtooth.rds")
tri_df <- read_rds("data/tri_df.rds")
vor_df <- read_rds("data/mesh_vor.rds")

data <- va_data(WTCPUE ~ 0 + factor(YEAR),
                data = goa_arrowtooth,
                theta = tri_df$asp90,
                beta = tri_df$slope,
                mesh = mesh)
pars <- va_pars(data, mesh)
map <- va_map(mappars = NULL, pars = pars)
random <- va_random(map)

obj <- va_obj(data, pars, map, random)

sim <- va_sim(obj)

vor_df |>
  st_sf() |>
  mutate(omega = sim$omega) |>
  ggplot() +
  geom_sf(aes(fill = omega), color = NA)
