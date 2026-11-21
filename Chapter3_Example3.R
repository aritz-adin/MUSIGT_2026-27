## -------------------------------------------------------------------------------------------------------------
library(dplyr)
library(gstat)
library(mapSpain)
library(sf)
library(tmap)


## -------------------------------------------------------------------------------------------------------------
stations <- read.csv("https://raw.githubusercontent.com/aritz-adin/MUSIGT_2026-27/main/data/manual_20200206.csv",
                     colClasses = c(CODE = "character"))
head(stations)


## -------------------------------------------------------------------------------------------------------------
stations_sf <- stations |>
  st_as_sf(coords=c("X","Y"), crs=25830)

