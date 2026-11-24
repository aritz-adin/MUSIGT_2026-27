## -----------------------------------------------------------------------------
## 1) Objectives: 
##    - Predict the maximum temp. (primary variable) in Navarre on 06/02/2010
##      using the minimum temp. (secondary variable) recorded at the same weather stations.
##    - Compare ordinary kriging and ordinary co-kriging using prediction maps,
##      their standard errors and cross-validation.
##    - Learn how to compare both methods fairly.
## -----------------------------------------------------------------------------
  
## -----------------------------------------------------------------------------
## 2) Packages 
## -----------------------------------------------------------------------------
library(dplyr)
library(gstat)
library(mapSpain)
library(sf)
library(tmap)


## -----------------------------------------------------------------------------
## 3) Description of the data 
## -----------------------------------------------------------------------------
stations <- read.csv("https://raw.githubusercontent.com/aritz-adin/MUSIGT_2026-27/main/data/manual_20100206.csv",
                     colClasses = c(CODE = "character"))
head(stations)

stations_sf <- stations |>
  st_as_sf(coords=c("X","Y"), crs=25830)


## Boundary of Navarre and prediction grid -------------------------------------
Navarra <- esp_get_prov(prov = "Navarra") |>
  st_transform(crs = 25830)

grid_sf <- st_make_grid(Navarra, cellsize = 1000, what = "centers") |>
  st_as_sf() |>
  st_filter(Navarra)

nrow(grid_sf)


## Map of weather stations and prediction grid nodes
tmap_mode("plot")

tm_shape(Navarra) + tm_borders(lwd=1.5) +
  tm_shape(grid_sf) + tm_dots(fill="grey", size=0.02) +
  tm_shape(stations_sf) + tm_dots(fill="red", size=0.5) +
  tm_add_legend(type = "symbols",
                labels = c("Grid nodes (1 km)", "Weather stations"),
                fill = c("grey", "red"),
                size = c(0.2, 0.5)) +
  tm_title("Prediction grid") +
  tm_layout(legend.position = c("left", "bottom"))


## Relationship between both variables -----------------------------------------
r <- cor(stations$TMAX, stations$TMIN)

plot(stations$TMIN, stations$TMAX, pch = 19,
     xlab = "TMIN (ºC)", ylab = "TMAX (ºC)",
     main = paste("Correlation:", round(r, 2)))
abline(lm(TMAX ~ TMIN, data = stations), col = "red", lwd = 2)


## -----------------------------------------------------------------------------
## 4) Semivariogram of the primary variable
## -----------------------------------------------------------------------------
v_tmax <- variogram(TMAX ~ 1, stations_sf, width = 10000)
v_tmax

vgm_tmax <- fit.variogram(v_tmax, model=vgm(psill=6, model="Sph",
                                            range=50000, nugget=1))
vgm_tmax

plot(v_tmax, vgm_tmax, pch = 19,
     xlab = "Distance (m)", ylab = "Semivariance",
     main = "Semivariogram of the maximum temperature")


## -----------------------------------------------------------------------------
## 5) Ordinary kriging (TMAX only) and co-kriging (TMAX and TMIN)
## -----------------------------------------------------------------------------
pred_krig <- krige(TMAX ~ 1, stations_sf, grid_sf, model=vgm_tmax, debug.level=0)


## Linear model of coregionalization
g <- gstat(NULL, id = "TMAX", formula = TMAX ~ 1, data = stations_sf)
g <- gstat(g,    id = "TMIN", formula = TMIN ~ 1, data = stations_sf)

v <- variogram(g, width = 10000)
lmc <- fit.lmc(v, g, model=vgm_tmax)
lmc$model

plot(v, lmc, pch = 19, col = "black",
     xlab = "Distance (m)", ylab = "Semivariance",
     scales = list(alternating = 1, tck = c(1, 0)))


## Co-kriging predictions 
pred_cokrig <- predict(lmc, grid_sf, debug.level=0)


## -----------------------------------------------------------------------------
## 6) Comparison of the predictions 
## -----------------------------------------------------------------------------
pred_sf <- grid_sf |>
  mutate(krig_pred   = pred_krig$var1.pred,
         krig_sd     = sqrt(pred_krig$var1.var),
         cokrig_pred = pred_cokrig$TMAX.pred,
         cokrig_sd   = sqrt(pred_cokrig$TMAX.var)) |>
  st_buffer(dist=500, endCapStyle="SQUARE")

tm_shape(pred_sf) +
  tm_polygons(fill = c("krig_pred", "cokrig_pred"),
              fill.scale  = tm_scale_continuous(values = "brewer.yl_or_rd"),
              fill.legend = tm_legend(title = "TMAX (ºC)"),
              fill.free = FALSE,    # common colour scale
              col = NULL) +         # no cell borders
  tm_facets_wrap(nrow=1) +
  tm_shape(Navarra) + tm_borders() +
  tm_shape(stations_sf) + tm_dots(fill="black", size=0.2) +
  tm_layout(panel.labels = c("Ordinary kriging", "Co-kriging"))

tm_shape(pred_sf) +
  tm_polygons(fill = c("krig_sd", "cokrig_sd"),
              fill.scale  = tm_scale_continuous(values = "brewer.blues"),
              fill.legend = tm_legend(title = "SD (ºC)"),
              fill.free = FALSE,
              col = NULL) +
  tm_facets_wrap(nrow=1) +
  tm_shape(Navarra) + tm_borders() +
  tm_shape(stations_sf) + tm_dots(fill="black", size=0.2) +
  tm_layout(panel.labels = c("Ordinary kriging", "Co-kriging"))


## -----------------------------------------------------------------------------
## 7) Cross validation 
## -----------------------------------------------------------------------------
# Ordinary kriging
cv_krig <- krige.cv(TMAX ~ 1, stations_sf, model=vgm_tmax, nfold=nrow(stations_sf),
                    debug.level=0)

# Co-kriging
cv_cokrig <- gstat.cv(lmc, nfold=nrow(stations_sf), remove.all=TRUE, debug.level=0)

cv_results <- data.frame(
  Method = c("Ordinary kriging", "Co-kriging"),
  MAE    = c(mean(abs(cv_krig$residual)), mean(abs(cv_cokrig$residual))),
  RMSE   = c(sqrt(mean(cv_krig$residual^2)), sqrt(mean(cv_cokrig$residual^2)))
)
cv_results |> 
  mutate(MAE = round(MAE, 3),
         RMSE=round(RMSE, 3))
