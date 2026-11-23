## -----------------------------------------------------------------------------
## 1) Objectives:
##    - Predict the concentration of lead (primary variable) in the Meuse
##      study area using zinc (secondary variable), in the typical situation
##      where the secondary variable has been observed at many more locations
##      than the primary variable.
##    - Compare ordinary kriging (lead only) and ordinary co-kriging (lead and
##      zinc) using prediction maps, their standard errors and cross-validation.
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## 2) Packages
## -----------------------------------------------------------------------------
library(dplyr)
library(gstat)
library(sf)
library(tmap)


## -----------------------------------------------------------------------------
## 3) Description of the data
## -----------------------------------------------------------------------------
data(meuse, package="sp")         # soil samples
data(meuse.grid, package="sp")    # prediction grid

meuse_sf <- meuse |>
  as_tibble() |>
  st_as_sf(coords=c("x","y"), crs=28992)

grid_sf <- meuse.grid |>
  as_tibble() |>
  st_as_sf(coords=c("x","y"), crs=28992)


## Select a random subset of 50 samples ----------------------------------------
set.seed(123)
idx <- sample(nrow(meuse_sf), 50)

lead_sf <- meuse_sf |> slice(idx)   # primary variable: 50 samples
zinc_sf <- meuse_sf                 # secondary variable: 155 samples


## -----------------------------------------------------------------------------
## 4) Ordinary kriging (lead only)
## -----------------------------------------------------------------------------
v_lead <- variogram(log(lead) ~ 1, lead_sf)
vgm_lead <- fit.variogram(v_lead, model=vgm(psill=0.5, model="Sph",
                                            range=1000, nugget=0.05))
vgm_lead

pred_krig <- krige(log(lead) ~ 1, lead_sf, grid_sf, model=vgm_lead,
                   debug.level=0)
head(pred_krig)


## -----------------------------------------------------------------------------
## 5) Ordinary co-kriging (lead and zinc)
## -----------------------------------------------------------------------------
## Direct and cross empirical semivariograms -----------------------------------
g <- gstat(NULL, id = "Pb", formula = log(lead) ~ 1, data = lead_sf)
g <- gstat(g,    id = "Zn", formula = log(zinc) ~ 1, data = zinc_sf)

v <- variogram(g)
v |> filter(dist==0)

v <- v |> filter(dist>0)


## Linear model of coregionalization -------------------------------------------
lmc <- fit.lmc(v, g, model=vgm_lead)
lmc$model


## Plot of the fitted model ----------------------------------------------------
plot(v, lmc, pch = 19, col = "black",
     xlab = "Distance (m)", ylab = "Semivariance",
     scales = list(alternating = 1, tck = c(1, 0)))


## Co-kriging predictions ------------------------------------------------------
pred_cokrig <- predict(lmc, grid_sf, debug.level=0)
head(pred_cokrig)


## -----------------------------------------------------------------------------
## 6) Comparison of the predictions
## -----------------------------------------------------------------------------
## Interactive maps of predictions and standard errors -------------------------
pred_sf <- grid_sf |>
  mutate(krig_pred   = pred_krig$var1.pred,
         krig_sd     = sqrt(pred_krig$var1.var),
         cokrig_pred = pred_cokrig$Pb.pred,
         cokrig_sd   = sqrt(pred_cokrig$Pb.var)) |>
  st_buffer(dist=20, endCapStyle="SQUARE")

Map.pred <- tm_shape(pred_sf) +
  tm_polygons(fill = c("krig_pred", "cokrig_pred"),
              fill.scale  = tm_scale_continuous(values = "brewer.yl_or_rd"),
              fill.legend = tm_legend(title = "log(lead)"),
              fill.free = FALSE,    # common colour scale
              col = NULL,           # no cell borders
              popup.vars = c("Kriging"    = "krig_pred",
                             "Co-kriging" = "cokrig_pred")) +
  tm_facets_wrap(nrow=1, sync=TRUE) +
  tm_shape(lead_sf) + tm_dots(fill="black", size=0.4) +
  tm_basemap("OpenStreetMap") +
  tm_layout(panel.labels = c("Ordinary kriging", "Co-kriging")) + 
  tm_view(use_WebGL = FALSE)

Map.sd <- tm_shape(pred_sf) +
  tm_polygons(fill = c("krig_sd", "cokrig_sd"),
              fill.scale  = tm_scale_continuous(values = "brewer.blues"),
              fill.legend = tm_legend(title = "SD (log scale)"),
              fill.free = FALSE,    # common colour scale
              col = NULL,           # no cell borders
              popup.vars = c("Kriging SE"    = "krig_sd",
                             "Co-kriging SE" = "cokrig_sd")) +
  tm_facets_wrap(nrow=1, sync=TRUE) +
  tm_shape(lead_sf) + tm_dots(fill="black", size=0.4) +
  tm_shape(meuse_sf |> slice(-idx)) + tm_dots(fill="grey40", size=0.3) + 
  tm_basemap("OpenStreetMap") +
  tm_layout(panel.labels = c("Ordinary kriging", "Co-kriging")) + 
  tm_view(use_WebGL = FALSE)

tmap_mode("view")
print(Map.pred)
print(Map.sd)


## -----------------------------------------------------------------------------
## 7) Cross-validation
## -----------------------------------------------------------------------------
# Ordinary kriging
cv_krig <- krige.cv(log(lead) ~ 1, lead_sf, model=vgm_lead, nfold=nrow(lead_sf),
                    debug.level=0)

# Co-kriging
cv_cokrig <- gstat.cv(lmc, nfold=nrow(lead_sf), remove.all=TRUE, debug.level=0)

cv_results <- data.frame(
  Method = c("Ordinary kriging", "Co-kriging"),
  MAE    = c(mean(abs(cv_krig$residual)), mean(abs(cv_cokrig$residual))),
  RMSE   = c(sqrt(mean(cv_krig$residual^2)), sqrt(mean(cv_cokrig$residual^2)))
)

cv_results |> 
  mutate(MAE = round(MAE, 3),
         RMSE = round(RMSE, 3))
