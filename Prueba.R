## -------------------------------------------------------------------------------------------------------------
library(dplyr)
library(gstat)
library(sf)


## -------------------------------------------------------------------------------------------------------------
data(meuse, package="sp")         # soil samples
data(meuse.grid, package="sp")    # prediction grid

meuse_sf <- meuse |>
  as_tibble() |>
  st_as_sf(coords=c("x","y"), crs=28992)

grid_sf <- meuse.grid |>
  as_tibble() |>
  st_as_sf(coords=c("x","y"), crs=28992)


## -------------------------------------------------------------------------------------------------------------
set.seed(123)
idx <- sample(nrow(meuse_sf), 50)

lead_sf <- meuse_sf |> slice(idx)   # primary variable: 50 samples
zinc_sf <- meuse_sf                 # secondary variable: 155 samples




## -------------------------------------------------------------------------------------------------------------
v_lead <- variogram(log(lead) ~ 1, lead_sf)
vgm_lead <- fit.variogram(v_lead, model=vgm(psill=0.5, model="Sph",
                                            range=1000, nugget=0.05))
vgm_lead

pred_krig <- krige(log(lead) ~ 1, lead_sf, grid_sf, model=vgm_lead)
head(pred_krig)


## -------------------------------------------------------------------------------------------------------------
g <- gstat(NULL, id = "Pb", formula = log(lead) ~ 1, data = lead_sf)
g <- gstat(g,    id = "Zn", formula = log(zinc) ~ 1, data = zinc_sf)

v <- variogram(g)
v |> filter(dist==0)


## -------------------------------------------------------------------------------------------------------------
v <- v |> filter(dist>0)


## -------------------------------------------------------------------------------------------------------------
lmc <- fit.lmc(v, g, model=vgm_lead)
lmc$model


## -------------------------------------------------------------------------------------------------------------
#| fig-width: 7
#| fig-height: 6
#| fig-cap: "Direct and cross empirical semivariograms of log-lead (Pb, 50 samples) and log-zinc (Zn, 155 samples), with the fitted linear model of coregionalization."

plot(v, lmc, pch = 19, col = "black",
     xlab = "Distance (m)", ylab = "Semivariance",
     scales = list(alternating = 1, tck = c(1, 0)))


## -------------------------------------------------------------------------------------------------------------
pred_cokrig <- predict(lmc, grid_sf)
head(pred_cokrig)






## -------------------------------------------------------------------------------------------------------------
pred_krig$se <- sqrt(pred_krig$var1.var)
pred_cokrig$se <- sqrt(pred_cokrig$Pb.var)




## -------------------------------------------------------------------------------------------------------------
# Ordinary kriging
cv_krig <- krige.cv(log(lead) ~ 1, lead_sf, model=vgm_lead, nfold=nrow(lead_sf))

# Co-kriging
cv_cokrig <- gstat.cv(lmc, nfold=nrow(lead_sf), remove.all=TRUE, debug.level=0)


## -------------------------------------------------------------------------------------------------------------
cv_results <- data.frame(
  Method = c("Ordinary kriging", "Co-kriging"),
  MAE    = c(mean(abs(cv_krig$residual)), mean(abs(cv_cokrig$residual))),
  RMSE   = c(sqrt(mean(cv_krig$residual^2)), sqrt(mean(cv_cokrig$residual^2)))
)
cv_results |> 
  mutate(MAE = round(MAE, 3),
         RMSE=round(RMSE, 3))

