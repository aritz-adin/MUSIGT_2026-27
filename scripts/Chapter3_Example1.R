## -----------------------------------------------------------------------------
## 1) Objectives:
##    - Explore two spatially correlated variables: the concentrations of
##      lead (primary variable) and zinc (secondary variable) in the Meuse
##      soil samples.
##    - Compute the direct semivariograms of both variables and their
##      cross-semivariogram.
##    - Fit a linear model of coregionalization, the joint model required
##      by co-kriging.
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## 2) Packages
## -----------------------------------------------------------------------------
library(dplyr)
library(gstat)
library(sf)


## -----------------------------------------------------------------------------
## 3) Description of the data
## -----------------------------------------------------------------------------
data(meuse, package="sp")        # soil samples
data(meuse.riv, package="sp")    # outline of the river

meuse_sf <- meuse |> 
  as_tibble() |> 
  st_as_sf(coords=c("x","y"), crs=28992)

river <- list(meuse.riv) |>
  st_polygon() |>
  st_sfc(crs = 28992)

head(meuse_sf)
head(river)


## Log-transformation of lead and zinc variables -------------------------------
par(mfrow = c(2, 2))
hist(meuse_sf$lead, main = "Lead", xlab = "ppm", col = "grey80")
hist(meuse_sf$zinc, main = "Zinc", xlab = "ppm", col = "grey80")
hist(log(meuse_sf$lead), main = "log(Lead)", xlab = "log(ppm)", col = "grey80")
hist(log(meuse_sf$zinc), main = "log(Zinc)", xlab = "log(ppm)", col = "grey80")


## -----------------------------------------------------------------------------
## 4) Semivariogram of the primary variable
## -----------------------------------------------------------------------------
## Empirical semivariogram -----------------------------------------------------
v_lead <- variogram(log(lead) ~ 1, meuse_sf)
v_lead


## Fitted theoretical model ----------------------------------------------------
vgm_lead <- fit.variogram(v_lead, model=vgm(psill=0.5, model="Sph",
                                            range=1000, nugget=0.05))
vgm_lead


## Plot of the empirical semivariogram and fitted model ------------------------
plot(v_lead, vgm_lead, pch = 19,
     xlab = "Distance (m)", ylab = "Semivariance",
     main = "Semivariogram of log-lead (Pb)")


## -----------------------------------------------------------------------------
## 5) Direct and cross-empirical semivariograms
## -----------------------------------------------------------------------------
g <- gstat(NULL, id = "Pb", formula = log(lead) ~ 1, data = meuse_sf)
g <- gstat(g,    id = "Zn", formula = log(zinc) ~ 1, data = meuse_sf)
g


## Compute the three empirical semivariograms ----------------------------------
v <- variogram(g)

v |> 
  filter(id=="Pb") |> 
  select(np, dist, gamma) |> 
  head()

v |> 
  filter(id=="Pb.Zn") |> 
  select(np, dist, gamma) |> 
  head()


## -----------------------------------------------------------------------------
## 6) Linear model of coregionalization
## -----------------------------------------------------------------------------
lmc <- fit.lmc(v, g, model=vgm_lead)
lmc$model


## Plot of the fitted model ----------------------------------------------------
plot(v, lmc, pch = 19, col = "black",
     xlab = "Distance (m)", ylab = "Semivariance",
     scales = list(alternating = 1, tck = c(1, 0)))
