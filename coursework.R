library(rgdal)
library(tmap)
library(tmaptools)
library(OpenStreetMap)
library(sf)
library(readr)
library(ggplot2)
library(sp)
library(knitr)
library(spdep)
library(gridExtra)
library(gstat)
library(GSIF)
library(rgeos)
library(spgwr)
library(spatialreg)
tmap_mode("plot")


#set wroking directory
#setwd <- (" ")

#load shapefile
#wards <- st_read("~/Desktop/new_data/boundary/London_Ward_CityMerged.shp")
all_data <- st_read("~/Desktop/new_dat/all_data.shp")
#summary(all_data)
#colnames(all_data)
#plot(wardcrime)

#create base map 
#bmap <- read_osm(wardcrime, type="osm-transport")

#create map of crime rates (sum) per ward jenks
ttm()
tm_shape(all_data, title="Motor Vehicle Theft")+
  tm_polygons(col="veh_theft", style="jenks")+
  tm_scale_bar()+
  tm_compass()+
  tm_layout(legend.bg.color = "white", title = "Motor Vehicle Theft",  title.bg.color = "white")

# Histogram
ggplot(data = all_data, aes(all_data$veh_theft)) + geom_histogram()                        



# Morans I (spatial lag(ward boarders) against incident rates)
nb <- poly2nb(all_data) #Build a neighborhood list from polygons
W <- nb2mat(nb, style="W")
colnames(W) <- rownames(W)
kable(W[1:10,1:10], digits=3, caption="First 10 rows and columns of W for London wards", booktabs=T)

W1 <- nb2listw(nb)  #Convert neighbours list to spatial weights list
moran(all_data$veh_theft, W1, n=length(W1$neighbours), S0=Szero(W1))
moran(all_data$veh_theft, W1)
moran.test(all_data$veh_theft, W1) #Test with randomization to see if this a true correlation 99%
moran.mc(all_data$veh_theft, W1, nsim=999)  #Test using Monte-Carlo simulation
moran.plot(all_data$veh_theft, W1, xlab = "Incident Rates", ylab = "Spatially Lagged Incident Rates")

#local Moran's
Ii <- localmoran(all_data$veh_theft, W1)
all_data$Ii <- Ii[,"Ii"]
tm_shape(all_data) + tm_polygons(col="Ii", style="quantile")
tm_shape(all_data) + tm_polygons(col="total_cars", style="quantile")


#Regression
vehicle.shp <- readOGR("~/Desktop/new_dat/all_data.shp", "all_data")

#Linear regression
vehicle.lm <- lm(veh_theft ~ bame + road + X16_24_unem + med_hous_1 + med_house_ +
                   crime_rate + cars_per_h + open_space + num_jobs, data=vehicle.shp@data)
#Residuals
vehicle.shp$lm.res <- residuals(vehicle.lm)
#Plot map
tm_shape(vehicle.shp)+tm_polygons("lm.res", palette="-RdBu", style="quantile")

#Plot residuals
ggplot(data=vehicle.shp@data, aes(lm.res)) + geom_histogram()
ggplot(data=vehicle.shp@data, aes(sample=lm.res)) + geom_qq() + geom_qq_line()

# Moran test
vehicle.W <- nb2listw(poly2nb(vehicle.shp))
lm.morantest(vehicle.lm, vehicle.W)
#This indicates significant autocorrelation in the residuals. 
#However, it does not tell us the type of autocorrelation present.

#Lagrange multiplier diagnostics for spatial dependence
#help to decide whether to use a spatail error or spatial lag model
lm.LMtests(vehicle.lm, vehicle.W, test="RLMlag")
lm.LMtests(vehicle.lm, vehicle.W, test="RLMerr")


#  num_jobs

# Spatial lag model - better ???
vehicle.lag <- lagsarlm(veh_theft ~ road + med_hous_1 + med_house_ + cars_per_h + 
                          open_space + no_car + area_sq_km + pop +
                          mean_age + road + not_uk_bor, 
                        data=vehicle.shp@data, listw=vehicle.W)
summary(vehicle.lag)
vehicle.shp$lag.res <- residuals(vehicle.lag)
tm_shape(vehicle.shp)+tm_polygons("lag.res", palette="-RdBu", style="quantile")




# Spatial error model
vehicle.err <- errorsarlm(veh_theft ~ road + med_hous_1 + med_house_ + cars_per_h + 
                            open_space + no_car + area_sq_km + pop +
                            mean_age + road + not_uk_bor, data=vehicle.shp@data, listw=vehicle.W)
summary(vehicle.err)
vehicle.shp$err.res <- residuals(vehicle.err)
vehicle.shp$err.fit <- exp(fitted.values(vehicle.err))
tm_shape(vehicle.shp)+tm_polygons("err.res", palette="-RdBu", style="quantile")




# Spatial Durbin model
vehicle.durbin <- lagsarlm(veh_theft ~ road + med_hous_1 + med_house_ + cars_per_h + 
                             open_space + no_car + area_sq_km + pop +
                             mean_age + road + not_uk_bor, 
                          data=vehicle.shp@data, listw=vehicle.W, type="mixed")
summary(vehicle.durbin)
vehicle.shp$durbin.res <- residuals(vehicle.durbin)
vehicle.shp$durbin.fit <- exp(fitted.values(vehicle.durbin))
tm_shape(vehicle.shp)+tm_polygons("durbin.res", palette="-RdBu", style="quantile")




# Geographically weighted regression

# Transform the polygons to EPSG:26986 projected CRS so bandwidth can be defined in metres.
vehicle.pr <- spTransform(vehicle.shp, CRS("+proj=lcc +lat_1=42.68333333333333 +lat_2=41.71666666666667 +lat_0=41 +lon_0=-71.5 +x_0=200000 +y_0=750000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#Bandwidth
bwG <- gwr.sel(veh_theft ~ road + med_hous_1 + med_house_ + cars_per_h + 
                 open_space + no_car + area_sq_km + pop +
                 mean_age + road + not_uk_bor, 
               data = vehicle.pr, gweight = gwr.Gauss, verbose = FALSE, longlat=FALSE)
bwG #2443.869


# https://data.cdrc.ac.uk/system/files/practical10_0.html
# GWR model
vehicle.gwr <- gwr(veh_theft ~ road + med_hous_1 + med_house_ + cars_per_h + 
                     open_space + no_car + area_sq_km + pop +
                     mean_age + road + not_uk_bor, 
            data = vehicle.pr, bandwidth = bwG, gweight = gwr.Gauss, hatmatrix = TRUE, longlat=FALSE)
vehicle.gwr
gwr.morantest(vehicle.gwr, vehicle.W)

# Plot errors
gwr.errors <-vehicle.gwr$SDF$gwr.e
vehicle.shp$gwr.res <- gwr.errors
tm_shape(vehicle.shp)+tm_polygons("gwr.res", palette="-RdBu", style="quantile")


