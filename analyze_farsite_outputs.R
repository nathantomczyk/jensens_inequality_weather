library(terra)
library(sf)
library(xml2)
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyterra)

# 1. Load your FARSITE outputs
# (Assuming you exported these as GeoTIFFs from FARSITE/FlamMap)
flame_length <- rast("./farsite_sims/flame_lengths.tif")
arr_time   <- rast("./farsite_sims/arrival_times.tif")

arr_day<-1+(arr_time+1560) %/% 1440

wind<-read_sf("./farsite_sims/wind_vectors.shp")

wind_raster <- rasterize(wind, flame_length, field = "mph", fun = mean)

raw_wind <- read.table("./farsite_sims/weather_test_5.wxs", 
                       skip = 3, 
                       header = TRUE)



raw_wind$date<-ymd(paste(raw_wind$Year,raw_wind$Mth,raw_wind$Day,sep="-"))

raw_wind_mean<- raw_wind %>% group_by(date) %>% summarize(mean_wind=mean(WindSpd))

raw_wind_mean$day<-1:nrow(raw_wind_mean)

#

lookup_matrix <- as.matrix(raw_wind_mean[, c("day", "mean_wind")])

#

mean_wind_raster <- classify(arr_day, lookup_matrix)

# Stack them for easy extraction
results_stack_high_res <- c(flame_length, wind_raster)
names(results_stack_high_res) <- c("flame_length", "wind_spd")

results_stack_low_res <- c(flame_length, mean_wind_raster)
names(results_stack_low_res) <- c("flame_length", "wind_spd")

# 2. Create the 200m Sampling Grid
# method="regular" creates a systematic grid
sample_pts_high <- spatSample(results_stack_high_res, 
                         size = 5000, # Approximate number of points
                         method = "regular", 
                         as.points = TRUE, 
                         na.rm = TRUE)

sample_pts_low <- spatSample(results_stack_low_res, 
                              size = 5000, # Approximate number of points
                              method = "regular", 
                              as.points = TRUE, 
                              na.rm = TRUE)

# Note: To get exactly 200m, you can also create a custom grid:
# grid_template <- rast(results_stack, res=200)
# sample_pts <- as.points(grid_template)

# 3. Extract data to a dataframe
df_samples_high <- as.data.frame(sample_pts_high, geom = "XY")

df_samples_low <- as.data.frame(sample_pts_low, geom = "XY")

ggplot()+geom_spatraster(data=flame_length)+scale_fill_viridis_c(name="Flame length (ft)")
ggsave(filename="map_of_flame_lengths.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

ggplot()+geom_spatraster(data=wind_raster)+scale_fill_viridis_c(name="Wind Speed (mph)")
ggsave(filename="map_of_wind_high_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

ggplot()+geom_spatraster(data=mean_wind_raster)+scale_fill_viridis_c(name="Wind Speed (mph)")
ggsave(filename="map_of_wind_low_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

ggplot(df_samples_high,aes(x=wind_spd,y=flame_length))+geom_point()+ylab("Flame Length (ft)")+xlab("Wind Speed (mph)")+theme_classic()+theme(text=element_text(size=20))+geom_smooth(method="lm")
ggsave(filename="effect_wind_high_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

high_resolution_model<-lm(flame_length~wind_spd,df_samples_high)
summary(high_resolution_model)

ggplot(df_samples_low,aes(x=wind_spd,y=flame_length))+geom_point()+ylab("Flame Length (ft)")+xlab("Wind Speed (mph)")+theme_classic()+theme(text=element_text(size=20))+geom_smooth(method="lm")
ggsave(filename="effect_wind_low_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

low_resolution_model<-lm(flame_length~wind_spd,df_samples_low)
summary(low_resolution_model)
