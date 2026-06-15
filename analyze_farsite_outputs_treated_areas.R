library(terra)
library(sf)
library(xml2)
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(betareg)

# 1. Load your FARSITE outputs
# (Assuming you exported these as GeoTIFFs from FARSITE/FlamMap)
flame_length <- rast("C:/Fire_sandbox/flame_lengths_treated.tif")
arr_time   <- rast("C:/Fire_sandbox/arrival_times_treated.tif")

arr_day<-1+(arr_time+1560) %/% 1440

wind<-read_sf("C:/Fire_sandbox/wind_vectors_treated.shp")

wind_raster <- rasterize(wind, flame_length, field = "mph", fun = mean)

raw_wind <- read.table("C:/Fire_sandbox/weather_test_5.wxs", 
                       skip = 3, 
                       header = TRUE)

raw_wind$date<-ymd(paste(raw_wind$Year,raw_wind$Mth,raw_wind$Day,sep="-"))

raw_wind_mean<- raw_wind %>% group_by(date) %>% summarize(mean_wind=mean(WindSpd))

raw_wind_mean$day<-1:nrow(raw_wind_mean)

#



lookup_matrix <- as.matrix(raw_wind_mean[, c("day", "mean_wind")])

#

mean_wind_raster <- classify(arr_day, lookup_matrix)

##### treatment

# 1. Define your threshold (e.g., 50% canopy cover)
fuel_conditions<-rast("C:/Fire_sandbox/Sandbox_landscape_larger_treated.tif")

threshold <- 50

treatments <- fuel_conditions$canopy < threshold


# Stack them for easy extraction
results_stack_high_res <- c(flame_length, wind_raster,treatments)
names(results_stack_high_res) <- c("flame_length", "wind_spd","treatment")

results_stack_low_res <- c(flame_length, mean_wind_raster,treatments)
names(results_stack_low_res) <- c("flame_length", "wind_spd","treatment")

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



# 3. Extract data to a dataframe
df_samples_high <- as.data.frame(sample_pts_high, geom = "XY")

df_samples_low <- as.data.frame(sample_pts_low, geom = "XY")

ggplot()+geom_spatraster(data=flame_length)+scale_fill_viridis_c(name="Flame length (ft)")
ggsave(filename="flame_length_treated.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

#ggplot()+geom_spatraster(data=wind_raster)+scale_fill_viridis_c(name="Wind Speed (mph)")

#ggplot()+geom_spatraster(data=mean_wind_raster)+scale_fill_viridis_c(name="Wind Speed (mph)")

###
max_flame<-max(df_samples_high$flame_length)

df_samples_high$flame_beta <- df_samples_high$flame_length / max_flame

# Ensure no absolute 0s or 1s (Beta distribution constraint)
df_samples_high$flame_beta <- ifelse(df_samples_high$flame_beta <= 0, 0.001, df_samples_high$flame_beta)
df_samples_high$flame_beta <- ifelse(df_samples_high$flame_beta >= 1, 0.999, df_samples_high$flame_beta)

df_samples_low$flame_beta <- df_samples_high$flame_length / max_flame

# Ensure no absolute 0s or 1s (Beta distribution constraint)
df_samples_low$flame_beta <- ifelse(df_samples_low$flame_beta <= 0, 0.001, df_samples_low$flame_beta)
df_samples_low$flame_beta <- ifelse(df_samples_low$flame_beta >= 1, 0.999, df_samples_low$flame_beta)


hr_model_beta<-betareg(flame_beta~wind_spd * treatment, 
                       data = df_samples_high)
summary(hr_model_beta)

lr_model_beta<-betareg(flame_beta~wind_spd * treatment, 
                       data = df_samples_low)
summary(lr_model_glm)


#####

plot_beta_predictions <- function(model, data, max_val = 180, title = "Predicted Flame Lengths (Beta Model)") {
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  
  # 1. Argument order safety check
  if(!is.numeric(max_val)) stop("max_val must be a number.")
  
  # 2. Create wind sequence
  wind_seq <- seq(min(data$wind_spd, na.rm = TRUE), 
                  max(data$wind_spd, na.rm = TRUE), 
                  length.out = 100)
  
  # 3. Get predictions and force a consistent column name using 'value'
  # We use regrid to ensure we are on the 0-1 response scale
  spec <- emmeans(model, ~ treatment | wind_spd, 
                  at = list(wind_spd = wind_seq)) %>%
    regrid(transform = "response")
  
  # 4. Convert to df and rename the estimate column to "prob" for consistency
  pred_grid <- as.data.frame(spec)
  
  # Identify the estimate column (it's the first one that isn't a factor or wind_spd)
  # In emmeans, the estimate is usually the column 'prob', 'response', or the variable name.
  # We'll find it by position or name.
  if("prob" %in% names(pred_grid)) {
    names(pred_grid)[names(pred_grid) == "prob"] <- "estimate"
  } else if ("response" %in% names(pred_grid)) {
    names(pred_grid)[names(pred_grid) == "response"] <- "estimate"
  } else {
    # If it's named after your variable (e.g. 'flame_beta'), it's usually column 3 or 4
    # This identifies the numeric column that isn't wind_spd
    names(pred_grid)[sapply(pred_grid, is.numeric) & names(pred_grid) != "wind_spd"][1] <- "estimate"
  }
  
  # 5. Rescale and find the CI columns (emmeans uses asymp.LCL/UCL for beta)
  pred_grid <- pred_grid %>%
    mutate(
      fit = estimate * max_val,
      lwr = asymp.LCL * max_val,
      upr = asymp.UCL * max_val
    )
  
  # 6. Create the plot
  ggplot(pred_grid, aes(x = wind_spd, y = fit, color = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(data = data, aes(y = flame_length), alpha = 0.05, size = 0.5) +
    geom_hline(yintercept = max_val, linetype = "dashed", color = "gray50") +
    scale_color_viridis_d(begin = 0.2, end = 0.8) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8) +
    labs(
      title = title,
      subtitle = paste("Beta regression scaled to 169 ft ceiling"),
      x = "Wind Speed (mph)",
      y = "Flame Length (ft)",
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal()
}

ggplot(df_samples_low,aes(x=wind_spd,y=flame_length,color=treatment))+geom_point()+geom_smooth()+scale_color_manual(values=c("firebrick","steelblue"))+theme_classic()+theme(text=element_text(size=20))+ylab("Flame Length (ft)")+xlab("Wind Speed (mph)")
ggsave(filename="effect_wind_and_treatment_high_res_new.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

ggplot(df_samples_high,aes(x=wind_spd,y=flame_length,color=treatment))+geom_point()+geom_smooth()+scale_color_manual(values=c("firebrick","steelblue"))+theme_classic()+theme(text=element_text(size=20))+ylab("Flame Length (ft)")+xlab("Wind Speed (mph)")
ggsave(filename="effect_wind_and_treatment_low_res_new.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)


ggplot(df_samples_high,aes(x=wind_spd,y=flame_length,color=treatment))+geom_point()+geom_smooth()

# Plot the High-Res model
plot_beta_predictions(hr_model_beta, df_samples_high,max_flame, "High Resolution")
ggsave(filename="effect_wind_and_treatment_high_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

# Plot the Low-Res model to compare the interaction slopes
plot_beta_predictions(lr_model_beta, df_samples_low,max_flame, "Low Resolution")
ggsave(filename="effect_wind_and_treatment_low_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

###############################################################################################################

estimate_stand_mortality <- function(flame_length_m) {
  # We convert Flame Length to a 'Severity Index'
  # Ponderosa stands typically transition to high mortality 
  # when flames exceed 2.5 to 3.5 meters (~8-12 ft).
  
  # Logistic parameters: 
  # k = steepness of the transition
  # x0 = the 'inflection point' where 50% mortality occurs
  k <- 1.5
  x0 <- 3.0 
  
  stand_mortality <- 1 / (1 + exp(-k * (flame_length_m - x0)))
  
  return(stand_mortality)
}

df_samples$mortaility<-estimate_stand_mortality(df_samples$`Farsite Flame Length`)

plot(df_samples$mean,df_samples$mortaility)

df_samples_high$mortality<-estimate_stand_mortality(df_samples_high$flame_length)

plot(df_samples_high$flame_length,df_samples_high$mortality)

ggplot(df_samples_high,aes(x=wind_spd,y=mortality,color=treatment))+geom_smooth()


df_samples_low$mortality<-estimate_stand_mortality(df_samples_low$flame_length)


ggplot(df_samples_low,aes(x=wind_spd,y=mortality,color=treatment))+geom_smooth()


df_samples_high$mortality_beta<-df_samples_high$mortality-0.001


hr_model_beta<-betareg(mortality_beta~wind_spd * treatment, 
                       data = df_samples_high)
summary(hr_model_beta)

df_samples_low$mortality_beta<-df_samples_low$mortality-0.001

lr_model_beta<-betareg(mortality_beta~wind_spd * treatment, 
                       data = df_samples_low)
summary(lr_model_glm)

######

plot_beta_predictions <- function(model, data, max_val = 1, title = "Predicted Flame Lengths (Beta Model)") {
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  
  # 1. Argument order safety check
  if(!is.numeric(max_val)) stop("max_val must be a number.")
  
  # 2. Create wind sequence
  wind_seq <- seq(min(data$wind_spd, na.rm = TRUE), 
                  max(data$wind_spd, na.rm = TRUE), 
                  length.out = 100)
  
  # 3. Get predictions and force a consistent column name using 'value'
  # We use regrid to ensure we are on the 0-1 response scale
  spec <- emmeans(model, ~ treatment | wind_spd, 
                  at = list(wind_spd = wind_seq)) %>%
    regrid(transform = "response")
  
  # 4. Convert to df and rename the estimate column to "prob" for consistency
  pred_grid <- as.data.frame(spec)
  
  # Identify the estimate column (it's the first one that isn't a factor or wind_spd)
  # In emmeans, the estimate is usually the column 'prob', 'response', or the variable name.
  # We'll find it by position or name.
  if("prob" %in% names(pred_grid)) {
    names(pred_grid)[names(pred_grid) == "prob"] <- "estimate"
  } else if ("response" %in% names(pred_grid)) {
    names(pred_grid)[names(pred_grid) == "response"] <- "estimate"
  } else {
    # If it's named after your variable (e.g. 'flame_beta'), it's usually column 3 or 4
    # This identifies the numeric column that isn't wind_spd
    names(pred_grid)[sapply(pred_grid, is.numeric) & names(pred_grid) != "wind_spd"][1] <- "estimate"
  }
  
  # 5. Rescale and find the CI columns (emmeans uses asymp.LCL/UCL for beta)
  pred_grid <- pred_grid %>%
    mutate(
      fit = estimate * max_val,
      lwr = asymp.LCL * max_val,
      upr = asymp.UCL * max_val
    )
  
  # 6. Create the plot
  ggplot(pred_grid, aes(x = wind_spd, y = fit, color = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(data = data, aes(y = flame_length), alpha = 0.05, size = 0.5) +
    geom_hline(yintercept = max_val, linetype = "dashed", color = "gray50") +
    scale_color_viridis_d(begin = 0.2, end = 0.8) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8) +
    labs(
      title = title,
      x = "Wind Speed (mph)",
      y = "Proportional mortality",
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal()
}


# Plot the High-Res model
plot_beta_predictions(hr_model_beta, df_samples_high,max_val = 1, "High Resolution")
ggsave(filename="effect_wind_and_treatment_high_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)

# Plot the Low-Res model to compare the interaction slopes
plot_beta_predictions(lr_model_beta, df_samples_low,max_val = 1, "Low Resolution")
ggsave(filename="effect_wind_and_treatment_low_res.tif",plot=last_plot(),units="in",height=8,width=8,dpi=600)
