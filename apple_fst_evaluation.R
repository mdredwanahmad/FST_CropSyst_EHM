################################################################################
# Fruit Surface Temperature (FST) model evaluation — Apple
# This script prepares observed FST and meteorological data from three apple
# sites (Sunrise, Quincy, and Snipes Rd.), runs the unified FST model, computes
# evaluation statistics, and generates all-day and peak-hour comparison figures.
################################################################################
library(tidyverse, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(viridisLite, warn.conflicts = FALSE)
library(viridis, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(readxl, warn.conflicts = FALSE)
library(rpart, warn.conflicts = FALSE)

### Set working directory
setwd("D:/Extreme_weather/FST_evaluation/")

################################################################################
# APPLE MODEL EVALUATION
################################################################################

### Sunrise orchard — August 2021
sunrise_aug_2021 <- fread("data/Sunrise_Aug2021_weather.csv",fill = TRUE, nrows = Inf)
sunrise_mod_FST <- fread("data/sunrise_combined_sensors.csv",fill = TRUE, nrows = Inf) %>% 
  dplyr::select(c("bin15","FST_max_madrm","FST_max_3sdrm")) %>% rename("Observed_FST"="FST_max_madrm")


### Combine Sunrise weather and observed fruit surface temperature data
Sunrise_Apple_combined_2021 <- left_join(sunrise_mod_FST,sunrise_aug_2021) %>% 
  # Convert meteorological variables to model units and add site/time variables
  mutate(Air_temp_C = (Air_temperature_15 - 32)*(5/9),
         Dew_point_c = (Dew_point - 32) * (5/9),
         Windspeed_ms = (Wind_Speed_2m*1.609*1000)/3600,
         lat = 47.3086,
         long = -120.0685,
         elevation = 267, 
         Station = 'Sunrise',
         Year = year(bin15),
         Jday = yday(bin15),
         Hour = hour(bin15),
         Year_day = paste0(Year, "_", Jday)) %>% na.omit() %>% ungroup() %>% 
  # Exclude days with incomplete records before model evaluation
  filter(!Jday %in% c(228,242)) %>% 
  # Run the unified FST model using the apple parameterization
  mutate(FST_simulated = unified_fruit_surface_temp(Air_temp_C, RH, Windspeed_ms,Jday,Solar_Radiation,Hour, lat, long, elevation,
                                                    fruit_diameter = 0.08,fruit_emissivity = 0.95, fruit_albedo = 0.6,r_fruit = 20000, 
                                                    day_id = Year_day))

### Quincy orchard — September 2021
# Prepare Quincy meteorological data
Quincy_weather<- fread("data/Quincy_Sep2021_weather.csv",fill = TRUE, nrows = Inf) %>% 
  mutate(Timestamps = mdy_hm(paste(Date,Hour))) %>% 
  arrange(Timestamps) %>% 
  # Retain meteorological variables required for FST simulation
  dplyr::select(c("Timestamps","Air_temperature_15","Dew_point","RH","Solar_Radiation","Wind_Speed_2m")) %>% 
  # Use a consistent timestamp variable name across evaluation sites
  rename('bin15'= "Timestamps")

# Prepare Quincy observed FST data
Quincy_FST <- fread("data/Quincy_combined_sensors.csv",fill = TRUE, nrows = Inf) %>% 
  ## select required variables and rename appropriately
  dplyr::select(c("bin15","FST_max_madrm","FST_max_3sdrm")) %>% rename("Observed_FST"="FST_max_madrm") %>% 
  mutate(bin15 = mdy_hm(bin15))



### Combine Quincy weather and observed FST data, then run the FST model
Quincy_Apple_combined_2021 <- left_join(Quincy_FST,Quincy_weather) %>% #left_join(Quincy_FST_v1,Quincy_sep_weather)
  ungroup() %>% 
  # Convert meteorological variables to model units and add site/time variables
  mutate(Air_temp_C = (Air_temperature_15 - 32)*(5/9),
         Dew_point_C = (Dew_point - 32) * (5/9),
         Windspeed_ms = (Wind_Speed_2m*1.609*1000)/3600,
         lat = 47.2255,
         long = -119.9568,
         elevation = 453, 
         Station = 'Quincy',
         Hour = hour(bin15),
         Jday = yday(bin15),
         Year = year(bin15),
         Year_day = paste0(Year, "_", Jday)) %>% na.omit() %>% ungroup() %>% 
  # Run the unified FST model using the apple parameterization
  mutate(FST_simulated = unified_fruit_surface_temp(Air_temp_C, RH, Windspeed_ms,Jday,Solar_Radiation,Hour, lat, long, elevation,
                                                    fruit_diameter = 0.08,fruit_emissivity = 0.95, fruit_albedo = 0.6,r_fruit = 20000, 
                                                    day_id = Year_day))


### Snipes Rd. commercial orchard — 2023
# Prepare Snipes Rd. observed FST and meteorological data
Snipes_apple <- fread("data/Snipes_combined_2023_control.csv",fill = TRUE, nrows = Inf) %>% 
  mutate(lat = 46.2946,
         long = -119.7250,
         elevation = 352,
         Station = 'Snipes_Rd.',
         RelHumidity = RelHumidity*100,
         RelHumidity_2 = RelHumidity_2*100,
         avg_dia_1 = avg_dia/1000,
         Year_day = paste0(Year, "_", Jday)) %>% ungroup() %>% 
  # Exclude days with incomplete records before model evaluation
  filter(!Jday %in% c(c(217,220,230,234))) %>% 
  # Variables ending in "_2" represent above-canopy measurements; these are used for the model evaluation
  mutate(FST_simulated = unified_fruit_surface_temp(AirTemp_2, RelHumidity_2, WindSpeed_2,Jday,Solar_2,Hour, lat, long, elevation,
                                                    fruit_diameter = avg_dia_1,fruit_emissivity = 0.95, fruit_albedo = 0.6,r_fruit = 20000, 
                                                    day_id = Year_day))



################################################################################
# EVALUATION FIGURES
################################################################################
##### 
# Create the output folder for evaluation figures if it does not already exist
if (!dir.exists("Figure")) {
  dir.create("Figure", recursive = TRUE)
}


### Representative time-series comparison for selected days at Snipes Rd.
# create Figure 5 (panel a)
data_2 <- Snipes_apple %>% filter(Jday %in% c(204:206))
P1<- ggplot(data_2)+
  geom_line( aes(x=bin15,y= FST_simulated,color = "Simulated apple FST °C"), group=1)+
  geom_line( aes(x=bin15,y= fruit_t25p,color = "Observed apple FST °C"), group=1)+
  scale_color_manual(values = c(
    "Observed apple FST °C" = "#00BA38","FST P95 °C"='purple',
    "Simulated apple FST °C" = "#0000FF", "FST top 20 mean °C" = "#FFD700", 'Air temperature °C'="#F8766D"))+
  labs(x= " ", y= "Apple fruit surface temperature °C", title = paste0("Snipes Apple"))+theme_bw()+ 
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 26), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 27,color = "black",angle = 0, hjust=0.5),
        axis.text.y = element_text(size = 27,color = "black",angle = 0, hjust=1),
        legend.position = "bottom", #c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_blank())

ggsave(P1, file=paste0("Figure/Snipes_apple.png"),width = 30, height = 22,dpi = 300, units = "cm")




################################################################################
# COMBINED APPLE DATASET AND EVALUATION METRICS
################################################################################
# Standardize variables from each site before combining datasets
# Sunrise: retain variables required for combined evaluation
sunrise_apple <- Sunrise_Apple_combined_2021 %>% 
  dplyr::select(c("bin15","Year","Hour","Jday","Observed_FST", "FST_max_3sdrm", "Station","FST_simulated","Air_temp_C",
                  "Windspeed_ms","RH","Solar_Radiation","lat","long","elevation","Station")) %>% 
  mutate(diameter = 0.08) %>% 
  rename("FST_max"="FST_max_3sdrm")

# Quincy: retain variables required for combined evaluation
quincy_apple <- Quincy_Apple_combined_2021 %>% 
  dplyr::select(c("bin15","Year","Hour","Jday","Observed_FST","FST_max_3sdrm", "Station","FST_simulated","Air_temp_C",
                  "Windspeed_ms","RH","Solar_Radiation","lat","long","elevation","Station")) %>% 
  mutate(diameter = 0.08) %>% 
  rename("FST_max"="FST_max_3sdrm")

# Snipes Rd.: retain variables required for combined evaluation and standardize variable names
Snipes_Rd <- Snipes_apple %>% 
  dplyr::select(c("bin15","Year","Hour","Jday","fruit_t25p","fruit_t20p", "Station","FST_simulated","AirTemp_2",
                  "WindSpeed_2","RelHumidity_2","Solar_2","lat","long","elevation","Station","avg_dia_1")) %>% 
  rename("Observed_FST"="fruit_t25p","FST_max"="fruit_t20p", "Air_temp_C"="AirTemp_2", "Windspeed_ms"="WindSpeed_2","RH"="RelHumidity_2","Solar_Radiation"="Solar_2", "diameter"="avg_dia_1")


### Define model-evaluation statistics
# Error is defined as simulated minus observed FST; therefore, positive mean
# bias indicates overprediction and negative mean bias indicates underprediction.
# The function reports RMSE, mean absolute error (stored as `mab`), mean bias,
# standard deviation of the errors, regression R², and squared correlation.
compute_stats <- function(obs, sim) {
  keep <- is.finite(obs) & is.finite(sim)
  obs <- obs[keep]; sim <- sim[keep]
  err <- sim - obs
  rmse <- sqrt(mean(err^2))
  mae  <- mean(abs(err))
  bias <- mean(err)
  sd_error <- sd(err)
  r2_lm <- summary(lm(sim ~ obs))$r.squared
  r2_corr <- cor(sim, obs)^2
  list(rmse = rmse, mab = mae, mb = bias, sd = sd_error, r2_lm = r2_lm, r2_corr = r2_corr)
}



# Combine the three apple evaluation sites into a single dataset
list_apple <- list(sunrise_apple,quincy_apple,Snipes_Rd)
combined_apple<- do.call(rbind,list_apple)

#### Station by station plot
#combined_apple_1 <- combined_apple %>% filter(Station=="Quincy", Hour %in% c(12:18))



# Compute model-performance statistics
Stats_apple <- compute_stats(combined_apple$Observed_FST, combined_apple$FST_simulated)
Stats_apple

# Format model-performance statistics for display on the figure
label_text_apple <- sprintf(
  "RMSE = %.2f °C\nMAB = %.2f °C\nMB = %.2f °C\nR²  = %.2f",
  Stats_apple$rmse, Stats_apple$mab, Stats_apple$mb, Stats_apple$r2_lm
)

# Create observed-versus-simulated FST comparison using 15-minute data
# create Figure 5 (panel c)
P1<- ggplot(combined_apple, aes(x = Observed_FST, y = FST_simulated))+
  geom_bin2d(bins = 70)+ #bins = 70
  scale_fill_continuous(type = "viridis")+ 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size=1) + 
  scale_x_continuous("Observed fruit surface temp. °C", breaks=seq(0,60, by=10))+
  scale_y_continuous("Simulated fruit surface temp. °C ", breaks=seq(0,60, by=10))+ 
  coord_equal(xlim = c(0, 60), ylim = c(0, 60)) +
  annotate("label", x = 0.7, y = 60,     # adjust position as you like
           hjust = 0, vjust = 1,
           label = label_text_apple,
           size = 9, linewidth = 0.25, fill = "white")+
  labs(x = " ", y = " ",
       title = "All day comparison")+theme_bw()+ 
  theme(
    axis.text = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.text = element_text(size = 20, color = "black"), 
    legend.title = element_text(size = 15, color = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.background = element_blank()
  )

ggsave(P1, file=paste0("Figure/combined_apple.png"), width = 30, height = 20,dpi = 300, units = "cm")



### Peak-hour evaluation (12:00–18:00, including only the 18:00 record)
peak_hours <- combined_apple %>% 
  filter(Hour %in% 12:17 | (Hour == 18 & format(bin15, "%H:%M:%S") == "18:00:00"))


# Compute model-performance statistics
Stats_apple <- compute_stats(peak_hours$Observed_FST, peak_hours$FST_simulated)
Stats_apple

# Format model-performance statistics for display on the figure
label_text_apple <- sprintf(
  "RMSE = %.2f °C\nMAB = %.2f °C\nMB = %.2f °C\nR²  = %.2f",
  Stats_apple$rmse, Stats_apple$mab, Stats_apple$mb, Stats_apple$r2_lm
)

# Create observed-versus-simulated FST comparison using 15-minute data
# create Figure 5 (panel e)
P1<- ggplot(peak_hours, aes(x = Observed_FST, y = FST_simulated))+
  geom_bin2d(bins = 70)+ #bins = 70
  scale_fill_continuous(type = "viridis")+ 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size=1) + 
  scale_x_continuous("Observed fruit surface temp. °C", breaks=seq(0,60, by=10))+
  scale_y_continuous("Simulated fruit surface temp. °C ", breaks=seq(0,60, by=10))+ 
  coord_equal(xlim = c(0, 60), ylim = c(0, 60)) +
  annotate("label", x = 0.7, y = 60,     # adjust position as you like
           hjust = 0, vjust = 1,
           label = label_text_apple,
           size = 9, linewidth = 0.25, fill = "white")+
  labs(x = " ", y = " ",
       title = "Combined peak hours (12:00 PM - 18:00 PM)")+theme_bw()+ 
  theme(
    axis.text = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.text = element_text(size = 20, color = "black"), 
    legend.title = element_text(size = 15, color = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.background = element_blank()
  )

ggsave(P1, file=paste0("Figure/combined_apple_peak_hours.png"), width = 30, height = 20,dpi = 300, units = "cm")

### Prepare the combined apple dataset for subsequent LHS/DELSA analyses
# Create the optional QA/QC output directory, if absent.
if (!dir.exists("Output_data")) {
  dir.create("Output_data", recursive = TRUE)
}

apple_combined_data <- combined_apple %>% 
  dplyr::select(-FST_simulated) %>% 
  mutate(day_id = paste0(Year,"_",Jday))

#write.csv(apple_combined_data, "Output_data/apple_combined_data.csv",row.names=F)
