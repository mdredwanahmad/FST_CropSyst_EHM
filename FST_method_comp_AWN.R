################################################################################
# COMPARISON OF FRUIT SURFACE TEMPERATURE (FST) ESTIMATION METHODS — AWN DATA
#
# Compare daily maximum FST estimated from hourly meteorological observations,
# reconstructed hourly meteorological inputs, and two single-hour approximations
# (the hour of maximum air temperature or maximum solar radiation).
# Estimates are generated for apple and blueberry parameterizations.
#
# The analysis uses AgWeatherNet (AWN) inputs associated with the study of
# Savalkar et al. (2024). The underlying study dataset is not distributed with
# this repository; contact the study's corresponding author for data access.
#
# Workflow: run the large-data section on Kamiak HPC, export daily summaries,
# then reload the summary file locally to create the comparison figure.
################################################################################

# Load packages used for data processing, modeling, and visualization
library(agclimtools,warn.conflicts = FALSE)
library(tidyverse, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(viridisLite, warn.conflicts = FALSE)
library(viridis, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)

################################################################################
# PART 1 — AWN PREPROCESSING AND FST CALCULATIONS (KAMIAK HPC)
################################################################################
# The following steps require the full AWN input files and auxiliary lookup
# tables at the HPC paths specified below.

# Read one AWN CSV file and retain its source filename for traceability.
read_AWN <- function(flnm) {
  fread(flnm, fill = TRUE, nrows = Inf) %>%
    mutate(filename = flnm)
}

# Reconstruct hourly relative humidity by linear interpolation between daily
# RH extremes, using reconstructed hourly air temperature as the index.
humidity_melodist<- function(Hmax, Hmin, Tmax,Tmin, tempH){
  Rh <- Hmax+ ((tempH-Tmin)/(Tmax - Tmin))* (Hmin-Hmax)
  return(Rh)
}

# Load the station-year inventory used to retain AWN records with complete data.
#"/home/mdredwanahmad.khan/Extreme_weather_data/Solar_red_test/AWN_Complete_data_years.csv" 
# "D:/Agrivolts/AWN_Complete_data_years.csv"
AWN_complete_years <- fread("/home/mdredwanahmad.khan/Extreme_weather_data/Solar_red_test/AWN_Complete_data_years.csv",
                            fill = TRUE, nrows = Inf) %>% 
  rename("Year" = "Years_have_complete_data")

#"D:/Extreme_weather/Blue_berries/Melodist/AWN/Output/mean_course_stations.csv"
# Load normalized monthly/hourly temperature profiles for reconstruction.
mean_course <- read.csv("/home/mdredwanahmad.khan/Extreme_weather_data/apple_growing_grids_WA/Input_data/mean_course_stations.csv") %>% 
  rename("Hour"= "hour")


# Read and combine all AWN station files; match records to the complete-year
# inventory and retain the meteorological variables needed for FST calculations.
#path = "D:/Extreme_weather/Energy_balance/Solar_rad/AWN/" 
AWN_data <-
  list.files(path = "/weka/scratch/user/mdredwanahmad.khan/20261012_001525/Clened_all_vars/",pattern = "^AWN_",
             full.names = T, recursive = T) %>%
  map_df(~read_AWN(.)) %>% 
inner_join(AWN_complete_years %>% dplyr::select("STATION_NAME","LATITUDE","LONGITUDE","Year")) %>% 
  dplyr::select(c("Year","Month","Day","Hour","Jday","STATION_NAME","LATITUDE","LONGITUDE", "elevation_m","Air_temp_C","Relative_humidity",
                  "Windspeed_ms","SOLAR_RAD_WM2","Dew_point_c")) %>% 
  # Restrict analysis to May–September (months 5–9).
  filter(Month %in% c(5:9)) %>%
  group_by(STATION_NAME,Year,Month,Day,Jday,LATITUDE,LONGITUDE) %>% 
  mutate(Tmax = max(Air_temp_C),
         Tmin = min(Air_temp_C),
         Srad_total = sum(SOLAR_RAD_WM2),
         Rmax = max(Relative_humidity),
         Rmin = min(Relative_humidity)) %>% ungroup() %>% 
  left_join(mean_course %>% dplyr::select(c("Hour","Month","normalized_temp"))) %>%
  # Disaggregate hourly air temperature from daily Tmin/Tmax and the monthly
  # normalized diurnal temperature profile; derive corresponding hourly RH.
  group_by(STATION_NAME,Year,Month,Day,Jday,LATITUDE,LONGITUDE) %>% 
  mutate(Temp_dis = ((Tmax-Tmin)*normalized_temp+Tmin),
         Rh_dis = humidity_melodist(Rmax,Rmin,Tmax,Tmin,Temp_dis)) %>% 
  ungroup() %>% 
  ## Disaggregate Solar radiation from daily mean radiation
  mutate(Clear_sky_rad = clear_sky(Hour, Jday,LATITUDE,LONGITUDE, -120, elevation_m)) %>% 
  group_by(STATION_NAME,Year,Month,Day,Jday,LATITUDE,LONGITUDE) %>% 
  mutate(SumCS = sum(Clear_sky_rad, na.rm = TRUE),
         CS_fr= Clear_sky_rad/SumCS,
         SR_hr = Srad_total*CS_fr) %>% ungroup() %>% 
  # Run FST simulations for two crops with observed and disaggregated hourly
  # meteorological conditions. 
  mutate(day_id = paste0(Year,"_",Jday),
         FST = unified_fruit_surface_temp(
           Air_temp_C, Relative_humidity, Windspeed_ms, Jday, SOLAR_RAD_WM2, Hour,
           LATITUDE, LONGITUDE, elevation_m,
           F_beam = 1, F_d = 1, F_a = 1, F_g = 0, F_r = 0,
           fruit_diameter = 0.016,
           r_fruit = 17958,
           day_id = day_id        
         ),
         FST_diss_berry = unified_fruit_surface_temp(
           Temp_dis, Rh_dis, Windspeed_ms, Jday, SR_hr, Hour,
           LATITUDE, LONGITUDE, elevation_m,
           F_beam = 1, F_d = 1, F_a = 1, F_g = 0, F_r = 0,
           fruit_diameter = 0.016,
           r_fruit = 17958,
           day_id = day_id        
         ),
         FST_apple = unified_fruit_surface_temp(
           Air_temp_C, Relative_humidity, Windspeed_ms, Jday, SOLAR_RAD_WM2, Hour,
           LATITUDE, LONGITUDE, elevation_m,
           F_beam = 1, F_d = 1, F_a = 1, F_g = 0, F_r = 0,
           fruit_diameter = 0.08, fruit_emissivity = 0.95, fruit_albedo = 0.6,
           method = "HS",r_fruit = 20000,
           day_id = day_id        
         ),
         FST__diss_apple = unified_fruit_surface_temp(
           Temp_dis, Rh_dis, Windspeed_ms, Jday, SR_hr, Hour,
           LATITUDE, LONGITUDE, elevation_m,
           F_beam = 1, F_d = 1, F_a = 1, F_g = 0, F_r = 0,
           fruit_diameter = 0.08, fruit_emissivity = 0.95, fruit_albedo = 0.6,
           method = "HS",r_fruit = 20000,
           day_id = day_id        
         ))



################################################################################
# DAILY FST MAXIMA AND SINGLE-HOUR APPROXIMATIONS
################################################################################
# Summarize simulated hourly FST to daily maxima and identify the hours of
# maximum FST, air temperature, and solar radiation. Also estimate FST using
# meteorological conditions at the hour of Tmax or solar-radiation maximum.
# Both single-hour estimates use the daily mean wind speed.
Daily_AWN <- AWN_data %>% 
  group_by(STATION_NAME,Year,Month,Day,Jday,day_id, LATITUDE, LONGITUDE, elevation_m) %>% 
  summarise(FST_max = max(FST),
            FST_max_diss = max(FST_diss_berry),
            FST_max_apple = max(FST_apple),
            FST_max_apple_diss = max(FST__diss_apple),
            FST_max_hour = Hour[which.max(FST)],
            FST_max_hour_apple = Hour[which.max(FST_apple)],
            hour_Tmax = Hour[which.max(Air_temp_C)],
            hour_sradmax = Hour[which.max(SOLAR_RAD_WM2)],
            Tmax = max(Air_temp_C, na.rm=T),
            Tmin = min(Air_temp_C, na.rm=T),
            Temp_at_Sr_max = Air_temp_C[which.max(SOLAR_RAD_WM2)],
            Srad_max = max(SOLAR_RAD_WM2, na.rm=T),
            Srad_at_Tmax = SOLAR_RAD_WM2[which.max(Air_temp_C)],
            wind= mean(Windspeed_ms),
            #Rh_min = min(Relative_humidity,na.rm = T),
            Rh_atTmax = Relative_humidity[which.max(Air_temp_C)],
            Rh_atSrmax = Relative_humidity[which.max(SOLAR_RAD_WM2)],
            ##SR_at_tmax_temp_max
            FST_Tmax_SR = unified_fruit_surface_temp(Tmax,Rh_atTmax,wind,Jday,Srad_at_Tmax,hour_Tmax,LATITUDE,LONGITUDE,elevation_m,
                                                     fruit_diameter = 0.016,r_fruit = 17958,day_id = day_id),
            ##SR_tmax_temp_at_SR_max
            FST_SRmax_T = unified_fruit_surface_temp(Temp_at_Sr_max,Rh_atSrmax,wind,Jday,Srad_max,hour_sradmax,LATITUDE,LONGITUDE,elevation_m,
                                                     fruit_diameter = 0.016,r_fruit = 17958,day_id = day_id),
            ### Apple model
            ##SR_at_tmax_temp_max
            FST_Tmax_SR_apple = unified_fruit_surface_temp(Tmax,Rh_atTmax,wind,Jday,Srad_at_Tmax,hour_Tmax,LATITUDE,LONGITUDE,elevation_m,
                                                           fruit_diameter = 0.08, fruit_emissivity = 0.95, fruit_albedo = 0.6,
                                                           method = "HS",r_fruit = 20000,
                                                           day_id = day_id ),
            ##SR_tmax_temp_at_SR_max
            FST_SRmax_T_apple = unified_fruit_surface_temp(Temp_at_Sr_max,Rh_atSrmax,wind,Jday,Srad_max,hour_sradmax,LATITUDE,LONGITUDE,elevation_m,
                                                           fruit_diameter = 0.08, fruit_emissivity = 0.95, fruit_albedo = 0.6,
                                                           method = "HS",r_fruit = 20000,
                                                           day_id = day_id ))

# Export daily FST summaries after completing the HPC run.
# This data set is provided in repository so the following section can be reproduced
#write.csv(Daily_AWN, "/home/mdredwanahmad.khan/Extreme_weather_data/Solar_red_test/calc_FST_daily_AWN.csv", row.names = F)

################################################################################
# PART 2 — LOCAL ANALYSIS OF DAILY HPC OUTPUT
################################################################################
# Run this section after exporting the daily summary from Kamiak.
# Reload the precomputed file for the methodological comparison and figure.
# This dataset is provided in repository, so the following analysis can be reproduced
Daily_AWN<- fread("data/calc_FST_daily_AWN.csv", fill=TRUE, nrows= Inf)

# Group the corresponding apple and blueberry variables by estimation
# method so that both crops receive consistent labels in the plot. 
Tmax_hour <- c("FST_Tmax_SR_b","FST_Tmax_SR_apple")
SRmax_hour <- c("FST_SRmax_T_b","FST_SRmax_T_apple")
Max_fst <- c('FST_max_b',"FST_max_apple")
Max_dis_fst <- c('FST_max_diss_b',"FST_max_apple_diss")

# Reshape the daily summary to long format and assign crop and method labels.
# The saved daily file uses *_b suffixes for blueberry estimates.
solar_rad_v1 <- Daily_AWN %>% #filter(Month %in% c(5:9)) %>% 
  dplyr::select(c('STATION_NAME','Year', 'Month', 'Day', 'Jday','LATITUDE', 'LONGITUDE',
                  'FST_max_b','FST_max_diss_b',"FST_max_apple","FST_max_apple_diss",
                  "FST_Tmax_SR_b","FST_SRmax_T_b","FST_Tmax_SR_apple","FST_SRmax_T_apple")) %>% ungroup() %>% 
  reshape2::melt(id = c('STATION_NAME','Year', 'Month', 'Day', 'Jday','LATITUDE', 'LONGITUDE')) %>% 
  mutate(fruit = ifelse(variable %in% c('FST_max_b','FST_max_diss_b','FST_Tmax_SR_b','FST_SRmax_T_b'),"blueberry","apple"),
         variable = case_when(variable %in% Tmax_hour ~ "FST at hourly Tmax",
                              variable %in% SRmax_hour ~ "FST at hourly SRmax",
                              variable %in% Max_fst ~ "Maximum FST (observed hourly)",
                              variable %in% Max_dis_fst ~ "Maximum FST (disaggregated hourly)"))

# Plot monthly distributions of daily FST estimates for each method and crop
# (Appendix Figure C02).
P1<- ggplot(solar_rad_v1, aes(x=as.character(Month),y= value, fill = variable))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom = 'errorbar') + facet_wrap(~fruit)+ 
  labs(x= " ", y= "Fruit surface temperature °C", title = paste0(" "))+ theme_bw()+
  scale_y_continuous("Fruit surface temperature °C",limits=c(0,60), breaks=seq(0,60, by=10))+
  scale_fill_manual(values = c(
    "FST at hourly Tmax" = "#619CFF",  "Maximum FST (disaggregated hourly)" = "#00BA38",       
    "FST at hourly SRmax" = "#FF0000","Maximum FST (observed hourly)" = "#0000FF",
    "FST_Tmax_SR_apple" = "#619CFF",  "FST_max_apple_diss" = "#00BA38",       
    "FST_SRmax_T_apple" = "#FF0000","FST_max_apple" = "#0000FF", "FST_max_apple_model" = "black"))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 22,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 22,color = "black",angle = 0, hjust=1),
        legend.position = "right", 
        legend.justification = c(0, 1))

ggsave(P1, file=paste0("Figure/FST_estimation_methods_1.png"), width = 30, height = 20,dpi = 300, units = "cm")



