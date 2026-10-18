################################################################################
# FRUIT SURFACE TEMPERATURE (FST) — SENSOR QA/QC AND DISTRIBUTION ANALYSIS
# Sites: Sunrise and Quincy apple orchards (Washington), and Corvallis blueberry.
# Workflow: prepare thermocouple data; characterize hourly FST distributions;
# screen unusually high sensor readings using 3-SD and MAD criteria, and finally choose based on MAD;
# compare original and screened FST summaries; and assess distribution normality.
# Note: Flagging individual observations and excluding entire sensor columns
# are distinct steps in the workflow below.
################################################################################

# Load required R packages.
library(tidyverse, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(viridisLite, warn.conflicts = FALSE)
library(viridis, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(reshape2,warn.conflicts = FALSE)
library(readxl,warn.conflicts = FALSE)
library(patchwork,warn.conflicts = FALSE)
library(e1071,warn.conflicts = FALSE)
library(GGally,warn.conflicts = FALSE)

################################################################################
# SUNRISE APPLE — SENSOR DATA PREPARATION AND QUALITY CONTROL
################################################################################
## Observed Sunrise Apple FST data 
sunrise_aug_FST2021_test <- fread("data/Honeycrisp_2021_Fruit_temperature_combined_sunrise.csv", fill = TRUE, nrows = Inf) %>% 
  mutate(Timestamp = mdy_hms(paste(date, time)),
         bin15 = floor_date(Timestamp, "15 minutes")) %>% 
  dplyr::select(c("bin15","Timestamp","1ch","2ch","3ch","4ch","Data_logger")) %>% 
  rename("channel_1"="1ch","channel_2"="2ch","channel_3"="3ch","channel_4"="4ch")

###### There were only two repeating rows/timesteps so took average across timesteps and data logger 
sunrise_15min <- sunrise_aug_FST2021_test %>%
  group_by(bin15, Data_logger) %>%
  summarise(
    channel_1 = mean(channel_1, na.rm = TRUE),
    channel_2 = mean(channel_2, na.rm = TRUE),
    channel_3 = mean(channel_3, na.rm = TRUE),
    channel_4 = mean(channel_4, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape thermocouple readings to long format for hourly distribution plots.
# Retain readings at or below 60 °C for this visualization.
sunrise_melt <- sunrise_15min %>% 
  reshape2::melt(id.vars = c("bin15", "Data_logger")) %>% 
  mutate(Hour = hour(bin15),
         Hour = factor(hour(bin15),levels = 0:23)) %>% filter(value <= 60)

# Load Sunrise meteorological observations.
sunrise_aug_2021 <- fread("data/Sunrise_Aug2021_weather.csv",fill = TRUE, nrows = Inf)
# Join FST and air temperature by 15-minute timestamp; calculate FST minus air temperature.
sunrise_combined <- left_join(sunrise_melt,sunrise_aug_2021 %>% dplyr::select(c("bin15","Air_temperature_15"))) %>% 
  mutate(Air_temp_C = (Air_temperature_15 - 32)*(5/9),
         diff = value-Air_temp_C)

# Show the distribution of individual thermocouple readings by hour of day.
# Create Appendix Figure C01 (panel a)
P1<- ggplot(sunrise_combined, aes(x=Hour,y= value))+ 
  geom_boxplot()+ stat_boxplot(geom='errorbar')+  
  labs(x= "Hour (0:23)", y=  "Apple fruit surface temperature °C ", title = paste0("Sunrise Apple"))+ theme_bw()+
  scale_fill_manual(values = c("apple" = "#D73027","blueberry" = "#4575B4"))+
  scale_x_discrete(drop = FALSE, breaks = c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22"))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 24,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 24,color = "black",angle = 0, hjust=1),
        legend.position = "bottom") 

ggsave(P1, file=paste0("Figure/Sunrise_FST_spread.png"), width = 30, height = 20,dpi = 300, units = "cm")


# Load the wide-format sensor dataset used for timestep-level QA/QC.
sunrise_wide<- read.csv("data/sunrise_wide.csv")

# Identify thermocouple columns for cross-sensor comparisons.
sensor_cols <- grep("^channel_",names(sunrise_wide),value = TRUE)

# Calculate cross-sensor statistics independently for each 15-minute timestamp.
# The 3-SD and MAD screens flag unusually high readings relative to other sensors.
sunrise_QC_long <- sunrise_wide %>%
  mutate(bin15 = ymd_hms(bin15)) %>% 
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "FST"
  ) %>%
  group_by(bin15) %>%
  mutate(
    FST_mean = mean(FST, na.rm = TRUE),
    FST_median = median(FST, na.rm = TRUE),
    FST_sd   = sd(FST, na.rm = TRUE),
    
    upper_3SD = FST_mean + 3 * FST_sd,
    FST_max = max(FST, na.rm=T),
    FST_max_3sd = FST_max-upper_3SD,
    FST_MAD = mad(FST,constant = 1.4826, na.rm = TRUE),
    
    z_score = (FST - FST_mean) / FST_sd,
    robust_z =(FST - FST_median) / FST_MAD,
    flag_3SD = FST > upper_3SD,
    
    flag_MAD = robust_z > 3.5
  ) %>%
  ungroup()

# Extract readings flagged by BOTH screening criteria for inspection (not removal here).
sunrise_qc <- sunrise_QC_long %>% filter(flag_MAD==TRUE & flag_3SD==TRUE)

# Create the optional QA/QC output directory, if absent.
if (!dir.exists("Output_data")) {
  dir.create("Output_data", recursive = TRUE)
}
#write.csv(sunrise_qc, 'Output_data/sunrise_quality_issue.csv',row.names = F)

# Sensor-exclusion lists identified from the respective 3-SD and MAD assessments.
# These lists are applied to subsequent maximum-FST summaries, not removed here.
sensors_3sd <- c("channel_4_O","channel_4_S","channel_4_C", "channel_3_P","channel_3_S","channel_2_S","channel_2_B","channel_1_S","channel_3_B","channel_1_J")
sensor_mad <- c("channel_1_O","channel_1_P","channel_1_S","channel_2_B","channel_2_P","channel_2_S","channel_3_A","channel_3_B","channel_3_P","channel_3_S","channel_4_C",
                "channel_4_I","channel_4_O","channel_4_S","channel_1_J")


# Extract all readings flagged by the MAD criterion for further review.
# we used the MAD based approach for final exclusion of flagged sesnors
sunrise_qc_mad <- sunrise_QC_long %>% filter(flag_MAD==TRUE)
#write.csv(sunrise_qc_mad, 'Output_data/sunrise_quality_issue_MAD.csv',row.names = F)

# Summarize each timestamp across all thermocouples, then recompute maximum FST
# after excluding sensor IDs on each screening list.
sunrise_wide1 <- sunrise_wide %>% 
  mutate(bin15 = ymd_hms(bin15)) %>% 
  rowwise() %>%
  mutate(
    # --- Statistics using ALL channel columns ---
    FST_mean   = mean(c_across(starts_with("channel_")), na.rm = TRUE),
    FST_median = median(c_across(starts_with("channel_")), na.rm = TRUE),
    FST_sd     = sd(c_across(starts_with("channel_")), na.rm = TRUE),
    FST_P75    = quantile(c_across(starts_with("channel_")), probs = 0.75, na.rm = TRUE),
    FST_P85    = quantile(c_across(starts_with("channel_")), probs = 0.85, na.rm = TRUE),
    FST_P95    = quantile(c_across(starts_with("channel_")), probs = 0.95, na.rm = TRUE),
    FST_max    = max(c_across(starts_with("channel_")), na.rm = TRUE),
    
    # --- Statistics EXCLUDING specific sensor vectors (Columns stay in dataset) ---
    # Example using sensors_3sd (change vector name or combine as needed)
    FST_max_3sdrm    = max(c_across(setdiff(starts_with("channel_"), any_of(sensors_3sd))), na.rm = TRUE),
    FST_max_madrm    = max(c_across(setdiff(starts_with("channel_"), any_of(sensor_mad))), na.rm = TRUE),
    
    # --- Difference calculations ---
    max_75    = FST_max - FST_P75,
    max_max3sd = FST_max - FST_max_3sdrm,
    max_maxmad = FST_max - FST_max_madrm) %>%
  ungroup()

#write.csv(sunrise_wide1,"Output_data/sunrise_wide_summary.csv", row.names = F)

################################################################################
# SUNRISE — NORMALITY AND DISTRIBUTION DIAGNOSTICS
################################################################################
# Apply the Shapiro–Wilk test to readings across sensors at each timestamp
# when at least three non-missing observations are available.
# A p-value below 0.05 indicates rejection of normality for that timestamp.
### Shapiro-wilk test
sunrise_normality <- sunrise_wide %>%
  rowwise() %>%
  mutate(n_sensor = sum(!is.na(c_across(all_of(sensor_cols)))),
         shapiro_p = {
           x <- c_across(all_of(sensor_cols))
           x <- x[!is.na(x)]
      
      if(length(x) >= 3) {
        shapiro.test(x)$p.value
      } else {
        NA_real_
      }},
    
    normal_005 = shapiro_p >= 0.05, 
    Hour = hour(bin15)
  ) %>%
  ungroup()

check_normal_hour <- sunrise_normality %>% filter(Hour %in% c(10:18),normal_005==FALSE)

#### During the daytime period (10:00–18:00), the Shapiro–Wilk test rejected normality for 373 of 504 15-minute sensor distributions (74.0%), indicating that the fruit-temperature distributions were frequently non-normal during the period of strongest solar heating.


# Generate histograms by clock time, pooling observations across dates
# to inspect the distribution of FST at each time of day.
sunrise_hist <- sunrise_QC_long %>% 
  dplyr::select(c("bin15","sensor","FST")) %>% na.omit() %>% 
  mutate(bin15 = ymd_hms(bin15),
         hms = format(bin15,"%H:%M:%S"),
         Hour = hour(bin15)) %>%  na.omit()

###### 
# Create the folder used for the distribution histograms, if needed.
if (!dir.exists("sunrise_dist/")) {
  dir.create("Figure/sunrise_dist", recursive = TRUE)
}

# create histogram for each timestamp across thermocouples/sensors
hms_count <- unique(sunrise_hist$hms)
for ( i in seq_along(hms_count)){
  data <- sunrise_hist %>% filter(hms == hms_count[i])
  
  P1<- ggplot(data, aes(x = FST)) +
    geom_histogram(
      bins = 12,
      boundary = 0
    ) +
    labs(
      x = "Fruit surface temperature (°C)",
      y = "Number of timesteps across thermocouples and days",
      title = paste0(hms_count[i])
    ) +
    theme_bw()
  
  clean_hms <- gsub(":", "-", as.character(hms_count[i]))
  
  ggsave(P1, file=paste0("Figure/sunrise_dist/Sunrise_",clean_hms,".png"), 
         width = 30, height = 20,dpi = 300, units = "cm")
}


# Calculate timestep-level cross-sensor skewness and classify its direction/magnitude.
sunrise_normality_1 <- sunrise_wide %>%
  rowwise() %>%
  mutate(
    skewness_FST = {
      x <- c_across(all_of(sensor_cols))
      x <- x[!is.na(x)]
      
      if(length(x) >= 3)
        e1071::skewness(x, type = 2)
      else
        NA_real_
    }
  ) %>%
  ungroup()

# add skewness columns
sunrise_normality_1 <- sunrise_normality_1 %>%
  mutate( bin15 = ymd_hms(bin15),
          Hour = hour(bin15),
    skew_direction = case_when(
      skewness_FST > 1 ~ "Strong positive skew",
      skewness_FST >= 0.5 ~ "Moderate positive skew",
      skewness_FST > -0.5 ~ "Approximately symmetric",
      skewness_FST > -1 ~ "Moderate negative skew",
      TRUE ~ "Strong negative skew"
    )
  ) %>% filter(Hour %in% c(10:18))

## About 58.5% of the daytime 15-minute thermocouple distributions show at least moderate skewness (\(|skewness| \ge 0.5\)), while only about 41.5% are approximately symmetric.

########During daytime hours (10:00–18:00), the Shapiro–Wilk test rejected normality for 74% of 15-minute timesteps. In addition, 58.5% of the timestep-specific thermocouple distributions exhibited moderate or strong skewness (\(|skewness| \ge 0.5\)), indicating that the distributions were frequently asymmetric as well as non-normal.

########################################
# QUINCY APPLE — SENSOR DATA PREPARATION AND QUALITY CONTROL

Quincy_combined<- read.csv("data/Quincy_combined_floor_15_timesteps.csv") %>% 
  rename("bin15"="Timesteps") %>% 
  mutate(
    # If a string lacks spaces/times (e.g. "2021-09-01"), append "00:00:00"
    bin15_clean = ifelse(nchar(bin15) == 10, paste0(bin15, " 00:00:00"), bin15),
    
    # Parse to datetime
    bin15 = ymd_hms(bin15_clean),
    Jday = yday(bin15)
  )

#### melt the dataset 
quincy_melt <- Quincy_combined %>% 
  reshape2::melt(id.vars = c('bin15','Jday'))

###### Quincy weather data 
Quincy_weather<- fread("data/Quincy_Sep2021_weather.csv",fill = TRUE, nrows = Inf) %>% 
  mutate(Timestamps = mdy_hm(paste(Date,Hour))) %>% 
  arrange(Timestamps) %>% 
  dplyr::select(c("Timestamps","Air_temperature_15","Dew_point","RH","Solar_Radiation","Wind_Speed_2m")) %>% 
  rename('bin15'= "Timestamps")

# Join observed FST with air temperature and compute FST–air differences.
quincy_melt_comb <- quincy_melt %>% 
  left_join(Quincy_weather %>% dplyr::select(c("bin15", "Air_temperature_15"))) %>% 
  mutate(value = as.numeric(value),
         Air_temp_C = (Air_temperature_15 - 32)*(5/9),
         diff = value - Air_temp_C,
         Hour = factor(hour(bin15),levels = 0:23))


# Plot the hourly distribution of Quincy thermocouple temperatures.
# Create Appendix Figure C01 (panel b)
P1<- ggplot(quincy_melt_comb, aes(x=Hour,y= value))+ 
  geom_boxplot()+ stat_boxplot(geom='errorbar')+  
  labs(x= "Hour (0:23)", y=  "Apple fruit surface temperature °C ", title = paste0("Quincy Apple"))+ theme_bw()+
  scale_fill_manual(values = c("apple" = "#D73027","blueberry" = "#4575B4"))+
  scale_x_discrete(drop = FALSE, breaks = c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22"))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 24,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 24,color = "black",angle = 0, hjust=1),
        legend.position = "bottom") 

# Save the Quincy hourly FST boxplot.
ggsave(P1, file=paste0("Figure/Quincy_FST_spread.png"), width = 30, height = 20,dpi = 300, units = "cm")


# Identify thermocouple columns for cross-sensor comparisons.
# Compute timestep-specific statistics and apply one-sided high-temperature screens.

themocouple_cols <- grep("^F", names(Quincy_combined), value = TRUE)
Quincy_QC_long <- Quincy_combined %>%
  pivot_longer(
    cols = all_of(themocouple_cols),
    names_to = "sensor",
    values_to = "FST"
  ) %>%
  group_by(bin15) %>%
  mutate(
    FST_mean = mean(FST, na.rm = TRUE),
    FST_median = median(FST, na.rm = TRUE),
    FST_sd   = sd(FST, na.rm = TRUE),
    
    upper_3SD = FST_mean + 3 * FST_sd,
    FST_max = max(FST, na.rm=T),
    FST_max_3sd = FST_max-upper_3SD,
    FST_MAD = mad(FST,constant = 1.4826, na.rm = TRUE),
    
    z_score = (FST - FST_mean) / FST_sd,
    robust_z =(FST - FST_median) / FST_MAD,
    flag_3SD = FST > upper_3SD,
    
    flag_MAD = robust_z > 3.5
  ) %>%
  ungroup()

# Extract readings flagged jointly (MAD and 3-SD) or by MAD alone for review.
quincy_qc <- Quincy_QC_long %>% filter(flag_MAD==TRUE & flag_3SD==TRUE)
quincy_qc_mad <- Quincy_QC_long %>% filter(flag_MAD==TRUE)
# Optional export of flagged records for manual inspection.
#write.csv(quincy_qc_mad, 'data/quincy_quality_issue_MAD.csv',row.names = F)
# Lists of thermocouples identified by the respective screening methods.
# Finally we used the MAD based approach
sensors_quincy_3sd <- c("F6","F18","F19")
sensors_quincy_mad <- c("F6","F9","F11","F16","F18","F22","F23")


# Calculate original and sensor-screened maximum FST at each timestamp.
Quincy_combined_1 <- Quincy_combined %>% 
  rowwise() %>%
  mutate(
    # --- Statistics using ALL F sensors ---
    FST_mean   = mean(c_across(all_of(themocouple_cols)), na.rm = TRUE),
    FST_median = median(c_across(all_of(themocouple_cols)), na.rm = TRUE),
    FST_sd     = sd(c_across(all_of(themocouple_cols)), na.rm = TRUE),
    FST_P75    = quantile(c_across(all_of(themocouple_cols)), probs = 0.75, na.rm = TRUE),
    FST_P85    = quantile(c_across(all_of(themocouple_cols)), probs = 0.85, na.rm = TRUE),
    FST_P95    = quantile(c_across(all_of(themocouple_cols)), probs = 0.95, na.rm = TRUE),
    FST_max    = max(c_across(all_of(themocouple_cols)), na.rm = TRUE),
    
    # --- Statistics EXCLUDING specific sensor vectors ---
    FST_max_3sdrm = max(c_across(all_of(setdiff(themocouple_cols, sensors_quincy_3sd))), na.rm = TRUE),
    FST_max_madrm = max(c_across(all_of(setdiff(themocouple_cols, sensors_quincy_mad))), na.rm = TRUE),
    
    # --- Difference calculations ---
    max_75     = FST_max - FST_P75,
    max_max3sd = FST_max - FST_max_3sdrm,
    max_maxmad = FST_max - FST_max_madrm
  ) %>%
  ungroup()

#write.csv(Quincy_combined_1, 'Output_data/Quincy_combined_1.csv',row.names = F)



######### 
################################################################################
# QUINCY — NORMALITY AND DISTRIBUTION DIAGNOSTICS
################################################################################
# Apply the Shapiro–Wilk test to readings across sensors at each timestamp
# when at least three non-missing observations are available.
# A p-value below 0.05 indicates rejection of normality for that timestamp.
quincy_normality <- Quincy_combined %>%
  rowwise() %>%
  mutate(n_sensor = sum(!is.na(c_across(all_of(themocouple_cols)))),
         shapiro_p = {
           x <- c_across(all_of(themocouple_cols))
           x <- x[!is.na(x)]
           
           if(length(x) >= 3) {
             shapiro.test(x)$p.value
           } else {
             NA_real_
           }},
         
         normal_005 = shapiro_p >= 0.05, 
         Hour = hour(bin15)
  ) %>%
  ungroup()

check_normal_hour_quincy <- quincy_normality %>% filter(normal_005==FALSE,Hour %in% c(10:18))
## 92/252 #287/708

#### During the daytime period (10:00–18:00), the Shapiro–Wilk test rejected normality for 92 of 252 15-minute sensor distributions (37.0%), indicating that the fruit-temperature distributions were frequently non-normal during the period of strongest solar heating.


# Create distributions by time of day, pooling measurements across dates.
quincy_hist <- Quincy_QC_long %>% 
  dplyr::select(c("bin15","sensor","FST")) %>% na.omit() %>% 
  mutate(hms = format(bin15,"%H:%M:%S"),
         Hour = hour(bin15)) %>%  na.omit()

# Create the output directory for Quincy distribution histograms.
if (!dir.exists("quincy_dist/")) {
  dir.create("Figure/quincy_dist", recursive = TRUE)
}

# Create histograms for normality check
hms_count <- unique(quincy_hist$hms)
for ( i in seq_along(hms_count)){
  data <- quincy_hist %>% filter(hms == hms_count[i])
  
  P1<- ggplot(data, aes(x = FST)) +
    geom_histogram(
      bins = 12,
      boundary = 0
    ) +
    labs(
      x = "Fruit surface temperature (°C)",
      y = "Number of timesteps across thermocouples and days",
      title = paste0(hms_count[i])
    ) +
    theme_bw()
  
  clean_hms <- gsub(":", "-", as.character(hms_count[i]))
  
  ggsave(P1, file=paste0("Figure/quincy_dist/quincy_",clean_hms,".png"), width = 30, height = 20,dpi = 300, units = "cm")
}


# Calculate timestep-level cross-sensor skewness and classify its direction/magnitude.

quincy_normality_1 <- Quincy_combined %>%
  rowwise() %>%
  mutate(
    skewness_FST = {
      x <- c_across(all_of(themocouple_cols))
      x <- x[!is.na(x)]
      
      if(length(x) >= 3)
        e1071::skewness(x, type = 2)
      else
        NA_real_
    }
  ) %>%
  ungroup()

# Add skewness class
quincy_normality_1 <- quincy_normality_1 %>%
  mutate( Hour = hour(bin15),
          skew_direction = case_when(
            skewness_FST > 1 ~ "Strong positive skew",
            skewness_FST >= 0.5 ~ "Moderate positive skew",
            skewness_FST > -0.5 ~ "Approximately symmetric",
            skewness_FST > -1 ~ "Moderate negative skew",
            TRUE ~ "Strong negative skew"
          )
  ) #%>% filter(Hour %in% c(10:18))
#130/252 #424/708


##mixed distributions; many timesteps are approximately symmetric/compatible with normality, but a substantial fraction are non-normal or skewed.

## Normality was evaluated across thermocouples independently at each 15-minute timestep using the Shapiro–Wilk test, supplemented by skewness and visual inspection of the distributions. At Sunrise, normality was rejected for 74% of daytime timesteps (10:00–18:00), and 58.5% of the daytime distributions exhibited moderate or strong skewness. At Quincy, normality was rejected for 36.5% of daytime timesteps, while 48.4% exhibited moderate or strong skewness. Thus, normality was not consistently supported across the two datasets, particularly during daytime solar-heating periods, with the strongest departures occurring at Sunrise.

### Because the normality assumption was not consistently satisfied and mean/standard deviation can themselves be influenced by extreme observations, a median- and MAD-based robust procedure was considered for identifying potentially anomalous thermocouple measurements



################################################################################
# CORVALLIS BLUEBERRY — THERMOCOUPLE QUALITY CONTROL
################################################################################
# Please obtain dataset first from Yang et al. (2021)
# The average of 16 bluberries (thermocouples) are used in Yang et al. (2021) model evaluation
# columns starts with "D0-E*" matches with their reported figure
# You can load whole dataset or just 16 columns we used in our analysis that matches with their figure
# I have used only 16 columns and timesteps columns
blueberry_all <- fread("data/blueberry_doe16_thermocouples.csv", fill = T, nrows = Inf)
blueberry_cols_doe <- grep("^D0-", names(blueberry_all), value = TRUE)

# Compute per-timestamp descriptive statistics and flag high-temperature outliers
# using a 3-SD upper threshold and a robust MAD-based score.
bluberry_long <- blueberry_all %>%
  pivot_longer(
    cols = all_of(blueberry_cols_doe),
    names_to = "sensor",
    values_to = "FST"
  ) %>%
  group_by(Time) %>%
  mutate(
    FST_mean = mean(FST, na.rm = TRUE),
    FST_median = median(FST, na.rm = TRUE),
    FST_sd   = sd(FST, na.rm = TRUE),
    
    upper_3SD = FST_mean + 3 * FST_sd,
    FST_max = max(FST, na.rm=T),
    FST_max_3sd = FST_max-upper_3SD,
    FST_MAD = mad(FST,constant = 1.4826, na.rm = TRUE),
    
    z_score = (FST - FST_mean) / FST_sd,
    robust_z =(FST - FST_median) / FST_MAD,
    flag_3SD = FST > upper_3SD,
    
    flag_MAD = robust_z > 3.5
  ) %>%
  ungroup()

# Extract MAD-flagged readings for review; this step does not remove them.
blueberry_qc <- bluberry_long %>% filter(flag_MAD==TRUE)
#write.csv(blueberry_qc,"Output_data/blueberry_qc_doe.csv", row.names = F)

# Thermocouple IDs selected for exclusion in the screened-maximum summary.
thermocouples_mad_doe <- c("D0-E-11","D0-E-12","D0-E-8","D0-E-6") #,"D0-E-4"

# Calculate original and sensor-screened maximum FST at each timestemp.
Blueberry_summary_1 <- blueberry_all %>% 
  rowwise() %>%
  mutate(
    # --- Statistics using ALL F sensors ---
    FST_mean   = mean(c_across(all_of(blueberry_cols_doe)), na.rm = TRUE),
    FST_median = median(c_across(all_of(blueberry_cols_doe)), na.rm = TRUE),
    FST_sd     = sd(c_across(all_of(blueberry_cols_doe)), na.rm = TRUE),
    FST_P75    = quantile(c_across(all_of(blueberry_cols_doe)), probs = 0.75, na.rm = TRUE),
    FST_P85    = quantile(c_across(all_of(blueberry_cols_doe)), probs = 0.85, na.rm = TRUE),
    FST_P95    = quantile(c_across(all_of(blueberry_cols_doe)), probs = 0.95, na.rm = TRUE),
    FST_max    = max(c_across(all_of(blueberry_cols_doe)), na.rm = TRUE),
    
    # --- Statistics EXCLUDING specific sensor vectors ---
    #FST_max_2sdrm = max(c_across(all_of(setdiff(blueberry_cols_doe, sensors_2sd_raspberry))), na.rm = TRUE),
    FST_max_madrm = max(c_across(all_of(setdiff(blueberry_cols_doe, thermocouples_mad_doe))), na.rm = TRUE),
    
    # --- Difference calculations ---
    max_75     = FST_max - FST_P75,
    max_maxmad = FST_max - FST_max_madrm) %>%
  ungroup()
#write.csv(Blueberry_summary_1,"Output_data/Blueberry_summary_1_doe.csv", row.names = F)

## This "Blueberry_summary_1" data set can be saved as "Tberry_full_observed.csv", this "FST_max" can be renamed as "Tmax_doE" or vice versa

# Didn't do any shapiro-wilk or skewness test here for blueberry dataset, as by now we decided to move with MAD based approach for exclusion of thermocouples and directly applied approach here for blueberry


##### blueberry box plots 
bluberry_long <- blueberry_all %>%
  pivot_longer(
    cols = all_of(blueberry_cols_doe),
    names_to = "sensor",
    values_to = "FST"
  ) %>% mutate(Time1 =  mdy_hm(Time),
               Hour = hour(Time1),
               Hour = factor(Hour, ordered = TRUE, labels = c(0:23))) %>% filter(FST <100)

# Compare the distribution of observed FST across hours of the day.
P1<- ggplot(bluberry_long, aes(x=Hour,y= FST))+ 
  geom_boxplot()+ 
  labs(x= "Hour (0:23)", y= "Blueberry fruit surface temperature °C ", title = paste0(" Corvallis Blueberry"))+ theme_bw()+
  scale_fill_manual(values = c("Baseline_scenario" = "#619CFF","OEC_scenario" = "#0000FF"))+
  scale_x_discrete(drop = FALSE, breaks = c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22"))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 24,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 24,color = "black",angle = 0, hjust=1),
        legend.position = "bottom",  
        legend.justification = c(0, 1))

ggsave(P1, file=paste0("Figure/blueberry_FST_spreads.png"), width = 30, height = 20,dpi = 300, units = "cm")




