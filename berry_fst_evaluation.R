################################################################################
# Blueberry Fruit Surface Temperature (FST) Model Evaluation
# This script processes observed blueberry FST and meteorological data,
# simulates FST using the unified energy balance model, evaluates model
# performance, and generates comparison figures for all-day and peak-hour
# conditions. The processed dataset is also prepared for DELSA analysis.
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


################################################################################
# BLUEBERRY DATA PREPARATION — CORVALLIS, OREGON
################################################################################

# Blueberry FST and meteorological data are from Yang et al. (2021).
# The original dataset is available upon request from the corresponding author.
# Raw data are not hosted in this GitHub repository.

# Load meteorological data for July 16–August 22 and retain the variables
# required for FST simulation.

# Load observed blueberry FST data and derive time-related variables.
# FST_max_madrm represents the maximum observed FST after MAD-based filtering.

# Load meteorological data from July 16 to August 22 and retain the variables
# required for FST simulation.
blueberry_weather <- read.csv("data/model_inputs_full_blueberry.csv") %>% 
  mutate(Timestamp = mdy_hm(Time)) %>% 
  dplyr::select(c("Timestamp","Jday","Hour","SWRad","Tair","RH","wind"))

# Load observed blueberry FST data and derive time-related variables.
# FST_max_madrm represents maximum observed FST after MAD-based filtering.
blueberry_observed <- read.csv("data/Tberry_full_observed.csv") %>% 
  mutate(Timestamp = mdy_hm(Time),
         Jday = yday(Timestamp),
         Hour = hour(Timestamp)) %>% 
  dplyr::select(c("Timestamp","Jday","Hour","Tberry_max_all","Tmax_doE","FST_max_madrm"))


# Merge observed FST and meteorological data using matching timestamps.
# Aggregate measurements into 15-minute intervals and calculate mean values
# for the variables required for model simulation and evaluation.
Comb_blueberry_data <- left_join(blueberry_weather, blueberry_observed) %>%
  # Extract calendar variables from timestamps
  mutate(Year   = year(Timestamp), 
         Month  = month(Timestamp),                 
         Day    = mday(Timestamp)) %>% 
  group_by(Year, Jday) %>%
  mutate(Dewpoint = min(Tair, na.rm = TRUE),
         # round to the nearest boundary time using floor_date()
         bin15 = floor_date(Timestamp, "15 minutes")) %>% #60 minutes #"15 minutes"
  group_by(bin15) %>% 
  # Calculate interval averages and retain the observation count
  summarise(
    across(c("Year","Month","Day","Hour","Jday","SWRad","Tair","RH","wind","Tmax_doE","FST_max_madrm"), ~ mean(.x, na.rm = TRUE)), 
    n_obs = dplyr::n(),                                      
    .groups = "drop") %>% 
  mutate(Year_day = paste0(Year, "_", Jday))



################################################################################
# BLUEBERRY FST SIMULATION
################################################################################
Blueberry_FST <- Comb_blueberry_data %>% 
  mutate(FST_ind_berry_hot = unified_fruit_surface_temp(Tair, RH, wind,Jday,SWRad,Hour, 44.5646, -123.262, 80,F_beam = 1, 
                                                    F_d = 1, F_a = 1, F_g = 0, F_r = 0, fruit_diameter = 0.016, 
                                                    fruit_emissivity = 0.97, fruit_albedo = 0.22,r_fruit = 17958, day_id = Year_day))




### Create Figure 5 (panel b)
data1 <- Blueberry_FST %>% filter(Jday %in% c(204:206))
P1<- ggplot(data1)+
  geom_line( aes(x=bin15,y= FST_ind_berry_hot,color = "Simulated blueberry FST °C"), group=1)+
  geom_line( aes(x=bin15,y= FST_max_madrm,color = "Observed blueberry FST °C"), group=1)+
  scale_color_manual(values = c(
    "Observed blueberry FST °C" = "#00BA38",
    "Simulated blueberry FST °C" = "#0000FF", 'Simulated berry GEBM_Rnew' = "#619CFF", "Simulated blueberry ind. FST °C" = "#C71585",
    "FST_berry_hot_rdnew"= "brown","Observed_Tberry_max" = "#FFD700" , "Simulated_Yang et al" = "#F8766D"))+
  labs(x= "  ", y= "Blueberry fruit surface temperature °C", title = paste0(" "))+theme_bw()+ 
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 26), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 27,color = "black",angle = 0, hjust=0.5),
        axis.text.y = element_text(size = 27,color = "black",angle = 0, hjust=1),
        legend.position = "bottom", 
        legend.justification = c(0, 1),     
        legend.background = element_blank())  

ggsave(P1, file=paste0("Figure/Blueberry_FST.png"), width = 30, height = 22,dpi = 300, units = "cm")

################################################################################
# MODEL PERFORMANCE EVALUATION
################################################################################

# Calculate evaluation statistics after excluding non-finite observations.
# Errors are defined as simulated minus observed FST.
#
# RMSE: Root Mean Square Error
# MAE: Mean Absolute Error (stored as 'mab' in the function)
# MB:   Mean Bias (positive = overprediction; negative = underprediction)
# SD:   Standard deviation of simulation errors
# R²:   Coefficient of determination from linear regression
# r2_corr: Squared Pearson correlation coefficient
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

################################################################################
# ALL-DAY MODEL EVALUATION
################################################################################

# Calculate performance statistics using all available 15-minute observations.
Stats_berry <- compute_stats(Blueberry_FST$FST_max_madrm, Blueberry_FST$FST_berry_hot_ind)
Stats_berry

# Format model-performance statistics for display within the figure.
label_text_berry <- sprintf(
  "RMSE = %.2f °C\nMAB = %.2f °C\nMB = %.2f °C\nR²  = %.2f",
  Stats_berry$rmse, Stats_berry$mab, Stats_berry$mb, Stats_berry$r2_lm
)

### Create Figure 5 (panel d)
P1<- ggplot(Blueberry_FST, aes(x = FST_max_madrm, y = FST_berry_hot_ind))+
  geom_bin2d(bins = 70)+ #bins = 70
  scale_fill_continuous(type = "viridis")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size=1) + 
  scale_x_continuous("Observed fruit surface temp. °C", breaks=seq(0,60, by=10))+
  scale_y_continuous("Simulated fruit surface temp. °C ", breaks=seq(0,60, by=10))+ 
  coord_equal(xlim = c(0, 60), ylim = c(0, 60)) +
  annotate("label",
           x = 0.7, y = 60,     # adjust position as you like
           hjust = 0, vjust = 1,
           label = label_text_berry,
           size = 9, 
           linewidth = 0.25, 
           fill = "white")+ #"NA" for transparent
  labs(x = " ", y = " ",
       title = "Blueberry all day ")+theme_bw()+
  theme(
    axis.text = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 22, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.text = element_text(size = 20, color = "black"), 
    legend.title = element_text(size = 15, color = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.background = element_blank()
  )


ggsave(P1, file=paste0("Figure/blueberry_FST_all_day.png"), width = 30, height = 20,dpi = 300, units = "cm")


################################################################################
# PEAK-HOUR MODEL EVALUATION
################################################################################
# Restrict the evaluation to peak-temperature hours (12:00–18:00),
# including only the 18:00 observation from the final hour.
data <- Blueberry_FST %>%
  filter(Hour %in% 12:17 | (Hour == 18 & format(bin15, "%H:%M:%S") == "18:00:00"))
# Calculate model-performance statistics for peak-hour observations.
Stats_berry <- compute_stats(data$FST_max_madrm, data$FST_berry_hot_ind)
Stats_berry

# Format model-performance statistics for display within the figure.
label_text_berry <- sprintf(
  "RMSE = %.2f °C\nMAB = %.2f °C\nMB = %.2f °C\nR²  = %.2f",
  Stats_berry$rmse, Stats_berry$mab, Stats_berry$mb, Stats_berry$r2_lm
)

### Create Figure 5 (panel f)
P1<- ggplot(data, aes(x = FST_max_madrm, y = FST_berry_hot_ind))+
  geom_bin2d(bins = 70)+ #bins = 70
  scale_fill_continuous(type = "viridis")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size=1) + 
  scale_x_continuous("Observed fruit surface temp. °C", breaks=seq(0,60, by=10))+
  scale_y_continuous("Simulated fruit surface temp. °C ", breaks=seq(0,60, by=10))+ 
  coord_equal(xlim = c(0, 60), ylim = c(0, 60)) +
  annotate("label",
           x = 0.7, y = 60,     # adjust position as you like
           hjust = 0, vjust = 1,
           label = label_text_berry,
           size = 9, 
           linewidth = 0.25, 
           fill = "white")+ 
  labs(x = " ", y = " ",
       title = "Blueberry Peak hours (12:00 PM to 6:00 PM)")+theme_bw()+
  theme(
    axis.text = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 22, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 22),
    legend.text = element_text(size = 20, color = "black"), 
    legend.title = element_text(size = 15, color = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.background = element_blank()
  )


ggsave(P1, file=paste0("Figure/blueberry_FST_Peak_hours.png"), width = 30, height = 20,dpi = 300, units = "cm")


################################################################################
# BLUEBERRY DATA PREPARATION FOR DELSA SENSITIVITY ANALYSIS
################################################################################

# Prepare the combined blueberry dataset for subsequent DELSA analysis.
# Add site coordinates, elevation, fruit diameter, and station identifier.
# Standardize variable names for consistency with the apple evaluation dataset.
blueberry_dataset <- Blueberry_FST %>% 
  mutate(lat = 44.5646, 
         long= -123.262,
         elevation = 80,
         diameter = 0.016,
         Station = "Corvalis") %>% 
  rename("day_id" ="Year_day", "Observed_FST"='FST_max_madrm') %>% 
  dplyr::select(c( "bin15","Year","Hour","Jday","Observed_FST","Station","Tair","wind","RH","SWRad","lat","long","elevation","diameter","day_id"))
#write.csv(blueberry_dataset,"Output_data/bleuberry_combined_data.csv",row.names=F)
