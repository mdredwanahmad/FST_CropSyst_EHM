################################################################################
# Irrigation Demand Analysis — Apple and Blueberry
#
# Compare annual and monthly ET-only baseline irrigation with scenarios that
# include overhead evaporative cooling (OEC), using historical simulations and
# RCP 4.5 / RCP 8.5 climate projections for the 2040s and 2080s.
# The script also summarizes qualifying OEC application records and simulated
# crop phenological dates, and exports comparison figures.
#
# Irrigation depths are reported in mm; 'baseline' denotes ET-only demand,
# while the OEC scenario combines normal irrigation and applied cooling water.
################################################################################

# Load packages used for data manipulation and visualization.
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
library(scales, warn.conflicts = FALSE)
library(ggh4x, warn.conflicts = FALSE) 
library(patchwork, warn.conflicts = FALSE)


# HISTORICAL ANNUAL IRRIGATION DEMAND
# Load annual historical baseline (ET-only) irrigation demand and create a grid-year key.
baseline_hist <- fread("data/irrigation_demand/Baseline_historical_from_annual.csv",fill = TRUE, nrows = Inf) %>% 
  dplyr::select(c("Scenario","timeline","Fruit","Grid","Year","Irrigation","Precipitation")) %>% unique() %>%
  mutate(identifier = paste0(Fruit,"_",Grid,"_",Year))


# Load daily historical overhead evaporative cooling (OEC) simulations.
combined_hist_OEC_daily<- fread("data/irrigation_demand/combined_hist_OEC.csv",fill = TRUE, nrows = Inf) %>% 
  # Convert Year and Doy of the irrigation event to normal Date format to find out Months
  mutate(Date = as.Date(DOY-1, origin=paste0(Year,"-01-01")),
         Month = format(Date,"%m"))


# Aggregate historical daily irrigation components to annual totals by fruit and grid.
# Total applied irrigation includes normal irrigation plus OEC applications.
combined_hist_OEC_annual <- combined_hist_OEC_daily %>% 
  group_by(Scenario,timeline,Fruit,Grid,Year) %>% 
  summarise(total_irrig_normal = sum(Irrig_normal1),
            total_evaporative = sum(IrrigCooling_mm),
            Total_irrig_normal_col = sum(IrrigNormal),
            Total_applied_irr = total_irrig_normal+total_evaporative)

# Match historical ET-only and ET+OEC annual demand; omit the 1979 year.
# For historical timeframes we used 1980 - 2022
historical_comp <- left_join(baseline_hist,combined_hist_OEC_annual, by = c("timeline", "Fruit", "Grid", "Year")) %>% 
  filter(Year != 1979) %>% 
  rename("Baseline_scenario" = "Irrigation","OEC_scenario"="Total_applied_irr")

# FUTURE IRRIGATION DEMAND (2040s AND 2080s)
# FUTURE BASELINE IRRIGATION
# Read annual future ET-only baseline irrigation simulations.
# Analysis periods: 2025–2055 (2040s) and 2065–2095 (2080s).
future_year <- c(2025:2055, 2065:2095)
# Retain specified years and create a unique climate-model/grid/crop/year key.
baseline_future <- fread("data/irrigation_demand/Baseline_future_from_annual.csv", fill = TRUE, nrows = Inf) %>% 
  dplyr::select(c("Scenario","timeline","RCP","GCM","Grid","Fruit","Year","Irrigation")) %>% 
  filter(Year %in% future_year) %>% 
  mutate(time = ifelse(Year %in% c(2025:2056),"2040s","2080s"),
         RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         identifier = paste0(RCP,"_",GCM,"_",Grid,"_",Fruit,"_",Year))


# Calculate average baseline irrigation across GCMs within each RCP/grid/year.
baseline_future_avg <- baseline_future %>% 
  group_by(Scenario,timeline,RCP,Grid,Fruit,Year,time) %>% 
  summarise(Irrig_mean = mean(Irrigation))


# Load daily future ET+OEC simulations under climate projections.
combined_future_OEC<- fread("data/irrigation_demand/combined_future_OEC_sunburn_dat.csv", fill = TRUE,nrows = Inf) %>% 
  # Convert Year and Doy of the irrigation evant to normal Date format to find out Months
  mutate(Date = as.Date(DOY-1, origin=paste0(Year,"-01-01")),
         Month = format(Date,"%m"))

# Sum daily normal irrigation and OEC application depths to annual demand.
# Keep only future model/grid/crop/year combinations present in baseline data.
OEC_future <- combined_future_OEC %>% 
  group_by(Scenario,timeline,GCM,RCP,Fruit,Grid,Year) %>% 
  summarise(total_irrig_normal = sum(Irrig_normal1),
            total_evaporative = sum(IrrigCooling_mm),
            Total_irrig_normal_col = sum(IrrigNormal),
            Total_applied_irr = total_irrig_normal+total_evaporative) %>% 
  dplyr::select(c("Scenario","timeline","RCP","GCM","Grid","Fruit","Year","total_irrig_normal","total_evaporative","Total_applied_irr")) %>% 
  filter(Year %in% future_year) %>% 
  mutate(time = ifelse(Year %in% c(2025:2055),"2040s","2080s"),
         RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         identifier = paste0(RCP,"_",GCM,"_",Grid,"_",Fruit,"_",Year)) %>% 
  filter(identifier %in% baseline_future$identifier)



# Join baseline and OEC annual demand for matching climate and spatial units.
future_comp <- left_join(baseline_future,OEC_future, by = c("timeline","RCP","GCM","Grid","Fruit","Year","time","identifier")) %>%  
  rename("Baseline_scenario" = "Irrigation","OEC_scenario"="Total_applied_irr")

# Reshape future scenarios into long format for comparable boxplots.
future_comp_v1 <- future_comp %>% #filter(Month %in% c(5:9)) %>% 
  dplyr::select(c("Year","timeline","Fruit","RCP","GCM","Grid","time","Baseline_scenario","OEC_scenario")) %>% ungroup() %>% 
  reshape2::melt(id = c("Year","timeline","Fruit","RCP","GCM","Grid","time"))


# Reshape historical scenarios and assign historical labels for plotting.
historical_comp_v1 <- historical_comp %>% 
  dplyr::select(c("Year","timeline","Fruit","Grid","Baseline_scenario","OEC_scenario")) %>% ungroup() %>% 
  reshape2::melt(id = c("Year","timeline","Fruit","Grid")) %>% 
  mutate(RCP = "historical",
         GCM = "hist.",
         time = "hist.") %>% 
  dplyr::select(c("Year","timeline","Fruit","RCP","GCM","Grid","time","variable","value"))

# Combine historical and future periods with an explicit chronological order.
dataframe_list <- list(historical_comp_v1,future_comp_v1)
combined_dataframe <- do.call(rbind,dataframe_list) %>% 
  mutate(time = factor(time, ordered = TRUE, levels =c("hist.","2040s","2080s")),
         variable1 = paste0(RCP,"_",variable)) #%>% filter(Grid %in% check_grid)


#### Take median across GCM and Years
# combined_dataframe dataframe can be used to see the whole spread
median_gcm_year <- combined_dataframe %>% 
  group_by(timeline,Fruit,RCP,Grid,time,variable,variable1) %>% 
  summarise(median=  median(value, na.rm=T))

# optional, save the dataframe if needed
# write.csv(combined_dataframe,"data/irrigation_mm.csv", row.names = F)

# Specify consistent fill colors for baseline and OEC under each climate scenario.
scenario_comb_palette <- c(
  "historical_Baseline_scenario"       = "#B3D9FF", # Lightest Blue
  "historical_OEC_scenario"   = "#FFCC80", # Light Orange/Gold
  
  "RCP 4.5_Baseline_scenario"     = "#66B2FF", # Medium Blue
  "RCP 4.5_OEC_scenario" = "#FF8C00", # Pure Orange
  
  "RCP 8.5_Baseline_scenario"     = "#004C99", # Deep Navy
  "RCP 8.5_OEC_scenario" = "#CC6600"  # Darker Burnt Orange
)

median_gcm_year <- combined_dataframe %>% 
  group_by(timeline,Fruit,RCP,Grid,time,variable,variable1) %>% 
  summarise(median=  median(value, na.rm=T))

# ANNUAL DEMAND COMPARISON FIGURE
# Compare annual ET-only and ET+OEC irrigation demand by crop, period, and RCP.
# Create Figure 8.
P1<- ggplot(median_gcm_year, aes(x=time,y= median, fill = variable1))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom='errorbar')+ facet_wrap(~Fruit)+ 
  labs(x= " ", y= "Agricultural water demand (mm)", title = paste0(" "))+ theme_bw()+
  scale_y_continuous("Median agricultural water demand (mm)",limits=c(0,900), breaks=seq(0,900, by=100))+
  scale_fill_manual(values = scenario_comb_palette)+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 20), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 22,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 22,color = "black",angle = 0, hjust=1),
        legend.position = "bottom", 
        legend.justification = c(0, 1))

## Save figure
ggsave(P1, file=paste0("Figure/Ir_demand_comp_v1.png"), width = 30, height = 20,dpi = 300, units = "cm")



# Summarize median annual demand and express changes relative to historical
# baseline ET-only demand for each fruit.
combined_median <- median_gcm_year %>% 
  group_by(time,Fruit,Grid,RCP, variable) %>% 
  summarise(median_val = median(value, na.rm=T)) %>% 
  pivot_wider(names_from = variable, values_from = median_val) %>% 
  group_by(Fruit, Grid) %>% 
  mutate(baseline_hist = Baseline_scenario[time == "hist." & RCP == "historical"],
         OEC_hist = OEC_scenario[time == "hist." & RCP == "historical"]) %>% ungroup() %>% 
  mutate(ET_change_future = ((Baseline_scenario - baseline_hist)/baseline_hist)*100,
         OEC_percentage_in_future = ((OEC_scenario - Baseline_scenario)/Baseline_scenario)*100,
         percentage_increase_total = ((OEC_scenario - OEC_hist)/OEC_hist)*100,
         OEC_increase_comp_hist_base = ((OEC_scenario - baseline_hist)/baseline_hist)*100,
         Fruit = factor(Fruit,ordered = TRUE, levels =c("apple","blueberry")))
#write.csv(combined_median, "Output_data/percentage_demand_change_grid.csv",row.names=F)



# Tabulate median irrigation demand by crop, RCP, and irrigation scenario.
median_ir_data <- median_gcm_year %>% 
  #group_by(Fruit,RCP,Grid, variable) %>%
  group_by(time, Fruit,RCP, variable) %>%
  summarise( #irrigation_demand_mean = mean(value, na.rm = T),
    irrigation_demand_mead = median(value, na.rm = T)) %>% 
  pivot_wider(names_from = variable,
              values_from = irrigation_demand_mead)




# OEC EVENT FREQUENCY — HISTORICAL AND FUTURE
# This filter identifies records with exactly 8.80 mm of cooling applied.
# Count records meeting the specified OEC application criterion (8.80 mm).
hist_freq <- combined_hist_OEC_daily %>% 
  mutate(identifier = paste0(Fruit,"_",Grid,"_",Year)) %>% 
  filter(IrrigCooling_mm==8.80, Year != 1979, identifier %in% baseline_hist$identifier) %>% 
  dplyr::select(Scenario,timeline,Grid,Fruit,Date,Year,Month,DOY,FST_C,IrrigCooling_mm,EffectiveCooling_mm,NonEffective_mm,IrrigNormal)


# Count qualifying OEC records by month for each simulation unit.
hist_freq1 <- hist_freq %>% 
  group_by(Scenario,timeline,Grid,Fruit,Year,Month) %>% 
  summarise(frequncy_irrig_evap = n(),
            time = "historical")

# Select future records meeting the same 8.80-mm OEC criterion.
future_freq <- combined_future_OEC %>% 
  filter(IrrigCooling_mm==8.80, Year %in% future_year) %>% 
  dplyr::select(Scenario,timeline,RCP,GCM,Grid,Fruit,Date,Year,Month,DOY,FST_C,IrrigCooling_mm,EffectiveCooling_mm,NonEffective_mm,IrrigNormal)

# Count qualifying OEC records by month for each simulation unit.
future_freq1 <- future_freq %>% 
  group_by(Scenario,timeline,RCP,GCM,Grid,Fruit,Year,Month) %>% 
  summarise(frequncy_irrig_evap = n())

# Standardize RCP and future-period labels (GCM-specific records remain separate).
future_freq_2 <- future_freq1 %>% 
  mutate(RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         time = ifelse(Year %in% c(2025:2055),"2040s","2080s"))


# Duplicate historical records for comparison against both future periods.
hist_freq_40s <- hist_freq1 %>% 
  mutate(RCP = "historical",
         GCM = "hist.",
         time = "2040s") %>% 
  dplyr::select("Scenario","timeline","RCP","GCM","Grid","Fruit","Year","Month","frequncy_irrig_evap","time")
hist_freq_80s <- hist_freq1 %>% 
  mutate(RCP = "historical",
         GCM = "hist.",
         time = "2080s") %>% 
  dplyr::select("Scenario","timeline","RCP","GCM","Grid","Fruit","Year","Month","frequncy_irrig_evap","time")

# Combine historical and future summaries into a common plotting dataset.
list_freq <- list(hist_freq_40s,hist_freq_80s,future_freq_2)
combine_frequency <- do.call(rbind,list_freq) 

combine_frequency <- combine_frequency %>% 
  mutate(
    Month = factor(as.numeric(as.character(Month)),
                   levels = 3:10,
                   labels = as.character(3:10))
  )


# Visualize monthly counts of qualifying OEC records by crop and scenario.
# Create Appendix Figure C03
P1<- ggplot(combine_frequency, aes(x=Month,y= frequncy_irrig_evap, fill = Fruit))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom='errorbar')+  facet_grid(time~RCP)+ #outlier.shape = NA
  labs(x= " ", y= "Number of evaporative cooling irrigation events", title = paste0(" "))+ theme_bw()+
  scale_fill_manual(values = c("apple" = "#D73027","blueberry" = "#4575B4"))+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous("Number of evaporative cooling irrigation events",limits=c(0,25), breaks=seq(0,25, by=5))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 20), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 22,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 22,color = "black",angle = 0, hjust=1),
        legend.position = "bottom")

ggsave(P1, file=paste0("Figure/Irrig_freq_combine.png"), width = 30, height = 20,dpi = 300, units = "cm")



# ANNUAL OEC-FREQUENCY SUMMARY
# These counts are based on selected OEC application records; their interpretation
# as distinct heat-stress days depends on one qualifying record per day.
# Sum historical monthly OEC-event counts to annual values.
hist_freq_annual <- hist_freq1 %>% 
  group_by(Scenario,timeline,Grid,Fruit,Year) %>% 
  summarise(annual_freq = sum(frequncy_irrig_evap)) %>% 
  mutate(time = 'hist.',
         GCM = 'hist.',
         RCP = 'hist.') %>% 
  dplyr::select("Scenario","timeline","RCP","time","GCM","Grid","Fruit","Year","annual_freq")

# Sum future monthly OEC-event counts to annual values.
future_freq_2 <- future_freq1 %>% 
  mutate(RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         time = ifelse(Year %in% c(2025:2055),"2040s","2080s"))
future_freq_annual <- future_freq_2 %>%
  group_by(Scenario,timeline,RCP,time,GCM,Grid,Fruit,Year) %>% 
  summarise(annual_freq = sum(frequncy_irrig_evap))

# Combine historical and future annual counts for plottin
annual_freq <- rbind(hist_freq_annual, future_freq_annual) %>% 
  mutate(time = factor(time, ordered = TRUE, levels =c("hist.","2040s","2080s")))

# Plot annual totals of qualifying OEC records by crop and climate scenario.
P1 <- ggplot(annual_freq, aes(x=time,y= annual_freq, fill = RCP))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom='errorbar')+  facet_wrap(~Fruit)+ #outlier.shape = NA
  labs(x= " ", y= "Number of evaporative cooling irrigation events", title = paste0(" "))+ theme_bw()+
  #scale_fill_manual(values = c("apple" = "#D73027","blueberry" = "#4575B4"))+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous("Seasonal heat stress days",limits=c(0,50), breaks=seq(0,50, by=10))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 20,face="bold"), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 20,color = "black",angle = 0, hjust=0.5, vjust = 0.5,face="bold"),
        axis.text.y = element_text(size = 20,color = "black",angle = 0, hjust=1,face="bold"),
        legend.position = "bottom")

ggsave(P1, file=paste0("Figure/seasonal_heatsress_days.png"), width = 30, height = 20,dpi = 300, units = "cm")



# MONTHLY IRRIGATION DEMAND COMPARISON — JUNE THROUGH AUGUST
# Monthly total applied demand is normal irrigation plus overhead cooling.
hist_monthly <- combined_hist_OEC_daily %>% 
  mutate(identifier = paste0(Fruit,"_",Grid,"_",Year),
         RCP = "hist.",
         GCM = "hist.",
         Month = as.numeric(Month)) %>% 
  filter(Year != 1979, identifier %in% baseline_hist$identifier) %>% 
  dplyr::select(Scenario,timeline,Grid,RCP,GCM,Fruit,Date,Year,Month,DOY,IrrigCooling_mm,Irrig_normal1) %>% 
  group_by(Scenario,timeline,Fruit,RCP,GCM,Grid,Year,Month) %>% 
  summarise(total_evaporative = sum(IrrigCooling_mm),
            Total_irrig_normal_col = sum(Irrig_normal1),
            Total_applied_irr = (total_evaporative+Total_irrig_normal_col))

# Aggregate future OEC and normal irrigation to monthly totals.
# Select future records meeting the same 8.80-mm OEC criterion.
future_monthly <- combined_future_OEC %>% 
  filter(Year %in% future_year) %>% 
  dplyr::select(Scenario,timeline,RCP,GCM,Grid,Fruit,Date,Year,Month,DOY,FST_C,IrrigCooling_mm,Irrig_normal1) %>% 
  group_by(Scenario,timeline,Fruit,RCP,GCM,Grid,Year,Month) %>% 
  summarise(total_evaporative = sum(IrrigCooling_mm),
            Total_irrig_normal_col = sum(Irrig_normal1),
            Total_applied_irr = (total_evaporative+Total_irrig_normal_col)) %>% 
  mutate(time = ifelse(Year %in% c(2025:2055),"2040s","2080s"),
         RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         Month = as.numeric(Month))

# Prepare historical reference values alongside the future monthly summaries.
# Duplicate historical records for comparison against both future periods.
hist_40s <- hist_monthly %>% 
  mutate(time = "2040s")
hist_80s <- hist_monthly %>% 
  mutate(time = "2080s")

# Combine historical and future summaries into a common plotting dataset.
list_monthly <- list(hist_40s,hist_80s,future_monthly)
combine_monthly <- do.call(rbind,list_monthly) %>% filter(Month %in% c(6:8)) %>% 
  rename("Baseline scenario"='Total_irrig_normal_col', "OEC scenario" = 'Total_applied_irr') %>% 
  melt(id.vars = c("Scenario","timeline","Fruit","RCP","GCM","Grid","Year","Month","total_evaporative","time" )) %>% 
  mutate(variable1 = paste0(RCP,"_",variable),
         Month = factor(as.numeric(as.character(Month)),
                        levels = 6:8,
                        labels = as.character(6:8)))

# Select apple monthly demand records.
apple <- combine_monthly %>% filter(Fruit== "apple")

# Compare ET-only and ET+OEC monthly irrigation distributions.
## Create Figure Appendix C05 (panel A)

P1<- ggplot(apple, aes(x=Month,y= value, fill = variable))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom='errorbar')+  facet_grid(time~RCP)+ #outlier.shape = NA
  labs(x= " ", y= "Monthly OEC percentage of total demand", title = paste0("A) apple"))+ theme_bw()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous("Monthly agricultural water demand (mm)",limits=c(0,400), breaks=seq(0,400, by=100))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 24,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 24,color = "black",angle = 0, hjust=1),
        legend.position = "bottom")

ggsave(P1, file=paste0("Figure/Apple_monthly1.png"), width = 30, height = 20,dpi = 300, units = "cm")


# Select blueberry monthly demand records.
blueberry <- combine_monthly %>% filter(Fruit== "blueberry")

# Compare ET-only and ET+OEC monthly irrigation distributions.
## Create Figure Appendix C05 (panel B)
P1<- ggplot(blueberry, aes(x=Month,y= value, fill = variable))+ 
  geom_boxplot(outlier.shape = NA)+ stat_boxplot(geom='errorbar')+  facet_grid(time~RCP)+ #outlier.shape = NA
  labs(x= " ", y= "Monthly OEC percentage of total demand", title = paste0("B) blueberry "))+ theme_bw()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous("Monthly agricultural water demand (mm)",limits=c(0,400), breaks=seq(0,400, by=100))+
  theme(plot.title = element_text(size = 20),axis.title = element_text(size = 22), 
        legend.text = element_text(size=18), 
        strip.text = element_text(size=20), legend.title = element_blank(),
        axis.text.x = element_text(size = 24,color = "black",angle = 0, hjust=0.5, vjust = 0.5),
        axis.text.y = element_text(size = 24,color = "black",angle = 0, hjust=1),
        legend.position = "bottom")

ggsave(P1, file=paste0("Figure/blueberry_monthly1.png"), width = 30, height = 20,dpi = 300, units = "cm")


# PHENOLOGICAL RESPONSES — HISTORICAL AND FUTURE
# Read historical crop-stage dates and compute timing and duration metrics.
baseline_season <- fread("data/irrigation_demand/Baseline_historical_from_seasonal.csv", fill = T, nrows = Inf) %>% 
  mutate(GCM = "hist.",
         RCP = "hist.",
         Year = year(Emergence_date),
         DOY_emergence = yday(as.Date(Emergence_date)),
         DOY_flowering = yday(as.Date(Flowering_date)),
         DOY_Maturity = yday(as.Date(Maturity_date)),
         DOY_harvest = yday(as.Date(Harvest_date)),
         DOY_dormancy_start = yday(as.Date(Inactive_dormancy_begin_date)),
         season_length = DOY_harvest - DOY_emergence, 
         emergence_dormancy = DOY_dormancy_start - DOY_emergence,
         time = "hist.") %>% filter(Year != 1979)


# Load future crop-stage dates and calculate day-of-year and season duration.
future_season <- fread("data/irrigation_demand/Baseline_future_from_seasonal.csv", fill = T, nrows = Inf) %>% 
  mutate(Year = year(Emergence_date),
         DOY_emergence = yday(as.Date(Emergence_date)),
         DOY_flowering = yday(as.Date(Flowering_date)),
         DOY_Maturity = yday(as.Date(Maturity_date)),
         DOY_harvest = yday(as.Date(Harvest_date)),
         DOY_dormancy_start = yday(as.Date(Inactive_dormancy_begin_date)),
         season_length = DOY_harvest - DOY_emergence, 
         emergence_dormancy = DOY_dormancy_start - DOY_emergence,
         RCP = ifelse(RCP == "RCP45","RCP 4.5","RCP 8.5"),
         time = ifelse(Year %in% c(2025:2055),"2040s","2080s")) %>% filter(Year %in% future_year)


# Combine phenological records and normalize dates to a common display year
# so shifts in day/month can be compared across historical and future periods.
combine_seasons <- rbind(baseline_season,future_season) %>% 
  # Some simulated emergence dates/budbreak dates shift into the preceding December for warming, making
  # the calculated season length negative; omit these records for plotting
  # (397 of 42,130 records, as documented in the original analysis).
  filter(season_length>0) %>% 
  mutate(emergence_date = as.Date(DOY_emergence, origin = "1999-12-31"),
         harvest_date   = as.Date(DOY_harvest, origin = "1999-12-31"),
         maturity_date = as.Date(DOY_Maturity, origin = "1999-12-31"),
         dormancy_date = as.Date(DOY_dormancy_start, origin = "1999-12-31"),
         time = factor(time, ordered = TRUE, levels =c("hist.","2040s","2080s"))) %>% 
  dplyr::select(c("Year","Fruit","Grid","timeline","Scenario","GCM","RCP","time","season_length","emergence_date","harvest_date","dormancy_date","emergence_dormancy")) %>% #"date_maturity"
  reshape2::melt(id=c("Year","Fruit","Grid","timeline","Scenario","GCM","RCP","time","season_length","emergence_dormancy")) %>% 
  mutate(date_normalised = as.Date(paste0("2024-",month(value),"-",day(value))))

# APPLE PHENOLOGY — budbreak, harvest, and dormancy 
# create Appendix Figure C04 (panel A)
apple <- combine_seasons %>% filter(Fruit == "apple") %>% 
  mutate(variable = case_when(variable == "emergence_date" ~ "Budbreak date",
                              variable == "harvest_date" ~ "Harvest date",
                              variable == "dormancy_date" ~ "Dormancy date"),
         variable = factor(variable, ordered = TRUE, levels =c("Budbreak date","Harvest date","Dormancy date")))

P1<- ggplot(apple, aes(x = time, y = date_normalised, fill = RCP)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_boxplot(geom = 'errorbar') + 
  #facet_grid(variable ~ Fruit)+
  facet_wrap( ~ variable, scales = "free_y", ncol = 3) + 
  #facet_wrap(Fruit ~ variable, scales = "free_y", ncol = 3) + 
  theme_bw() + 
  # Use ggh4x to set specific scales for each column of facets
  facetted_pos_scales(y = list(
      # Emergence panels (Column 1)
      variable == "Budbreak date" ~ scale_y_date(
        limits = as.Date(c("2024-01-01", "2024-05-15")), 
        date_labels = "%b %d"),
      # Harvest panels (Column 2)
      variable == "Harvest date" ~ scale_y_date(
        limits = as.Date(c("2024-06-01", "2024-11-01")), 
        date_labels = "%b %d"),
      # Dormancy panels (Column 3)
      variable == "Dormancy date" ~ scale_y_date(
        limits = as.Date(c("2024-10-01", "2024-12-31")), 
        date_labels = "%b %d"))) +
  labs(x= " ", y= " ", title = paste0("A) apple "))+
  scale_fill_manual(values = c("hist." = "#F8766D", "RCP 4.5" = "#00BA38", "RCP 8.5" = "#619CFF")) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 22, face = "bold"), 
    legend.text = element_text(size = 18), 
    strip.text = element_text(size = 18), 
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.position = "bottom")

### Save plot
ggsave(P1, file=paste0("Figure/Phenology_apple.png"), width = 30, height = 20,dpi = 300, units = "cm")



# BLUEBERRY PHENOLOGY — budbreak, harvest, and dormancy 
# create Appendix Figure C04 (panel B)
blueberry <- combine_seasons %>% filter(Fruit == "blueberry") %>% 
  mutate(variable = case_when(variable == "emergence_date" ~ "Budbreak date",
                              variable == "harvest_date" ~ "Harvest date",
                              variable == "dormancy_date" ~ "Dormancy date"),
         variable = factor(variable, ordered = TRUE, levels =c("Budbreak date","Harvest date","Dormancy date")))

P1<- ggplot(blueberry, aes(x = time, y = date_normalised, fill = RCP)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_boxplot(geom = 'errorbar') + 
  #facet_grid(variable ~ Fruit)+
  facet_wrap( ~ variable, scales = "free_y", ncol = 3) + 
  #facet_wrap(Fruit ~ variable, scales = "free_y", ncol = 3) + 
  theme_bw() + 
  # Use ggh4x to set specific scales for each column of facets
  facetted_pos_scales(y = list(
    # Emergence panels (Column 1)
    variable == "Budbreak date" ~ scale_y_date(
      limits = as.Date(c("2024-01-01", "2024-05-15")), 
      date_labels = "%b %d"),
    # Harvest panels (Column 2)
    variable == "Harvest date" ~ scale_y_date(
      limits = as.Date(c("2024-06-01", "2024-11-01")), 
      date_labels = "%b %d"),
    # Dormancy panels (Column 3)
    variable == "Dormancy date" ~ scale_y_date(
      limits = as.Date(c("2024-10-01", "2024-12-31")), 
      date_labels = "%b %d"))) +
  labs(x= " ", y= " ", title = paste0("B) blueberry "))+
  scale_fill_manual(values = c("hist." = "#F8766D", "RCP 4.5" = "#00BA38", "RCP 8.5" = "#619CFF")) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18), 
    strip.text = element_text(size = 18), 
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.position = "bottom")

### Save plot
ggsave(P1, file=paste0("Figure/Phenology_blueberry.png"), width = 30, height = 20,dpi = 300, units = "cm")
