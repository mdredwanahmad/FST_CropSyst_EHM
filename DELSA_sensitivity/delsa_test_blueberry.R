library(tidyverse, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(viridisLite, warn.conflicts = FALSE)
library(viridis, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)


######### Energy balance model ######### 
unified_fruit_surface_temp <- function(Tair, RH, u_wind, Jday, SWRad, hour, LAT, LONG, elev,
                                       F_beam = 1, F_d = 1, F_a = 1, F_g = 0, F_r = 0,
                                       fruit_diameter = 0.1, fruit_emissivity = 0.97, fruit_albedo = 0.22, 
                                       r_fruit = 17958, method = "HS",day_id = NULL) {
  
  # coerce
  # Convert all input variables to numeric value
  to_num <- function(x) as.numeric(x)
  Tair  <- to_num(Tair);  RH <- to_num(RH);  u_wind <- pmax(to_num(u_wind), 0.01)
  Jday  <- to_num(Jday);  SWRad <- pmax(to_num(SWRad), 0)
  hour  <- to_num(hour);  LAT <- to_num(LAT); LONG <- to_num(LONG); elev <- to_num(elev)
  
  ## recycle to common length and make sure each variable length 1
  n <- length(Tair)
  recycle <- function(x) if (length(x) == n) x else rep(x[1], n)
  RH    <- recycle(RH); u_wind <- recycle(u_wind); Jday <- recycle(Jday)
  SWRad <- recycle(SWRad); hour <- recycle(hour); LAT <- recycle(LAT)
  LONG  <- recycle(LONG); elev <- recycle(elev)
  
  ## checking if data is finite, to force to make input data finite first
  # this will stop function if it finds input data is not finite
  # if (any(!is.finite(c(Tair, RH, u_wind, Jday, SWRad, hour, LAT, LONG, elev))))
  #   stop("Non-finite inputs detected.")
  
  # constants and global variables
  # atmospheric pressure based on elevation , eqn 3.7 (Campbell and Norman, 2000)
  exponent <- -(elev/8200)
  Patm <- 101.3*exp(exponent)
  R <- 287.058  ## gas constant, J kg-1 K-1
  lat  <- LAT  * pi/180  #latitude in radians
  # longitude (degrees west) negative
  long <- - as.numeric(LONG)        #absolute value of longitude for longitudinal correction
  stand_meridian <- 120 # standard meridian (degrees west for Pacific time)
  ## longitude corrections gives the time difference in hours between your location and standard meridians, 
  ## a negative value means west of the standard meridian and solar noon happens earlier than standard meridians and +(ve) means solar noon later
  long_correction <- (long - stand_meridian) / (360/24) #adjusting local solar time to match a standard meridian
  Linke <- 2.5 # Linke turbidity factor (Hofierka and Suri, 2002)
  sigma <- 5.67e-8
  earth_emissivity <- 0.95
  
  ## view factors for calculating incoming radiations
  # (Campbell and Norman, 2000; Evans, 2004)
  view_factor_beam <- F_beam
  view_factor_diff <- F_d
  view_factor_lwa  <- F_a
  view_factor_lwg  <- F_g
  view_factor_r    <- F_r
  
  ### Variables for computing fluxes
  es <- 0.611 * exp(17.3*Tair / (Tair + 237.3))   # saturated vapor pressure, kPa (FAO 56, eqn - 11)
  ea <- es * RH / 100                             # kPa
  lat_heat_vap  <- 2500 - 2.2274 * Tair           # kJ/kg
  thermal_cond_air <- 0.02624 * ((Tair + 273.15)/300)^0.8646
  air_spec_heat <- 1002.5 + 275e-6 * (Tair + 273.15 - 200)^2
  air_density   <- Patm * 1000 / pmax(R * (Tair + 273.15), 1e-9)
  
  dynamic_visc_air <- 1.458e-6 * (Tair + 273.15)^1.5 / (Tair + 273.15 + 110.4)
  kinetic_visc_air <- dynamic_visc_air / air_density
  Reynold_number   <- u_wind * fruit_diameter / pmax(kinetic_visc_air, 1e-12)
  Nusselt_number   <- 0.34 * pmax(Reynold_number, 1e-12)^0.6
  ## To cover the free flow condition when wind speed is 0, Nusselt number for free flow is 2 (Zhen et al 2024, Whitket 1972, P. mendez & Wikipedia)
  Nusselt_number_final <- pmax(Nusselt_number, 2.0) 
  thermal_diff_air <- thermal_cond_air / (air_density * air_spec_heat) #eqn 8.4 (Campbell and Norman, 2000)
  r_h <- (1/1.5) * fruit_diameter / pmax(thermal_diff_air * Nusselt_number_final, 1e-12)
  
  # vapor pressure conductance
  conductance_evap <- {
    bound_layer_thickness <- (2.8 * (fruit_diameter / u_wind)^0.5 + 0.25 / u_wind) / 1000
    r_bl <- bound_layer_thickness / pmax(thermal_diff_air, 1e-9)
    1 / pmax(r_fruit + r_bl, 1e-9)  # m/s
  }
  es_blueberry <- function(Tb) 0.611 * exp(17.3*Tb / (Tb + 237.3))  # Saturation vapor pressue of fruit surface temperature, kPa
  
  # ---- solar geometry (HS) -------------------------------------------------
  #################################
  ## Extraterrestrial radiation ###
  #################################
  
  ## Hofierka and Suri (2002) method
  day_angle <- 2*pi/365.25 * Jday #Hofierka and Suri, (2002)  equation 3
  time_adj  <- (279.575 + 0.986 * Jday) * pi/180 #Campbell and Norman, (2000) page 169
  #eqn_of_time #Campbell and Norman, (2000) equation 11.4
  eqn_of_time <- ((-104.7*sin(time_adj)+596.2*sin(2*time_adj)+4.3*sin(3*time_adj)-12.7*sin(4*time_adj)-
                     429.3*cos(time_adj)-2*cos(2*time_adj)+19.3*cos(3*time_adj))/3600)
  #eqn_of_time #Campbell and Norman, (2000) equation 11.3
  solar_time <- { h <- hour - 0.5; out <- h - long_correction + eqn_of_time; ifelse(out > 24, out - 24, out) } ## check
  hour_angle <- 0.261799 * (solar_time - 12) #Hofierka and Suri, (2002) 
  # Extraterrestrial irradiance normal to the solar beam, Hofierka and Suri, (2002)  eqn 1
  extra_rad_normal <- 1367 * (1 + 0.03344 * cos(day_angle - 0.048869))
  # Solar declination angle , Hofierka and Suri, (2002) Eqn 15
  sol_decl_angle   <- asin(0.3978 * sin(day_angle - 1.4 + 0.0355 * sin(day_angle - 0.0489)))
  # Solar elevation angle , Hofierka and Suri, (2002) Eqn 13
  sol_elev_angle   <- asin(cos(lat)*cos(sol_decl_angle)*cos(hour_angle) + sin(lat)*sin(sol_decl_angle))
  # Solar elevation angle , Hofierka and Suri, (2002) Eqn 33
  extra_rad <- ifelse(sin(sol_elev_angle) < 0, 0, extra_rad_normal * sin(sol_elev_angle))
  
  # ---- daily grouping key --------------------------------------------------
  # group it daily so that Julian day, Year, other variables doesn't confuse
  if (is.null(day_id)) day_id <- Jday
  day_id <- as.character(if (length(day_id) == n) day_id else rep(day_id[1], n))
  
  ##################################
  ## Direct and diffuse radiation ##
  ##################################
  
  ### calculation of direct rad 
  eps <- 1e-6
  rad_df <- data.frame(day_id, Jday = Jday, Rs = SWRad, Ra = extra_rad)
  # considering only day light hours when extra terrestrial radiation higher than at least 10
  # this physical guard prevents exploding KT and KB calculation for the night or evening time
  #rad_df <- subset(rad_df, Ra > 10) 
  d_Rs   <- stats::aggregate(Rs ~ day_id, data = rad_df, FUN = function(x) mean(x, na.rm=TRUE))
  d_Ra   <- stats::aggregate(Ra ~ day_id, data = rad_df, FUN = function(x) mean(x, na.rm=TRUE))
  d_J    <- stats::aggregate(Jday ~ day_id, data = rad_df, FUN = function(x) x[1])  # representative Jday
  
  daily <- Reduce(function(a,b) merge(a,b,by="day_id",all=TRUE), list(d_Rs,d_Ra,d_J))
  daily$K_T <- daily$Rs / pmax(daily$Ra, eps)
  daily$K_T <- daily$Rs / pmax(daily$Ra, eps)
  # As Kt is the ratio between observed or measured solar radiation (Rs) to extraterrestrial radiation (Ra). It should be within 0 to 1 (Allen et al 1998)
  # Perhaps clamping can be done till 1.1 to allow some noise from sensors or disaggregation errors
  #daily$K_T <- pmin(pmax(daily$K_T, 0), 1)
  # This deviation can be happened around morning or evening time when solar radiation has a small value based on the dataset accuracy
  #K_B = ratio of direct radation to extraterrestrial radiation (Vignola and McDaniels, 1986) #(Vignola and McDaniels, 1984)
  daily$K_B <- ifelse(
    daily$K_T <= 0.175,
    0.125 * daily$K_T^2,
    0.013 - 0.175*daily$K_T + 0.52*daily$K_T^2 + 1.03*daily$K_T^3 +
      (0.038*daily$K_T - 0.13*daily$K_T^2) * sin(2*pi/365 * (daily$Jday - 20))
  )
  daily$direct_rad <- daily$K_B * daily$Ra
  
  idx <- match(day_id, daily$day_id)
  direct_rad <- daily$direct_rad[idx] / pmax(daily$Rs[idx], eps) * SWRad
  ## Safeguard to prevent exploding direct radiation, if KT and KB calculation explodes for evening or night hours
  direct_rad <- pmin(pmax(direct_rad, 0), SWRad)
  ###
  
  ## Correct the radiation
  elev_angle <- sol_elev_angle
  direct_rad_corr <- ifelse(elev_angle > 0.36652, direct_rad / sin(elev_angle), direct_rad / sin(0.34907))
  ## Diffuse radiation
  diffuse_rad <- pmax(0, SWRad - direct_rad)   # clamp ≥ 0
  ## reflected radiation
  albedo <- 0.22 ## grass albedo
  reflect_sw <- view_factor_r * albedo * (diffuse_rad + direct_rad)
  ## Corrected Solar radiation
  SWRad_corr <- view_factor_beam * direct_rad_corr + view_factor_diff * diffuse_rad
  
  ## Longwave radation ##
  ### Ratio of real-to-clear sky global radiation ############
  # Corrected altitude angle, Hofierka and Suri, (2002), Eq. 6
  delta <- 0.061359*(0.1594+1.123*sol_elev_angle+0.065656*sol_elev_angle^2) /
    (1+28.9344*sol_elev_angle+277.3971*sol_elev_angle^2)
  sol_elev_angle_corr <- sol_elev_angle + delta
  
  # Air mass, Hofierka and Suri, (2002), Eq. 5
  air_mass <- exp(-elev / 8434.5) / (sin(sol_elev_angle_corr) + 0.50572*(sol_elev_angle_corr + 6.07995)^-1.6364)
  m <- air_mass
  # Rayleigh optical thickness, Hofierka and Suri, (2002), Eq. 8}
  optical_thickness <-ifelse(m<=20, 1/(6.6296+1.7513*m-0.1202*m^2+0.0065*m^3-0.00013*m^4), 1/(10.4+0.718*m)) # Rayleigh optical thickness, H.S. (2002), Eq. 8}
  rayleigh <- optical_thickness
  
  # clear-sky beam irradiance on a horizontal surface, Hofierka and Suri, (2002), Eqs. 4, 10
  clear_sky_beam <- pmax(0, extra_rad_normal * exp(-0.8662 * Linke * air_mass * rayleigh) * sin(sol_elev_angle))
  
  #transmission_function, Hofierka and Suri, (2002), Eqn. 22
  T_n <- (-0.015843 + 0.030543*Linke + 0.0003797*Linke^2)
  
  #sol_alt_function, Hofierka and Suri, (2002), Eqn. 23 & 24
  F_d <- {
    A1_p <- 0.26463 - 0.061581*Linke + 0.0031408*Linke^2
    A1 <- ifelse(A1_p*T_n < 0.0022, 0.0022/T_n, A1_p)
    A2 <- 2.04020 + 0.018945*Linke - 0.011161*Linke^2
    A3 <- -1.3025 + 0.039231*Linke + 0.0085079*Linke^2
    A1 + A2*sin(sol_elev_angle) + A3*sin(sol_elev_angle)^2
  }
  
  # Clear sky diffuse component on a horizontal surface, Hofierka and Suri, (2002), Eqn. 21
  clear_sky_diff <- pmax(0, extra_rad_normal * T_n * F_d)
  # clear-sky global irradiance on a horizontal surface, Hofierka and Suri, (2002), Eq. 36
  clear_sky_rad  <- clear_sky_beam + clear_sky_diff
  
  ## clear sky index
  cloud_df <- data.frame(day_id, Rs0 = clear_sky_rad, Rs = SWRad)
  dC <- stats::aggregate(cloud_df[, c("Rs0","Rs")], list(day_id = cloud_df$day_id),
                         function(x) mean(x, na.rm=TRUE))
  #csi <- pmin(1, dC$Rs / pmax(dC$Rs0, 1e-6))
  csi <- dC$Rs / pmax(dC$Rs0, 1e-6)
  clear_sky_index <- csi[ match(day_id, dC$day_id) ]
  clear_sky_index[!is.finite(clear_sky_index)] <- 0
  
  # Calculate longwave radiation
  # FAO-56 ((Allen et al., 1998)), Eq. 39, net outgoing long wave radiation
  # net_longwave_rad <- sigma * (0.34 - 0.14 * sqrt(pmax(ea, 0))) *
  #   (1.35 * clear_sky_index - 0.35) * (Tair + 273.16)^4
  cloud_term <- pmax(1.35 * clear_sky_index - 0.35, 0)
  net_longwave_rad <- sigma * (0.34 - 0.14 * sqrt(pmax(ea, 0))) * cloud_term * (Tair + 273.16)^4
  
  longwave_ground  <- earth_emissivity * sigma * (Tair + 273.15)^4
  longwave_atmosphere <- longwave_ground - net_longwave_rad
  longwave_absorbed  <- fruit_emissivity * (view_factor_lwa * longwave_atmosphere + view_factor_lwg * longwave_ground)
  
  ### Accounting energy fluxes
  ## Total incomimng solar radiation
  incoming_const <- (SWRad_corr + reflect_sw) * (1 - fruit_albedo) + longwave_absorbed
  
  ## Per element energy balance for solver to find out fruit surface temperature
  # latent heat
  latent_i <- function(i, Tb) - (lat_heat_vap[i] * 1000) * conductance_evap[i] * air_density[i] *
    (es_blueberry(Tb) - ea[i]) / Patm[i]
  # sensible heat
  sensible_i <- function(i, Tb) - air_spec_heat[i] * air_density[i] * (Tb - Tair[i]) / r_h[i]
  # energy radiattted from berry
  rad_i      <- function(   Tb) - fruit_emissivity * sigma * (Tb + 273.15)^4
  
  ## Energy balance for individual value
  EB_i <- function(i, Tb) latent_i(i, Tb) + sensible_i(i, Tb) + rad_i(Tb) + incoming_const[i]
  
  ## solve function
  solve_one <- function(i) {
    # coarse bracket to find the sign change within a pair and then find the root within the sign changed value
    grid <- c(-10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 65, 75)
    # find the values of EB for each coarse grid value
    vals <- sapply(grid, function(z) EB_i(i, z))
    # keep teh finite value to check sign change in the next step
    id <- which(is.finite(vals))
    # check if it has at least two valid value exist to find sign change
    if (length(id) >= 2) {
      # start loops through consecutive valid points
      for (k in seq_len(length(id)-1)) {
        # check if the function chnage sign between two neighbouring points
        if (sign(vals[id[k]]) != sign(vals[id[k+1]])) {
          #if it find a sign change in line 271 (vals) then it will iterate in line 280 to find FST value within the limit till EB=0
          # otherwise it go the next fall back loops and will try to find a value within -3 to 65 using same logic (line 285 to 288) and same uniroot function
          # if none works then it will fall back to the next bisection method to provide a solution
          # but it mostly solve in the first step in line #280 or line #288
          return(uniroot(function(z) EB_i(i, z), c(grid[id[k]], grid[id[k+1]]))$root)
        }
      }
    }
    # simple bounds check
    # if the model don't find sign change in previous step, then passing a practical range with more computationally laborious approach to find the root
    lo <- -5; hi <- 65
    f_lo <- EB_i(i, lo); f_hi <- EB_i(i, hi)
    # checking if the low and high value have finite EB value and have any sign change
    if (is.finite(f_lo) && is.finite(f_hi) && sign(f_lo) != sign(f_hi)) {
      uniroot(function(z) EB_i(i, z), c(lo, hi))$root
    } else {
      # minimize the absolute value of EB using optimize function to minimize the residual of EB
      # optimize a base r function based on brent root algorithm, its based on golden-section search and parabolic interpolation
      optimize(function(z) abs(EB_i(i, z)), c(lo, hi))$minimum
    }
  }
  
  ## Find optimized fruit surface temperature
  Tb1 <- vapply(seq_len(n), solve_one, numeric(1))
  ifelse(is.finite(Tb1), Tb1, NA_real_)
}



##########################
#### Load dataset  #######
#Delsa_data <- fread("D:/Extreme_weather/Paper/DELSA_SENSITIVITY/DELSA_LHS_with_local_sensitivity_ind_blueberry_crpparams.csv", fill = T, nrows = Inf)
#blueberry_combined_data<- fread("D:/Extreme_weather/Paper/DELSA_SENSITIVITY/bleuberry_combined_data.csv", fill = T, nrows = Inf)
#Delsa_data <- fread("/home/mdredwanahmad.khan/Extreme_weather_data/DELSA/DELSA_LHS_with_local_sensitivity_ind_blueberry_crpparams.csv", fill = T, nrows = Inf)
#blueberry_combined_data<- fread("/home/mdredwanahmad.khan/Extreme_weather_data/DELSA/bleuberry_combined_data.csv", fill = T, nrows = Inf)


# Wrapper for one DELSA parameter combination
run_fst_model <- function(weather, params) {
  
  # Make sure only one parameter combination is passed
  stopifnot(nrow(params) == 1)
  
  p <- params[1, ]
  
  weather %>%
    mutate(
      
      # -------------------------------------------------
      # Run identification
      # -------------------------------------------------
      run_id = p$run_id,
      sample_id = p$sample_id,
      
      # -------------------------------------------------
      # Perturb meteorological inputs
      # -------------------------------------------------
      
      # Temperature bias is in °C --> additive
      Air_temp_run_C = Tair + p$air_temp_bias_C,
      
      # RH bias 
      RH_run = RH + p$rh_bias_pct,
      
      # Keep RH physically between 0 and 100%
      RH_run = pmin(100, pmax(0, RH_run)),
      
      # Wind and radiation are expressed as multipliers
      Windspeed_run_ms = wind * p$wind_speed_multiplier,
      
      
      # Prevent negative values if ever generated
      Windspeed_run_ms = pmax(0, Windspeed_run_ms),
      
      
      # -------------------------------------------------
      # Run FST model
      # -------------------------------------------------
      Sim_FST = unified_fruit_surface_temp(
        
        Air_temp_run_C,
        RH_run,
        Windspeed_run_ms,
        Jday,
        SWRad,
        Hour,
        lat,
        long,
        elevation,
        
        # Radiation exposure/view factors
        F_d    = p$diffuse_sw_sky_view,
        F_a    = p$longwave_sky_view,
        F_r    = p$ground_reflected_sw_view,
        
        # Fruit characteristics
        fruit_diameter   = p$fruit_diameter_m,
        fruit_emissivity = p$fruit_emissivity,
        fruit_albedo     = p$fruit_albedo,
        r_fruit = p$canopy_resistance_s_m,
        
        day_id = day_id
      ),
      
      
      # DELSA bookkeeping
      base_sample_id = p$base_sample_id,
      run_type = p$run_type,
      perturbed_parameter = p$perturbed_parameter
    )
}



## FST model run with the DELSA sensitivity parameters 
blueberry_simulations <- purrr::map_dfr(
  seq_len(nrow(Delsa_data)),
  ~ run_fst_model(
    weather = blueberry_combined_data,
    params  = Delsa_data[.x, ] )) %>% mutate(Bias = Sim_FST-Observed_FST) %>% dplyr::select(c("bin15","run_id","Observed_FST","Sim_FST","Bias"))

# write the data back from DELSA analysis
write.csv(blueberry_simulations,"/home/mdredwanahmad.khan/Extreme_weather_data/DELSA/blueberry_simulations.csv", row.names=FALSE)