### Script to estimate uncertainty in discharge calculations
### Use pressure transducer error (0.1%) to estimate uncertainty
### in associated discharge calculations

### 1 Aug. 2025 - A. Hounshell
###################################################################
# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Create data folder within directory structure
if(file.exists("Data")==FALSE){
  dir.create(file.path(wd, "Data"))
}

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate,zoo,rLakeAnalyzer)
#####################################################################
### Load in Inflow data ----
# Weir discharge/temperature - https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=202&revision=8
# Last Downloaded: 01 March 23
# Using data from 2017-2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/8/cc045f9fe32501138d5f4e1e7f40d492"
#infile1 <- paste0(getwd(),"/Data/Inflow_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/Inflow_2013_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))%>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & DateTime < as.POSIXct("2022-01-01"))

# Create new dataframe with pressures only
diff_wvwa <- inflow %>%
    select(DateTime, WVWA_Pressure_psi, WVWA_Baro_pressure_psi, WVWA_Pressure_psia) %>%
    mutate(WVWA_Preassure_psi_low = WVWA_Pressure_psi - 0.001 * WVWA_Pressure_psi,
           WVWA_Preassure_psi_high = WVWA_Pressure_psi + 0.001 * WVWA_Pressure_psi,
           WVWA_Baro_pressure_psi_low = WVWA_Baro_pressure_psi - 0.001 * WVWA_Baro_pressure_psi,
           WVWA_Baro_pressure_psi_high = WVWA_Baro_pressure_psi + 0.001 * WVWA_Baro_pressure_psi) %>%
    mutate(WVWA_Pressure_psia_low = WVWA_Preassure_psi_low - WVWA_Baro_pressure_psi_low,
           WVWA_Pressure_psia_high = WVWA_Preassure_psi_high - WVWA_Baro_pressure_psi_high)

# Gut-check plot
ggplot(diff_wvwa, aes(x = DateTime)) +
  geom_line(aes(y = WVWA_Pressure_psia_low, color="Low"), size = 0.5) +
  geom_line(aes(y = WVWA_Pressure_psia_high, color="High"), size = 0.5) +
  geom_line(aes(y = WVWA_Pressure_psia, color="Actual"), size = 0.5) +
  labs(title = "WVWA Pressure with Uncertainty Bounds",
       x = "Date",
       y = "Pressure (psia)") +
  theme_minimal()

#####################################################################
## Calculate flow rates following: https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=202&revision=8
## Inflow_Aggregation_EDI_2021.R

### CALCULATE THE FLOW RATES AT INFLOW ### #(MEL 2018-07-06)
#####################################################################
#flow3 <- flow2*0.70324961490205 - 0.1603375 + 0.03048  # Response: Height above the weir (m); Distance between pressure sensor and weir lip: 0.1298575 m = 0.18337 psi
#pressure*conversion factor for head in m - distance from tip of transducer to lip of weir + distance from tip of transducer to pressure sensor (eq converted to meters)
#flow4 <- (0.62 * (2/3) * (1.1) * 4.43 * (flow3 ^ 1.5) * 35.3147) # Flow CFS - MEL: I have not changed this; should be rating curve with area of weir
#flow_final <- flow4*0.028316847                                  # Flow CMS - just a conversion factor from cfs to cms
#####################################################################

# BUT FIRST! We must separate out the data for different rating curves/time periods!
# separate the dataframe into pre and post v-notch weir to apply different equations
# Rectangular weir:
diff_pre <- diff_wvwa[diff_wvwa$DateTime< as.POSIXct('2019-06-07 01:00:00'),]

# Original installation of V-notch weir
diff_post <- diff_wvwa %>% 
  filter(DateTime >= as.POSIXct('2019-06-07 01:00:00') & DateTime <= as.POSIXct('2020-08-24 00:00:00'))

# V-notch weir post blow-out but before final pressure transducer location!
diff_aug20 <- diff_wvwa %>% 
  filter(DateTime > as.POSIXct('2020-08-24 00:00:00') & DateTime <= as.POSIXct('2020-09-02 15:00:00'))

# V-notch weir post blow-out with FINAL location of the pressure transducer
diff_sep20 <- diff_wvwa %>% 
  filter(DateTime > as.POSIXct('2020-09-02 15:00:00'))

#####################################################################
# the old weir equations are taken directly from MEL's Inflow Aggregation script
# Use for pressure data prior to 2019-06-06: see notes above for description of equations
# NOTE: Pressure_psia < 0.185 calculates -flows (and are automatically set to NA's)
# THIS WILL NOT CHANGE EVER AGAIN! KEEP AS IS FOR RECTANGULAR WEIR
diff_pre <- diff_pre %>% 
    mutate(flow1 = (WVWA_Pressure_psia)*0.70324961490205 - 0.1603375 + 0.03048,
            flow1_low = (WVWA_Pressure_psia_low)*0.70324961490205 - 0.1603375 + 0.03048,
            flow1_high = (WVWA_Pressure_psia_high)*0.70324961490205 - 0.1603375 + 0.03048) %>% 
    mutate(flow_cfs = (0.62 * (2/3) * (1.1) * 4.43 * (flow1 ^ 1.5) * 35.3147),
           flow_cfs_low = (0.62 * (2/3) * (1.1) * 4.43 * (flow1_low ^ 1.5) * 35.3147),
           flow_cfs_high = (0.62 * (2/3) * (1.1) * 4.43 * (flow1_high ^ 1.5) * 35.3147)) %>% 
    mutate(Flow_cms = flow_cfs*0.028316847,
           Flow_cms_low = flow_cfs_low*0.028316847,
           Flow_cms_high = flow_cfs_high*0.028316847) %>%
    select(DateTime, WVWA_Pressure_psia, WVWA_Pressure_psia_low, WVWA_Pressure_psia_high, 
            Flow_cms, Flow_cms_low, Flow_cms_high)

# Make flow NA when psi <= 0.184 (distance between pressure sensor and bottom of weir = 0.1298575 m = 0.18337 psi)
# Technically already completed above, but double check here
diff_pre$Flow_cms = ifelse(diff_pre$WVWA_Pressure_psia < 0.184, NA, diff_pre$Flow_cms)
diff_pre$Flow_cms_low = ifelse(diff_pre$WVWA_Pressure_psia_low < 0.184, NA, diff_pre$Flow_cms_low)
diff_pre$Flow_cms_high = ifelse(diff_pre$WVWA_Pressure_psia_high < 0.184, NA, diff_pre$Flow_cms_high)

#####################################################################
diff_post <- diff_post %>% 
    mutate(head = ((65.501*WVWA_Pressure_psia)-9.849)/100,
           head_low = ((65.501*WVWA_Pressure_psia_low)-9.849)/100,
           head_high = ((65.501*WVWA_Pressure_psia_high)-9.849)/100) %>% 
    mutate(Flow_cms = 2.391 * (head^2.5),
           Flow_cms_low = 2.391 * (head_low^2.5),
           Flow_cms_high = 2.391 * (head_high^2.5)) %>%
    select(DateTime, WVWA_Pressure_psia, WVWA_Pressure_psia_low, WVWA_Pressure_psia_high, 
            Flow_cms, Flow_cms_low, Flow_cms_high)

# According to gage height vs. pressure relationship, pressure < 0.010 should be flagged
diff_post$Flow_cms = ifelse(diff_post$WVWA_Pressure_psia <= 0.150, NA, diff_post$Flow_cms)
diff_post$Flow_cms_low = ifelse(diff_post$WVWA_Pressure_psia_low <= 0.150, NA, diff_post$Flow_cms_low)
diff_post$Flow_cms_high = ifelse(diff_post$WVWA_Pressure_psia_high <= 0.150, NA, diff_post$Flow_cms_high)

#####################################################################
diff_aug20 <- diff_aug20 %>% 
    mutate(head = ((55.556*WVWA_Pressure_psia)-1.8333)/100,
           head_low = ((55.556*WVWA_Pressure_psia_low)-1.8333)/100,
           head_high = ((55.556*WVWA_Pressure_psia_high)-1.8333)/100) %>% 
    mutate(Flow_cms = 2.391 * (head^2.5),
           Flow_cms_low = 2.391 * (head_low^2.5),
           Flow_cms_high = 2.391 * (head_high^2.5)) %>%
    select(DateTime, WVWA_Pressure_psia, WVWA_Pressure_psia_low, WVWA_Pressure_psia_high, 
            Flow_cms, Flow_cms_low, Flow_cms_high)

# According to rating curve, pressure < 0.033 should be removed
diff_aug20$Flow_cms = ifelse(diff_aug20$WVWA_Pressure_psia <= 0.033, NA, diff_aug20$Flow_cms)
diff_aug20$Flow_cms_low = ifelse(diff_aug20$WVWA_Pressure_psia_low <= 0.033, NA, diff_aug20$Flow_cms_low)
diff_aug20$Flow_cms_high = ifelse(diff_aug20$WVWA_Pressure_psia_high <= 0.033, NA, diff_aug20$Flow_cms_high)

#####################################################################
diff_sep20 <- diff_sep20 %>% 
    mutate(head = ((69.896*WVWA_Pressure_psia)-0.447)/100,
           head_low = ((69.896*WVWA_Pressure_psia_low)-0.447)/100,
           head_high = ((69.896*WVWA_Pressure_psia_high)-0.447)/100) %>%
  mutate(Flow_cms = 2.391 * (head^2.5),
         Flow_cms_low = 2.391 * (head_low^2.5),
         Flow_cms_high = 2.391 * (head_high^2.5)) %>%
    select(DateTime, WVWA_Pressure_psia, WVWA_Pressure_psia_low, WVWA_Pressure_psia_high, 
            Flow_cms, Flow_cms_low, Flow_cms_high)

# According to rating curve, pressure < 0.006 should be removed!
# This will need to be checked/updated each year as the rating curve changes
diff_sep20$Flow_cms = ifelse(diff_sep20$WVWA_Pressure_psia <= 0.006, NA, diff_sep20$Flow_cms)
diff_sep20$Flow_cms_low = ifelse(diff_sep20$WVWA_Pressure_psia_low <= 0.006, NA, diff_sep20$Flow_cms_low)
diff_sep20$Flow_cms_high = ifelse(diff_sep20$WVWA_Pressure_psia_high <= 0.006, NA, diff_sep20$Flow_cms_high)

# and put it all back together
diff <- rbind(diff_pre, diff_post, diff_aug20, diff_sep20)

## Plot as a gut-check
ggplot(diff, aes(x = DateTime)) +
  geom_line(aes(y = Flow_cms, color="Actual"), size = 0.5) +
  geom_line(aes(y = Flow_cms_low, color="Low"), size = 0.5, alpha=0.5) +
  geom_line(aes(y = Flow_cms_high, color="High"), size = 0.5, alpha=0.5) +
  labs(title = "WVWA Flow with Uncertainty Bounds",
       x = "Date",
       y = "Flow (cms)") +
  theme_minimal()

#####################################################################
## Estimate uncertainty from WVWA flow calculations - across all years
## for simplicity
diff %>%
    mutate(year = year(DateTime),
            doy = yday(DateTime)) %>%
    group_by(year, doy) %>%
    summarise(Flow_cms = mean(Flow_cms, na.rm = TRUE),
              Flow_cms_low = mean(Flow_cms_low, na.rm = TRUE),
              Flow_cms_high = mean(Flow_cms_high, na.rm = TRUE)) %>%
  mutate(Flow_cms_error = Flow_cms_high - Flow_cms_low,
         Flow_cms_error_pct = (Flow_cms_error / Flow_cms) * 100) %>%
    select(Flow_cms_error, Flow_cms_error_pct) %>%
    summarise(mean_error = mean(Flow_cms_error, na.rm = TRUE),
              mean_error_pct = mean(Flow_cms_error_pct, na.rm = TRUE))

## Daily mean error: ~5%

#####################################################################
## For completeness, check VT error, too!
## Includes an internal atmospheric pressure correction
## Assume same error as WVWA pressure transducer (0.1%)
diff_vt <- inflow %>%
    select(DateTime, VT_Pressure_psia) %>%
    mutate(VT_Pressure_psia_low = VT_Pressure_psia - 0.001 * VT_Pressure_psia,
           VT_Pressure_psia_high = VT_Pressure_psia + 0.001 * VT_Pressure_psia)

#####################################################################
# VT data for rectangular weir
# Used same equation as above following original Inflow Preparation script
# THIS WILL NEVER CHANGE! - WITH THE RECTANGULAR WEIR
VT_pre <- diff_vt[diff_vt$DateTime < as.POSIXct('2019-06-07 01:00:00'),]  
VT_pre <- VT_pre[complete.cases(VT_pre),]
VT_pre <- VT_pre %>% 
    mutate(head = (VT_Pressure_psia)*0.70324961490205 - 0.1603375 + 0.03048,
            head_low = (VT_Pressure_psia_low)*0.70324961490205 - 0.1603375 + 0.03048,
            head_high = (VT_Pressure_psia_high)*0.70324961490205 - 0.1603375 + 0.03048) %>% 
    mutate(flow_cfs = (0.62 * (2/3) * (1.1) * 4.43 * (head ^ 1.5) * 35.3147),
            flow_cfs_low = (0.62 * (2/3) * (1.1) * 4.43 * (head_low ^ 1.5) * 35.3147),
            flow_cfs_high = (0.62 * (2/3) * (1.1) * 4.43 * (head_high ^ 1.5) * 35.3147)) %>% 
    mutate(VT_Flow_cms = flow_cfs*0.028316847,
           VT_Flow_cms_low = flow_cfs_low*0.028316847,
           VT_Flow_cms_high = flow_cfs_high*0.028316847) %>% 
    select(DateTime, VT_Pressure_psia, VT_Pressure_psia_low, VT_Pressure_psia_high, 
            VT_Flow_cms, VT_Flow_cms_low, VT_Flow_cms_high) 

# Make flow as NA when psi <= 0.184 (distance between pressure sensor and bottom of weir = 0.1298575 m = 0.18337 psi)
# Technically already completed above, but double check here
VT_pre$VT_Flow_cms = ifelse(VT_pre$VT_Pressure_psia < 0.184, NA, VT_pre$VT_Flow_cms)
VT_pre$VT_Flow_cms_low = ifelse(VT_pre$VT_Pressure_psia_low < 0.184, NA, VT_pre$VT_Flow_cms_low)
VT_pre$VT_Flow_cms_high = ifelse(VT_pre$VT_Pressure_psia_high < 0.184, NA, VT_pre$VT_Flow_cms_high)

#####################################################################
# VT data for v-notch weir, pre blow-out (June 2019 - August 2020)
VT_post <- diff_vt %>% 
  filter(DateTime >= as.POSIXct('2019-06-07 01:00:00') & DateTime <= as.POSIXct('2020-08-24 15:15:00'))

VT_post <- VT_post %>% 
    mutate(head = ((70.64*VT_Pressure_psia)-5.663)/100,
           head_low = ((70.64*VT_Pressure_psia_low)-5.663)/100,
           head_high = ((70.64*VT_Pressure_psia_high)-5.663)/100) %>% 
    mutate(VT_Flow_cms = 2.391*(head^2.5),
           VT_Flow_cms_low = 2.391*(head_low^2.5),
           VT_Flow_cms_high = 2.391*(head_high^2.5)) %>% 
    select(DateTime, VT_Pressure_psia, VT_Pressure_psia_low, VT_Pressure_psia_high, 
            VT_Flow_cms, VT_Flow_cms_low, VT_Flow_cms_high) 

# Using rating curve, pressure at gage height = 0; pressure = 0.080
VT_post$VT_Flow_cms = ifelse(VT_post$VT_Pressure_psia <= 0.080, NA, VT_post$VT_Flow_cms)
VT_post$VT_Flow_cms_low = ifelse(VT_post$VT_Pressure_psia_low <= 0.080, NA, VT_post$VT_Flow_cms_low)
VT_post$VT_Flow_cms_high = ifelse(VT_post$VT_Pressure_psia_high <= 0.080, NA, VT_post$VT_Flow_cms_high)

#####################################################################
VT_aug20 <- diff_vt %>% 
  filter(DateTime >= as.POSIXct('2020-08-24 15:30:00') & DateTime <= as.POSIXct('2020-09-02 14:15:00'))

VT_aug20 <- VT_aug20 %>% 
    mutate(head = ((58.14*VT_Pressure_psia)+2.9302)/100,
            head_low = ((58.14*VT_Pressure_psia_low)+2.9302)/100,
            head_high = ((58.14*VT_Pressure_psia_high)+2.9302)/100) %>% 
    mutate(VT_Flow_cms = 2.391*(head^2.5),
           VT_Flow_cms_low = 2.391*(head_low^2.5),
           VT_Flow_cms_high = 2.391*(head_high^2.5)) %>% 
        select(DateTime, VT_Pressure_psia, VT_Pressure_psia_low, VT_Pressure_psia_high, 
            VT_Flow_cms, VT_Flow_cms_low, VT_Flow_cms_high) 

# Using rating curve, pressure at gage height = 0; pressure = -0.050
VT_aug20$VT_Flow_cms = ifelse(VT_aug20$VT_Pressure_psia <= 0, NA, VT_aug20$VT_Flow_cms)
VT_aug20$VT_Flow_cms_low = ifelse(VT_aug20$VT_Pressure_psia_low <= 0, NA, VT_aug20$VT_Flow_cms_low)
VT_aug20$VT_Flow_cms_high = ifelse(VT_aug20$VT_Pressure_psia_high <= 0, NA, VT_aug20$VT_Flow_cms_high)

#####################################################################
VT_sep20 <- diff_vt %>% filter(DateTime >= as.POSIXct('2020-09-02 14:30:00'))

VT_sep20 <- VT_sep20 %>% 
    mutate(head =((70.919*VT_Pressure_psia)+6.114)/100,
            head_low = ((70.919*VT_Pressure_psia_low)+6.114)/100,
            head_high = ((70.919*VT_Pressure_psia_high)+6.114)/100) %>% 
    mutate(VT_Flow_cms = 2.391*(head^2.5),
           VT_Flow_cms_low = 2.391*(head_low^2.5),
           VT_Flow_cms_high = 2.391*(head_high^2.5)) %>% 
    select(DateTime, VT_Pressure_psia, VT_Pressure_psia_low, VT_Pressure_psia_high, 
           VT_Flow_cms, VT_Flow_cms_low, VT_Flow_cms_high)

# Use rating curve, pressure at gage height = 0; pressure = -0.086
VT_sep20$VT_Flow_cms = ifelse(VT_sep20$VT_Pressure_psia <= 0, NA, VT_sep20$VT_Flow_cms)
VT_sep20$VT_Flow_cms_low = ifelse(VT_sep20$VT_Pressure_psia_low <= 0, NA, VT_sep20$VT_Flow_cms_low)
VT_sep20$VT_Flow_cms_high = ifelse(VT_sep20$VT_Pressure_psia_high <= 0, NA, VT_sep20$VT_Flow_cms_high)

VTinflow <- rbind(VT_pre, VT_post, VT_aug20, VT_sep20)

## Plot as a gut-check
ggplot(VTinflow, aes(x = DateTime)) +
  geom_line(aes(y = VT_Flow_cms, color="Actual"), size = 0.5) +
  geom_line(aes(y = VT_Flow_cms_low, color="Low"), size = 0.5, alpha=0.5) +
  geom_line(aes(y = VT_Flow_cms_high, color="High"), size = 0.5, alpha=0.5) +
  labs(title = "VT Flow with Uncertainty Bounds",
       x = "Date",
       y = "Flow (cms)") +
  theme_minimal()

#####################################################################
## Estimate uncertainty from VT flow calculations - across all years
## for simplicity
VTinflow %>%
  mutate(VT_Flow_cms_error = VT_Flow_cms_high - VT_Flow_cms_low,
         VT_Flow_cms_error_pct = (VT_Flow_cms_error / VT_Flow_cms) * 100) %>%
    select(VT_Flow_cms_error, VT_Flow_cms_error_pct) %>%
    summarise(mean_error = mean(VT_Flow_cms_error, na.rm = TRUE),
              mean_error_pct = mean(VT_Flow_cms_error_pct, na.rm = TRUE))

## Mean error: ~0.5%