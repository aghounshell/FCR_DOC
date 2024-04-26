###############################################################################

### Script to model/calculate 'whole-ecosystem' DOC processing (epi and hypo)
### Focusing on: 2017-2021 (5 years of whole-ecosystem manipulations)

### Last modified: 1 Mar. 2024, A. Hounshell
## Code review and updating EDI calls

### Last modified: 12 Apr. 2024, A. Hounshell
## Use bathy data from EDI
## Updates following DWH code-review (thanks, Dexter!!)

###############################################################################
# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,rMR,lme4,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics,truncnorm)

###############################################################################
## Load in various datastreams for box model
# Import Lake Analyzer thermocline results to determine median thermocline depth
# See: LakeAnalyzer_thermo.R - from 2016-2021
# Combining both CTD and YSI data to calculate thermocline depth from cast data
# collected at the catwalk
la_results <- read.csv("./Data/rev_FCR_results_LA.csv") %>% 
  select(datetime,thermo.depth) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = datetime)

### Load in Inflow data ----
# Weir discharge/temperature - https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=202&revision=8
# Last Downloaded: 01 March 23
# Using data from 2017-2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/8/cc045f9fe32501138d5f4e1e7f40d492"
#infile1 <- paste0(getwd(),"/Data/Inflow_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

inflow <- read.csv("./Data/Inflow_2013_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## If flows are 'below detection' (i.e., Flow_Flag = 3 or 13)
min_flow_cms <- min(inflow$WVWA_Flow_cms,na.rm=TRUE)

inflow <- inflow %>% 
  mutate(WVWA_Flow_cms = ifelse(WVWA_Flag_Flow %in% c(3,13),min_flow_cms,WVWA_Flow_cms),
         VT_Flow_cms = ifelse(VT_Flag_Flow %in% c(3,13),min_flow_cms,VT_Flow_cms))

# Use VT inflow to back-fill any missing WVWA inflow (where possible!)
# Plot WVWA vs. VT inflow - gut check!
ggplot(inflow,aes(x=WVWA_Flow_cms,y=VT_Flow_cms))+
  geom_point()

# Find relationship between WVWA and VT flow
# Use to fill in any missing WVWA values with VT flow (where possible)
flow_lm <- lm(WVWA_Flow_cms ~ VT_Flow_cms, data = inflow)

inflow <- inflow %>% 
  mutate(WVWA_Flow_cms = ifelse(is.na(WVWA_Flow_cms), flow_lm$coefficients[2]*VT_Flow_cms+flow_lm$coefficients[1], WVWA_Flow_cms)) %>% 
  mutate(flow_diff = abs(VT_Flow_cms - WVWA_Flow_cms))

# Average inflow by day and constrain to study period: 2017-2021
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_at(vars("WVWA_Flow_cms"),funs(mean(.,na.rm=TRUE),sd)) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & DateTime < as.POSIXct("2022-01-01"))

## Save daily inflow for ARIMA modeling
write.csv(inflow_daily, "./Data/inflow_daily.csv",row.names=FALSE)

# Create daily timeseries - missing data = NA
inflow_daily_full <- as.data.frame(seq(as.POSIXct("2017-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
inflow_daily_full <- inflow_daily_full %>% 
  dplyr::rename(DateTime = `seq(as.POSIXct("2017-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
inflow_daily_full <- left_join(inflow_daily_full, inflow_daily,by="DateTime")

# Calculate total variance - daily sd + difference between WVWA and VT inflow (0.002 cms)
# This will be used as inflow uncertainty for modeling
inflow_daily_full <- inflow_daily_full %>% 
  mutate(total_sd = sqrt((sd^2)+((mean(inflow$flow_diff,na.rm=TRUE))^2)))

# If total sd is NA for that day, then assume the SD is equal to the mean SD for the entire studyp period
inflow_daily_full <- inflow_daily_full %>% 
  mutate(total_sd = ifelse(is.na(total_sd),mean(inflow_daily_full$total_sd,na.rm=TRUE),sd))

# Create matrix of random variables from mean and total_sd
inflow_model_input <- as.data.frame(matrix(data = NA, ncol=1000,nrow=length(inflow_daily_full$DateTime)))

# Use rlnorm to constrain distribution to be >0
# See here: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/comment-page-1/
for (i in 1:length(inflow_daily_full$DateTime)){
  if(is.na(inflow_daily_full$mean[i])){
    inflow_model_input[i,1:1000] <- NA
  } else{
    # Find location and shape for rlnorm using the mean and sd
    m <- inflow_daily_full$mean[i]
    s <- inflow_daily_full$total_sd[i]/2
    location <- log(m^2 / sqrt(s^2 + m^2))
    shape <- sqrt(log(1 + (s^2 / m^2)))

    # Calculate distribution for each day
    inflow_model_input[i,1:1000] <- rlnorm(n=1000,mean=location,sd=shape)
  }
}

### Load DOC data ----
# From EDI: https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=199&revision=10
# Last downloaded: 1 Mar 2024
# Using data from 2017-2021
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/10/aa2ccc23688fc908f9d61cb217210a3d" 
#infile1 <- paste0(getwd(),"/Data/chemistry_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & DateTime < as.POSIXct("2022-01-01"))

# Separate into inflow
chem_100 <- chem %>% 
  filter(Site == 100) %>% 
  drop_na(DOC_mgL)

# And Site 50
chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  drop_na(DOC_mgL)

## Plot DOC concentrations from station 200, when available
## For SI
wetlands_time <- chem %>% 
  filter(Site == 200) %>% 
  drop_na(DOC_mgL) %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75)+
  geom_point(size=2)+
  xlab("")+
  ylab(expression(paste("Wetland DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

wetlands_box <- chem %>% 
  filter(Site == 200) %>% 
  drop_na(DOC_mgL) %>% 
  ggplot(mapping=aes(y=DOC_mgL))+
  geom_boxplot()+
  xlab("")+
  ylab(expression(paste("Wetland DOC (mg L"^-1*")")))+
  xlim(-1,1)+
  ylim(0,8)+
  theme_classic(base_size = 15)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggarrange(wetlands_time, ggarrange(wetlands_box, ncol = 2, labels = c("B."),font.label=list(face="plain",size=15)), 
          nrow = 2, labels = "A.", font.label=list(face="plain",size=15)) 

ggsave("./Figs/SI_Wetlands_DOC.jpg",width=7,height=7,units="in",dpi=320)

## Merge and plot measured inflow vs. measured DOC concentration
doc_100 <- left_join(chem_100,inflow_daily,by="DateTime") %>% 
  mutate(year = year(DateTime))

docvinflow <- doc_100 %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  ggplot(mapping=aes(x=log10(mean),y=DOC_mgL,color=as.factor(year)))+
  geom_point(size=1.5)+
  geom_smooth(method='lm')+ 
  scale_x_continuous(name=expression(paste("log"[10]*"(Inflow) (m"^3*"s"^-1*")")), breaks=c(-4,-3,-2,-1), label=c(0.0001,0.001,0.01,0.1))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  ylim(0,5)+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank())+
  theme(legend.position="top")

docvinflow_all <- doc_100 %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  ggplot(mapping=aes(x=log10(mean),y=DOC_mgL))+
  geom_point(size=1.5)+
  geom_smooth(method='lm')+ 
  scale_x_continuous(name=expression(paste("log"[10]*"(Inflow) (m"^3*"s"^-1*")")), breaks=c(-4,-3,-2,-1), label=c(0.0001,0.001,0.01,0.1))+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  ylim(0,5)+
  theme_classic(base_size=15)

combine <- ggarrange(docvinflow,docvinflow_all,nrow=1,ncol=2,labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

summary(doc_lm)
ggsave("./Figs/SI_DOCvInflow.png",combine,dpi=800,width=11,height=5)

doc_lm <- lm(DOC_mgL ~ mean, data = doc_100)

# Load in bathymetric data from EDI
# From EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1254.1
# Last downloaded: 12 Apr 2024
#i nUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184" 
# infile1 <- paste0(getwd(),"/Data/bathy_summary_stats.csv")
# download.file(inUrl1,infile1,method="curl")
bathy_data <- read.csv("./Data/bathy_summary_stats.csv", header=T) %>% 
  filter(Reservoir=="FCR")

vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0),
                         "Vol_m3" = c(sum(bathy_data$Volume_layer_L[1:4])/1000, #0-1 m (0.1 m)
                                     sum(bathy_data$Volume_layer_L[5:8])/1000, #1.3-2.3 m (1.6 m)
                                     sum(bathy_data$Volume_layer_L[9:14])/1000, #2.6-4.1 m (3.8 m)
                                     sum(bathy_data$Volume_layer_L[15:19])/1000, #4.4-5.6 m (5 m)
                                     sum(bathy_data$Volume_layer_L[20:23])/1000, #5.9-6.8 m (6.2 m)
                                     sum(bathy_data$Volume_layer_L[24:28])/1000, #7.1-8.3 m (8 m)
                                     sum(bathy_data$Volume_layer_L[29:31])/1000)) #8.7-9.3 m (9 m)

## Save for use in Env_ARIMA
write.csv(vol_depths, "./Data/vol_depths.csv",row.names=FALSE)

### Calculate DOC processing for Epi and Hypo -----
# First, separate by depth and convert to mass for each depth layer
# Remove time points if there are no observations in the epi, mid, or hypo layers
doc_box <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fill = NA, values_fn = mean, names_prefix = "DOC_") %>% 
  filter(!if_all(c(DOC_0.1, DOC_1.6), is.na)) %>% 
  filter(!if_all(c(DOC_3.8,DOC_5,DOC_6.2), is.na)) %>% 
  filter(!if_all(c(DOC_8,DOC_9),is.na))

# Create full timeseries
doc_box_full <- as.data.frame(seq(as.POSIXct("2017-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
doc_box_full <- doc_box_full %>% 
  dplyr::rename(DateTime = `seq(as.POSIXct("2017-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
doc_box_full <- left_join(doc_box_full,doc_box,by="DateTime")

doc_box_full <- doc_box_full %>% 
  mutate(DOC_0.1 = na.fill(na.approx(DOC_0.1,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_1.6 = na.fill(na.approx(DOC_1.6,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_3.8 = na.fill(na.approx(DOC_3.8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_5 = na.fill(na.approx(DOC_5,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_6.2 = na.fill(na.approx(DOC_6.2,na.rm=FALSE),"extend")) %>%   
  mutate(DOC_8 = na.fill(na.approx(DOC_8,na.rm=FALSE),"extend")) %>% 
  mutate(DOC_9 = na.fill(na.approx(DOC_9,na.rm=FALSE),"extend"))

## Create 'boot-strapped' values from [DOC] based on DOC MDL
# Create 3D array of random variables from the measured DOC and the MDL (0.11)
# For each timestep and each depth
doc_lake_model_input <- array(data = NA, dim=c(length(doc_box_full$DateTime),1000,7))

for (i in 1:length(doc_box_full$DateTime)){
  for (k in 1:7){
    # Find location and shape for rlnorm using the mean and sd
    m <- doc_box_full[i,k+1]
    s <- 0.11/2
    location <- log(m^2 / sqrt(s^2 + m^2))
    shape <- sqrt(log(1 + (s^2 / m^2)))
    
    doc_lake_model_input[i,1:1000,k] <- rlnorm(n=1000,mean=location,sd=shape)
  }
}

## Create boot-strapped parameters for volume - assuming volume +/-5%
vol_model_input <- as.data.frame(matrix(data=NA, ncol=1000,nrow=7))

for (i in 1:7){
  # Find location and shape for rlnorm using the mean and sd
  m <- vol_depths$Vol_m3[i]
  s <- vol_depths$Vol_m3[i]*0.05
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  
  # Calculate distribution for each depth
  vol_model_input[i,1:1000] <- rlnorm(n=1000,mean=location,sd=shape)
}

## Create boot-strapped parameters for DOC inflow - using MDL for the SD
doc_inflow_full <- as.data.frame(seq(as.POSIXct("2017-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
doc_inflow_full <- doc_inflow_full %>% 
  rename(DateTime = `seq(as.POSIXct("2017-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
doc_inflow_full <- left_join(doc_inflow_full,chem_100,by="DateTime")

doc_inflow_full <- doc_inflow_full %>% 
  select(DateTime,DOC_mgL) %>% 
  mutate(DOC_mgL = na.fill(na.approx(DOC_mgL,na.rm=FALSE),"extend"))

doc_inflow_input <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_inflow_full$DateTime)))

for (i in 1:length(doc_inflow_full$DateTime)){
  # Find location and shape for rlnorm using the mean and sd
  m <- doc_inflow_full$DOC_mgL[i]
  s <- 0.11/2
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  
  # Calculate distribution for each depth
  doc_inflow_input[i,1:1000] <- rlnorm(n=1000,mean=location,sd=shape)
}

###############################################################################
## Determine location of the epi and hypo for each time point
# Add in thermocline depth information
thermocline_depth <- as.data.frame(seq(as.POSIXct("2017-01-01",tz="EST"),as.POSIXct("2021-12-31",tz="EST"),by="days"))
thermocline_depth <- thermocline_depth %>% 
  dplyr::rename(DateTime = `seq(as.POSIXct("2017-01-01", tz = "EST"), as.POSIXct("2021-12-31", tz = "EST"), by = "days")`)
thermocline_depth <- left_join(thermocline_depth,la_results,by="DateTime")

thermocline_depth <- thermocline_depth %>% 
  select(DateTime,thermo.depth) %>% 
  mutate(thermo.depth = na.fill(na.approx(thermo.depth, na.rm=FALSE),"extend"))

# Use thermocline to find the location of the epi and hypo
thermocline_depth <- thermocline_depth %>% 
  mutate(thermo.depth = round(thermo.depth,digits=1)) %>% 
  mutate(epi_bottom_depth_m = ifelse(thermo.depth > 9.0, 9.0,
                                     ifelse(thermo.depth > 7.0, 8.0,
                                            ifelse(thermo.depth > 6.0, 6.2,
                                                   ifelse(thermo.depth > 4.4, 5.0,
                                                          ifelse(thermo.depth > 3.0, 3.8,
                                                                 ifelse(thermo.depth > 1.6, 1.6,
                                                                        ifelse(thermo.depth > 0.0, 0.1, NA)))))))) %>% 
  mutate(hypo_top_depth_m = ifelse(thermo.depth <= 0.0, 0.1,
                                   ifelse(thermo.depth <= 1.6, 1.6,
                                          ifelse(thermo.depth <= 3.0, 3.8,
                                                 ifelse(thermo.depth <= 4.4, 5.0,
                                                        ifelse(thermo.depth <= 6.0, 6.2,
                                                               ifelse(thermo.depth <= 7.0, 8.0,
                                                                      ifelse(thermo.depth <= 9.0, 9.0, NA))))))))

###############################################################################
### Calculate volume weighted epi and hypo concentration using variable thermocline from model above
doc_wgt <- left_join(doc_box,thermocline_depth,by="DateTime")

doc_wgt <- doc_wgt %>% 
  mutate(doc_epi_mgL = ifelse(epi_bottom_depth_m == 0.1, DOC_0.1*vol_depths$Vol_m3[1]/sum(vol_depths$Vol_m3[1]), 
                              ifelse(epi_bottom_depth_m == 1.6, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2]))/sum(vol_depths$Vol_m3[1:2]),
                                     ifelse(epi_bottom_depth_m == 3.8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                            ifelse(epi_bottom_depth_m == 5.0, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                   ifelse(epi_bottom_depth_m == 6.2, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                          ifelse(epi_bottom_depth_m == 8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA))))))) %>% 
  mutate(doc_hypo_mgL = ifelse(hypo_top_depth_m == 1.6, ((DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]), 
                               ifelse(hypo_top_depth_m == 3.8, ((DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                      ifelse(hypo_top_depth_m == 5, ((DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                             ifelse(hypo_top_depth_m == 6.2, ((DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                    ifelse(hypo_top_depth_m == 8, ((DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                           ifelse(hypo_top_depth_m == 9, (DOC_9*vol_depths$Vol_m3[7])/sum(vol_depths$Vol_m3[7]),NA)))))))

doc_wgt <- full_join(doc_wgt,chem_100, by = "DateTime")

doc_wgt <- doc_wgt %>% 
  select(DateTime, doc_epi_mgL, doc_hypo_mgL, DOC_mgL) %>% 
  arrange(DateTime) %>% 
  pivot_longer(!DateTime, names_to = "Loc", values_to = "DOC_mgL")

doc_wgt <- doc_wgt %>% 
  mutate(Loc = ifelse(Loc == "doc_epi_mgL", "Epi",
                      ifelse(Loc == "doc_hypo_mgL", "Hypo", 
                             ifelse(Loc== "DOC_mgL", "Inflow", NA)))) %>% 
  drop_na()

# Separate by year, too
doc_wgt <- doc_wgt %>% 
  mutate(year = year(DateTime))

## Save for ARIMA modeling
write.csv(doc_wgt, "./Data/EpiHypo_DOC.csv",row.names=FALSE)

## Plot vol weighted epi and hypo [DOC] and inflow [DOC]
vw_epi <- doc_wgt %>% 
  filter(Loc == "Epi") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#7EBDC2")+
  geom_point(size=2,color="#7EBDC2")+
  xlab("")+
  ylab(expression(paste("Epi VW DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

vw_hypo <- doc_wgt %>% 
  filter(Loc == "Hypo") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#393E41")+
  geom_point(size=2,color="#393E41")+
  xlab("")+
  ylab(expression(paste("Hypo VW DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

inflow_doc <- doc_wgt %>% 
  filter(Loc == "Inflow") %>% 
  ggplot(mapping=aes(x=DateTime,y=DOC_mgL))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75,color="#F0B670")+
  geom_point(size=2,color="#F0B670")+
  xlab("")+
  ylab(expression(paste("Inflow DOC (mg L"^-1*")")))+
  ylim(0,8)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

ggarrange(vw_epi,vw_hypo,inflow_doc,nrow=3,ncol=1,labels = c("A.", "B.", "C."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig3_VW_DOC.jpg",width=7,height=10,units="in",dpi=320)

all_box <- doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  ggplot(mapping=aes(x=Loc,y=DOC_mgL,fill=Loc))+
  geom_boxplot(size=0.8,alpha=0.5)+
  scale_fill_manual(breaks=c('Epi','Hypo','Inflow'),values=c("#7EBDC2","#393E41","#F0B670"))+
  xlab("")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  theme_classic(base_size = 15)+
  theme(legend.position = "none")

## Plot boxplots by year
year_box <- doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>%
  ggplot(mapping=aes(x=as.character(year),y=DOC_mgL,fill=Loc))+
  geom_boxplot(size=0.8,alpha=0.5)+
  scale_fill_manual(breaks=c('Epi','Hypo','Inflow'),values=c("#7EBDC2","#393E41","#F0B670"))+
  xlab("Year")+
  ylab(expression(paste("DOC (mg L"^-1*")")))+
  ylim(0,8)+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(year_box, ggarrange(all_box, ncol = 2, labels = c("B."),font.label=list(face="plain",size=15)), 
          nrow = 2, labels = "A.", font.label=list(face="plain",size=15)) 

ggsave("./Figs/SI_DOC_Boxplots.png",dpi=800,width=8,height=7)

## Calculate stats for SI:
doc_stats <- doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  group_by(Loc) %>% 
  summarise(min = min(DOC_mgL),
            max = max(DOC_mgL),
            median = median(DOC_mgL),
            mean = mean(DOC_mgL),
            sd = sd(DOC_mgL))

doc_stats_year <- doc_wgt %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  group_by(Loc, year) %>% 
  summarise(min = min(DOC_mgL),
            max = max(DOC_mgL),
            median = median(DOC_mgL),
            mean = mean(DOC_mgL),
            sd = sd(DOC_mgL))

write.csv(doc_stats_year,'./Figs/doc_stats_year.csv')

###############################################################################
## Have bootstrapped inputs (n=1000) for:
# In lake DOC concentrations for all depths (doc_lake_model_input)
# Inflow DOC concentrations (doc_inflow_input)
# Inflow rates (inflow_model_input)
# Volume of layers (vol_model_input)
# Then also have location of the thermocline (and therefore, of the epi and hypo)

## Then calculate processing

### Calculate total mass of DOC at each depth ###
doc_lake_mass <- array(data = NA, dim=c(length(doc_box_full$DateTime),1000,7)) # Mass of DOC at each depth (DOC_mgL * Vol_m3 = DOC_g)

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    for (k in 1:7){
      doc_lake_mass[i,j,k] <- vol_model_input[k,j]*doc_lake_model_input[i,j,k]
    }
  }
}

## Calculate mass of inflow - for hypo and epi
doc_inflow_mass <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_inflow_full$DateTime)))

for (j in 1:1000){
  for (i in 1:length(doc_inflow_full$DateTime)){
    doc_inflow_mass[i,j] <- as.numeric(doc_inflow_input[i,j])*as.numeric(inflow_model_input[i,j])*60*60*24
  }
}

### Calculate outflow from epi and hypo ###
# Epi outflow = inflow * [DOC] at 0.1 m = g/d
doc_epi_mass_outflow <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    doc_epi_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,1]*60*60*24
  }
}

# 'Outflow' to Epi from Hypo: concentration * outflow = mass/day
doc_hypo_mass_outflow <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Mass of DOC outflow from the hypo

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,2]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,3]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,4]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,5]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,6]*60*60*24
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9) {
      doc_hypo_mass_outflow[i,j] <- inflow_model_input[i,j]*doc_lake_model_input[i,j,7]*60*60*24
    }
  }
}

### Calculate [DOC]/dt for each time point ###

## First calculate mass of epi and hypo at any time point
# Include changes due to the thermocline
doc_epi_mass <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Mass of epi at each time point

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    if(thermocline_depth$epi_bottom_depth_m[i] == 0.1){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 1.6){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]+doc_lake_mass[i,j,2]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 3.8){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]+doc_lake_mass[i,j,2]+doc_lake_mass[i,j,3]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 5){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]+doc_lake_mass[i,j,2]+doc_lake_mass[i,j,3]+doc_lake_mass[i,j,4]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 6.2){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]+doc_lake_mass[i,j,2]+doc_lake_mass[i,j,3]+doc_lake_mass[i,j,4]+doc_lake_mass[i,j,5]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 8){
      doc_epi_mass[i,j] <- doc_lake_mass[i,j,1]+doc_lake_mass[i,j,2]+doc_lake_mass[i,j,3]+doc_lake_mass[i,j,4]+doc_lake_mass[i,j,5]+doc_lake_mass[i,j,6]
    }
  }
}

doc_hypo_mass <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Mass of hypo at each time point

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,2]+doc_lake_mass[i,j,3]+doc_lake_mass[i,j,4]+doc_lake_mass[i,j,5]+doc_lake_mass[i,j,6]+doc_lake_mass[i,j,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,3]+doc_lake_mass[i,j,4]+doc_lake_mass[i,j,5]+doc_lake_mass[i,j,6]+doc_lake_mass[i,j,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,4]+doc_lake_mass[i,j,5]+doc_lake_mass[i,j,6]+doc_lake_mass[i,j,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,5]+doc_lake_mass[i,j,6]+doc_lake_mass[i,j,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,6]+doc_lake_mass[i,j,7]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9){
      doc_hypo_mass[i,j] <- doc_lake_mass[i,j,7]
    }
  }
}

## Then calculate changing volume: 
epi_vol <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Volume of epi for each time point

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    if(thermocline_depth$epi_bottom_depth_m[i] == 0.1){
      epi_vol[i,j] <- vol_model_input[1,j]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 1.6){
      epi_vol[i,j] <- vol_model_input[1,j]+vol_model_input[2,j]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 3.8){
      epi_vol[i,j] <- vol_model_input[1,j]+vol_model_input[2,j]+vol_model_input[3,j]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 5){
      epi_vol[i,j] <- vol_model_input[1,j]+vol_model_input[2,j]+vol_model_input[3,j]+vol_model_input[4,j]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 6.2){
      epi_vol[i,j] <- vol_model_input[1,j]+vol_model_input[2,j]+vol_model_input[3,j]+vol_model_input[4,j]+vol_model_input[5,j]
    } else if(thermocline_depth$epi_bottom_depth_m[i] == 8){
      epi_vol[i,j] <- vol_model_input[1,j]+vol_model_input[2,j]+vol_model_input[3,j]+vol_model_input[4,j]+vol_model_input[5,j]+vol_model_input[6,j]
    }
  }
}

hypo_vol <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Volume of hypo for each time point

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    if(thermocline_depth$hypo_top_depth_m[i] == 1.6){
      hypo_vol[i,j] <- vol_model_input[2,j]+vol_model_input[3,j]+vol_model_input[4,j]+vol_model_input[5,j]+vol_model_input[6,j]+vol_model_input[7,j]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 3.8){
      hypo_vol[i,j] <- vol_model_input[3,j]+vol_model_input[4,j]+vol_model_input[5,j]+vol_model_input[6,j]+vol_model_input[7,j]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 5){
      hypo_vol[i,j] <- vol_model_input[4,j]+vol_model_input[5,j]+vol_model_input[6,j]+vol_model_input[7,j]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 6.2){
      hypo_vol[i,j] <- vol_model_input[5,j]+vol_model_input[6,j]+vol_model_input[7,j]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 8){
      hypo_vol[i,j] <- vol_model_input[6,j]+vol_model_input[7,j]
    } else if(thermocline_depth$hypo_top_depth_m[i] == 9){
      hypo_vol[i,j] <- vol_model_input[7,j]
    }
  }
}

## Now calculate DOC/dt for epi and hypo for each timepoint
doc_dt_epi <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Change in DOC/dt for the epi

for (j in 1:1000){
  for (i in 1:(length(doc_box_full$DateTime)-1)){
    doc_dt_epi[i+1,j] <- (doc_epi_mass[i+1,j]-doc_epi_mass[i,j])
  }
}

doc_dt_hypo <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Change in DOC/dt for the hypo

for (j in 1:1000){
  for (i in 1:(length(doc_box_full$DateTime)-1)){
    doc_dt_hypo[i+1,j] <- (doc_hypo_mass[i+1,j]-doc_hypo_mass[i,j])
  }
}

### Calculate entrainment from hypo to epi for each time point ###
# Loosely following FCR_DOCModel_edited_19May17 from Carey et al. 2018

# Entrainment
#if Entr is positive, then epi is getting bigger and hypo is getting smaller; 
#if Entr is negative, then hypo is getting bigger and epi is getting smaller

# Double check for calculating changes in thermocline depth
entr <- as.data.frame(matrix(data=NA, ncol=1, nrow=length(doc_box_full$DateTime)))
for (i in 2:length(doc_box_full$DateTime)){
  entr[i,1] <- thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[i-1]
}

doc_entr <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # Entrainment for each time point

for (i in 2:length(doc_box_full$DateTime)){
  if (entr$V1[i]==0){
    doc_entr[i,] = 0
  } else if (entr$V1[i]==1.2){
    doc_entr[i,] = doc_lake_mass[i,,5]
  }
}

### SOMETHING IS GOING WRONG HERE?? NOT 'CALCULATING' ENTRAINMENT TERM FOR CERTAIN DEPTHS

## ID entraiment based on changes in thermocline depth between timepoints
## Then determine the mass of DOC that moved between the epi and hypo if the thermocline depth changed

for (i in 2:length(doc_box_full$DateTime)){
  if(thermocline_depth$epi_bottom_depth_m[i] == thermocline_depth$epi_bottom_depth_m[(i-1)]){
    doc_entr[i,] <- 0
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -7.9){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,2]+doc_lake_mass[i,,3]+doc_lake_mass[i,,4]+doc_lake_mass[i,,5]+doc_lake_mass[i,,6])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -4.9){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,2]+doc_lake_mass[i,,3]+doc_lake_mass[i,,4])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -3.4){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,3]+doc_lake_mass[i,,4])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -3.0){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,5]+doc_lake_mass[i,,6])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -2.4){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,4]+doc_lake_mass[i,,5])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -2.2){
    doc_entr[i,] <- -doc_lake_mass[i,,3]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.8){
    doc_entr[i,] <- -doc_lake_mass[i,,6]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.5){
    doc_entr[i,] <- -doc_lake_mass[i,,2]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -1.2){
    doc_entr[i,] <- -doc_lake_mass[i,,5]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.2){
    doc_entr[i,] <- doc_lake_mass[i,,5]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.5){
    doc_entr[i,] <- sum(doc_lake_mass[i,,2])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 1.8){
    doc_entr[i,] <- doc_lake_mass[i,,6]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 2.2){
    doc_entr[i,] <- doc_lake_mass[i,,3]
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 3.0){
    doc_entr[i,] <- (doc_lake_mass[i,,5]+doc_lake_mass[i,,6])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 3.4){
    doc_entr[i,] <- (doc_lake_mass[i,,3]+doc_lake_mass[i,,4])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 6.4){
    doc_entr[i,] <- (doc_lake_mass[i,,3]+doc_lake_mass[i,,4]+doc_lake_mass[i,,5]+doc_lake_mass[i,,6])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == -3.7){
    doc_entr[i,] <- -1*(doc_lake_mass[i,,2]+doc_lake_mass[i,,3])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 3.7){
    doc_entr[i,] <- (doc_lake_mass[i,,2]+doc_lake_mass[i,,3])
  } else if(thermocline_depth$epi_bottom_depth_m[i]-thermocline_depth$epi_bottom_depth_m[(i-1)] == 6.1){
    doc_entr[i,] <- (doc_lake_mass[i,,2]+doc_lake_mass[i,,3]+doc_lake_mass[i,,4]+doc_lake_mass[i,,5])
  }
}


###############################################################################
## Then calculate DOC internal loading:

# Calculate processing for epi
doc_epi_process_g <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # DOC epi processing for each time point
doc_epi_process_mgL <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))

p = 0.74 # From Carey et al. 2018 - percentage of discharge to the epi vs. hypo

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    doc_epi_process_g[i,j] = doc_dt_epi[i,j]-(doc_inflow_mass[i,j]*p)-(doc_hypo_mass_outflow[i,j]*(1-p))+doc_epi_mass_outflow[i,j]-doc_entr[i,j]
  }
}

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    doc_epi_process_mgL[i,j] = (doc_dt_epi[i,j]-(doc_inflow_mass[i,j]*p)-(doc_hypo_mass_outflow[i,j]*(1-p))+doc_epi_mass_outflow[i,j]-doc_entr[i,j])/epi_vol[i,j]
  }
}

doc_hypo_process_g <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))
doc_hypo_process_mgL <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    doc_hypo_process_g[i,j] = doc_dt_hypo[i,j]-(doc_inflow_mass[i,j]*(1-p))+(doc_hypo_mass_outflow[i,j]*(1-p))+doc_entr[i,j]
  }
}

for (j in 1:1000){
  for (i in 1:length(doc_box_full$DateTime)){
    doc_hypo_process_mgL[i,j] = (doc_dt_hypo[i,j]-(doc_inflow_mass[i,j]*(1-p))+(doc_hypo_mass_outflow[i,j]*(1-p))+doc_entr[i,j])/hypo_vol[i,j]
  }
}

### Average across model runs and calculate sd for reach of the various inputs and for processing
doc_inputs_g <- as.data.frame(matrix(data=NA,nrow=length(doc_box_full$DateTime),ncol=20))

colnames(doc_inputs_g) <- c('mean_doc_inflow_g',
                            'sd_doc_inflow_g',
                            'mean_doc_hypo_outflow_g',
                            'sd_doc_hypo_outflow_g',
                            'mean_doc_dt_hypo_g',
                            'sd_doc_dt_hypo_g',
                            'mean_doc_entr_g',
                            'sd_doc_entr_g',
                            'mean_doc_dt_epi_g',
                            'sd_doc_dt_epi_g',
                            'mean_doc_epi_outflow_g',
                            'sd_doc_epi_outflow_g',
                            'mean_doc_epi_process_g',
                            'sd_doc_epi_process_g',
                            'mean_doc_epi_process_mgL',
                            'sd_doc_epi_process_mgL',
                            'mean_doc_hypo_process_g',
                            'sd_doc_hypo_process_g',
                            'mean_doc_hypo_process_mgL',
                            'sd_doc_hypo_process_mgL')

for (i in 1:length(doc_box_full$DateTime)){
  doc_inputs_g$mean_doc_inflow_g[i] = mean(as.numeric(doc_inflow_mass[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_inflow_g[i] = sd(as.numeric(doc_inflow_mass[i,]),na.rm=TRUE)
  
  doc_inputs_g$mean_doc_hypo_outflow_g[i] = mean(as.numeric(doc_hypo_mass_outflow[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_outflow_g[i] = sd(as.numeric(doc_hypo_mass_outflow[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_dt_hypo_g[i] = mean(as.numeric(doc_dt_hypo[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_dt_hypo_g[i] = sd(as.numeric(doc_dt_hypo[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_entr_g[i] = mean(as.numeric(doc_entr[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_entr_g[i] = sd(as.numeric(doc_entr[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_dt_epi_g[i] = mean(as.numeric(doc_dt_epi[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_dt_epi_g[i] = sd(as.numeric(doc_dt_epi[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_epi_outflow_g[i] = mean(as.numeric(doc_epi_mass_outflow[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_outflow_g[i] = sd(as.numeric(doc_epi_mass_outflow[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_epi_process_g[i] = mean(as.numeric(doc_epi_process_g[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_process_g[i] = sd(as.numeric(doc_epi_process_g[i,]),na.rm=TRUE)
  
  doc_inputs_g$mean_doc_epi_process_mgL[i] = mean(as.numeric(doc_epi_process_mgL[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_epi_process_mgL[i] = sd(as.numeric(doc_epi_process_mgL[i,]),na.rm=TRUE)

  doc_inputs_g$mean_doc_hypo_process_g[i] = mean(as.numeric(doc_hypo_process_g[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_process_g[i] = sd(as.numeric(doc_hypo_process_g[i,]),na.rm=TRUE)
  
  doc_inputs_g$mean_doc_hypo_process_mgL[i] = mean(as.numeric(doc_hypo_process_mgL[i,]),na.rm=TRUE)
  doc_inputs_g$sd_doc_hypo_process_mgL[i] = sd(as.numeric(doc_hypo_process_mgL[i,]),na.rm=TRUE)
}

DateTime <- doc_box_full %>% 
  select(DateTime)

final_doc_inputs_g <- cbind(DateTime,doc_inputs_g)

final_doc_inputs_g <- na.omit(final_doc_inputs_g)

final_doc_inputs_g <- final_doc_inputs_g %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"))

## Save final model output
write.csv(final_doc_inputs_g, "./Data/final_doc_inputs.csv",row.names=FALSE)

## Load in final model output (as needed!)
final_doc_inputs_g <- read.csv("./Data/final_doc_inputs.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Constrain to sampling days
doc_model_timepoints <- left_join(doc_box,final_doc_inputs_g,"DateTime")

###############################################################################
## Let's make some plots!
## Plot epi and hypo processing
doc_proc_17 <- ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_epi_process_g/1000-sd_doc_epi_process_g/1000,ymax=mean_doc_epi_process_g/1000+sd_doc_epi_process_g/1000,fill="Epi"),alpha=0.5)+
  geom_ribbon(final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_g/1000-sd_doc_hypo_process_g/1000,ymax=mean_doc_hypo_process_g/1000+sd_doc_hypo_process_g/1000,fill="Hypo"),alpha=0.5)+
  geom_line(final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_epi_process_g/1000,color="Epi"),size=1)+
  geom_line(final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_hypo_process_g/1000,color="Hypo"),size=1)+
  geom_point(doc_model_timepoints,mapping=aes(x=DateTime,y=mean_doc_epi_process_g/1000,color="Epi"),size=2.2)+
  geom_point(doc_model_timepoints,mapping=aes(x=DateTime,y=mean_doc_hypo_process_g/1000,color="Hypo"),size=2.2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-06-01"),as.POSIXct("2017-11-15"))+
  #ylim(-16,13)+
  ylab(expression(paste("DOC (kg d"^-1*")")))+
  xlab("2017")+
  theme_classic(base_size=15)+
  guides(fill="none")+
  theme(legend.title=element_blank())

doc_proc_18 <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=2.2)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=2.2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2018-06-01"),as.POSIXct("2018-11-15"))+
  ylim(-16,13)+
  ylab(expression(paste("DOC (mg L"^-1*" d"^-1*")")))+
  xlab("2018")+
  theme_classic(base_size=15)+
  guides(fill="none")+
  theme(legend.title=element_blank())

doc_proc_19 <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=2.2)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=2.2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2019-06-01"),as.POSIXct("2019-11-15"))+
  ylim(-16,13)+
  ylab(expression(paste("DOC (mg L"^-1*" d"^-1*")")))+
  xlab("2019")+
  theme_classic(base_size=15)+
  guides(fill="none")+
  theme(legend.title=element_blank())

doc_proc_20 <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=2.2)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=2.2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2020-06-01"),as.POSIXct("2020-11-15"))+
  ylim(-16,13)+
  ylab(expression(paste("DOC (mg L"^-1*" d"^-1*")")))+
  xlab("2020")+
  theme_classic(base_size=15)+
  guides(fill="none")+
  theme(legend.title=element_blank())

doc_proc_21 <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dotted",color="darkgrey")+
  geom_hline(yintercept=0,linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="Epi"),alpha=0.5)+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="Hypo"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=1)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="Epi"),size=2.2)+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="Hypo"),size=2.2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2021-06-01"),as.POSIXct("2021-11-15"))+
  ylim(-16,13)+
  ylab(expression(paste("DOC (mg L"^-1*" d"^-1*")")))+
  xlab("2021")+
  theme_classic(base_size=15)+
  guides(fill="none")+
  theme(legend.title=element_blank())

ggarrange(doc_proc_17,doc_proc_18,doc_proc_19,doc_proc_20,doc_proc_21,ncol=1,nrow=5,common.legend=TRUE,
          labels = c("A.", "B.", "C.", "D.","E."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig5_EpiHypo_DOC_Processing.jpg",width=10,height=14,units="in",dpi=320)

## Plot inputs/outputs to the Epi DOC model
epi_inflow <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=((mean_doc_inflow_g*0.74/1000)-(sd_doc_inflow_g*0.74/1000)),ymax=((mean_doc_inflow_g*0.74/1000)+(sd_doc_inflow_g*0.74/1000)),fill="Inflow"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.74/1000),color="Inflow"))+
  geom_point(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.74/1000),color="Inflow"))+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_hypo_outflow_g*0.26/1000)-(sd_doc_hypo_outflow_g*0.26/1000),ymax=(mean_doc_hypo_outflow_g*0.26/1000)+(sd_doc_hypo_outflow_g*0.26/1000),fill="Hypo Inflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Hypo Inflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Hypo Inflow"))+
  ylab(expression(paste("External DOC (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Inflow","Hypo Inflow"), values=c("#F0B670","#393E41"))+
  scale_fill_manual(breaks=c("Inflow","Hypo Inflow"),values=c("#F0B670","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank(),legend.position = "top")

epi_outflow <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_epi_outflow_g/1000)-(sd_doc_epi_outflow_g/1000),ymax=(mean_doc_epi_outflow_g/1000)+(sd_doc_epi_outflow_g/1000),fill="Epi Outflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_outflow_g/1000,color="Epi Outflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_outflow_g/1000,color="Epi Outflow"))+
  ylab(expression(paste("DOC Outflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Epi Outflow"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Epi Outflow"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.position = "none")

epi_change <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_dt_epi_g/1000-sd_doc_dt_epi_g/1000,ymax=mean_doc_dt_epi_g/1000+sd_doc_dt_epi_g/1000,fill="Epi DOC/dt"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_dt_epi_g/1000,color="Epi DOC/dt"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_dt_epi_g/1000,color="Epi DOC/dt"))+
  geom_ribbon(mapping=aes(x=DateTime,ymax=-(mean_doc_entr_g-sd_doc_entr_g)/1000,ymin=-(mean_doc_entr_g+sd_doc_entr_g)/1000,fill="Entr"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=-mean_doc_entr_g/1000,color="Entr"))+
  geom_point(mapping=aes(x=DateTime,y=-mean_doc_entr_g/1000,color="Entr"))+
  ylab(expression(paste("DOC/dt (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Epi DOC/dt","Entr"), values=c("#7EBDC2","#E7804B"))+
  scale_fill_manual(breaks=c("Epi DOC/dt","Entr"),values=c("#7EBDC2","#E7804B"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank(),legend.position = "top")

epi_internal <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_epi_process_g/1000-sd_doc_epi_process_g/1000,ymax=mean_doc_epi_process_g/1000+sd_doc_epi_process_g/1000,fill="Epi Internal"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_epi_process_g/1000,color="Epi Internal"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_epi_process_g/1000,color="Epi Internal"))+
  ylab(expression(paste("Internal DOC (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Epi Internal"), values=c("#7EBDC2"))+
  scale_fill_manual(breaks=c("Epi Internal"),values=c("#7EBDC2"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.position = "none")

ggarrange(epi_inflow,epi_outflow,epi_change,epi_internal,nrow=4,ncol=1,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/SI_Epi_model.jpg",width=9,height=12,units="in",dpi=320)

## Plot hypo model inputs/outputs
hypo_inflow <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=((mean_doc_inflow_g*0.26/1000)-(sd_doc_inflow_g*0.26/1000)),ymax=((mean_doc_inflow_g*0.26/1000)+(sd_doc_inflow_g*0.26/1000)),fill="Inflow"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.26/1000),color="Inflow"))+
  geom_point(mapping=aes(x=DateTime,y=(mean_doc_inflow_g*0.26/1000),color="Inflow"))+
  ylab(expression(paste("External DOC (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Inflow"), values=c("#F0B670"))+
  scale_fill_manual(breaks=c("Inflow"),values=c("#F0B670"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.position = "none")

hypo_outflow <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=(mean_doc_hypo_outflow_g*0.26/1000)-(sd_doc_hypo_outflow_g*0.26/1000),ymax=(mean_doc_hypo_outflow_g*0.26/1000)+(sd_doc_hypo_outflow_g*0.26/1000),fill="Hypo Outflow"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Hypo Outflow"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_outflow_g*0.26/1000,color="Hypo Outflow"))+
  ylab(expression(paste("DOC Outflow (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Hypo Outflow"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("Hypo Outflow"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.position = "none")

hypo_change <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymin=mean_doc_dt_hypo_g/1000-sd_doc_dt_hypo_g/1000,ymax=mean_doc_dt_hypo_g/1000+sd_doc_dt_hypo_g/1000,fill="Hypo DOC/dt"),alpha=0.50)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_dt_hypo_g/1000,color="Hypo DOC/dt"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_dt_hypo_g/1000,color="Hypo DOC/dt"))+
  geom_ribbon(mapping=aes(x=DateTime,ymax=(mean_doc_entr_g-sd_doc_entr_g)/1000,ymin=(mean_doc_entr_g+sd_doc_entr_g)/1000,fill="Entr"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_entr_g/1000,color="Entr"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_entr_g/1000,color="Entr"))+
  ylab(expression(paste("DOC/dt (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Hypo DOC/dt","Entr"), values=c("#393E41","#E7804B"))+
  scale_fill_manual(breaks=c("Hypo DOC/dt","Entr"),values=c("#393E41","#E7804B"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.title=element_blank(),legend.position = "top")

hypo_internal <- ggplot(final_doc_inputs_g)+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_ribbon(mapping=aes(x=DateTime,ymax=(mean_doc_hypo_process_g-sd_doc_hypo_process_g)/1000,ymin=(mean_doc_hypo_process_g+sd_doc_hypo_process_g)/1000,fill="Hypo Internal"),alpha=0.5)+
  geom_line(mapping=aes(x=DateTime,y=mean_doc_hypo_process_g/1000,color="Hypo Internal"))+
  geom_point(mapping=aes(x=DateTime,y=mean_doc_hypo_process_g/1000,color="Hypo Internal"))+
  ylab(expression(paste("Internal DOC (kg d"^-1*")")))+
  xlab("")+
  scale_color_manual(breaks=c("Hypo Internal"), values=c("#393E41"))+
  scale_fill_manual(breaks=c("Hypo Internal"),values=c("#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  guides(fill="none")+
  theme_classic(base_size=15)+
  theme(legend.position = "none")

ggarrange(hypo_inflow,hypo_outflow,hypo_change,hypo_internal,nrow=4,ncol=1,labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/SI_Hypo_Model.jpg",width=9,height=12,units="in",dpi=320)

###############################################################################
### Some numbers/calculations for the MS!
## Min, Max, Mean/SD for DOC internal sources
## Constrain to June-Turnover
## Epi:
final_doc_inputs_g %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(6,7,8,9,10)) %>% 
  summarise(min_epi = min(mean_doc_epi_process_mgL),
            max_epi = max(mean_doc_epi_process_mgL),
            med_epi = median(mean_doc_epi_process_mgL),
            mean_epi = mean(mean_doc_epi_process_mgL),
            sd_epi = sd(mean_doc_epi_process_mgL))

## Hypo:
final_doc_inputs_g %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(6,7,8,9,10)) %>% 
  summarise(min_hypo = min(mean_doc_hypo_process_mgL),
            max_hypo = max(mean_doc_hypo_process_mgL),
            med_hypo = median(mean_doc_hypo_process_mgL),
            mean_hypo = mean(mean_doc_hypo_process_mgL),
            sd_hypo = sd(mean_doc_hypo_process_mgL))

