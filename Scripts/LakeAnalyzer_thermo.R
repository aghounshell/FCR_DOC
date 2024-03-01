# Script to calculate thermocline depth with rLakeAnalyzer
# A Hounshell, 30 May 2021

# Updated on 12 Aug 2022 to include 2021 data

###############################################################################

# Clear workspace
rm(list = ls())

# Set working directory
wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,lubridate,zoo)

#devtools::install_github("GLEON/rLakeAnalyzer")

# CTD and YSI casts - combine together for most complete time-period
#need to import CTD observations from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/11/d771f5e9956304424c3bc0a39298a5ce" 
#infile1 <- paste0(getwd(),"/CTD_final_2013_2020.csv")
#download.file(inUrl1,infile1,method="curl")

ctd <- read.csv('./Data/CTD_2013_2021.csv') %>% #read in observed CTD data, which has multiple casts on the same day (problematic for comparison)
  filter(Reservoir=="FCR") %>%
  mutate(Date = as.POSIXct(strptime(Date, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:PAR_umolm2s)

ctd_50 <- ctd %>% 
  filter(Site==50) %>% 
  dplyr::rename(time = Date)

# Import YSI observations from EDI
#inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/198/8/07ba1430528e01041435afc4c65fbeb6"
#infile1 <- paste0(getwd(),"/YSI_PAR_profiles_2013-2020.csv")
#download.file(inUrl1,infile1,method="curl")

ysi <- read_csv('./Data/YSI_PAR_profiles_2013-2021.csv') %>% 
  filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d",tz='EST'))) %>% 
  select(Reservoir:pH)

ysi_50 <- ysi %>% 
  filter(Site==50) %>% 
  dplyr::rename(time = DateTime)

# Combine CTD and YSI data for site 50
# Select unique dates from both CTD and YSI casts
ysi_date_list <- as.data.frame(unique(as.Date(ysi_50$time)))
names(ysi_date_list)[1] <- "time"
ysi_date_list$ysi_fcr <- rep(-99,length(ysi_date_list$time))
ysi_date_list <- ysi_date_list %>% 
  arrange(time)

ctd_date_list <- as.data.frame(unique(as.Date(ctd_50$time)))
names(ctd_date_list)[1] <- "time"
ctd_date_list$ctd_fcr <- rep(-99,length(ctd_date_list$time))
ctd_date_list <- ctd_date_list %>% 
  arrange(time)

# Combine Unique dates list by date
fcr_dates <- merge(ysi_date_list, ctd_date_list, by="time", all.x=TRUE, all.y=TRUE)

### Merge data CTD and YSI datasets for FCR
fcr_merge <- merge(ctd_50, ysi_50, by="time", all.x=TRUE, all.y=TRUE)

# Find where there are Na values in the CTD data: need to do it for each column
ctd_fcr_na <- is.na(fcr_merge$Depth_m.x)
fcr_merge$Depth_m.x[ctd_fcr_na] <- fcr_merge$Depth_m.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$Temp_C.x)
fcr_merge$Temp_C.x[ctd_fcr_na] <- fcr_merge$Temp_C.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$DO_mgL.x)
fcr_merge$DO_mgL.x[ctd_fcr_na] <- fcr_merge$DO_mgL.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$DO_pSat)
fcr_merge$DO_pSat[ctd_fcr_na] <- fcr_merge$DOSat[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$Cond_uScm.x)
fcr_merge$Cond_uScm.x[ctd_fcr_na] <- fcr_merge$Cond_uScm.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$PAR_umolm2s.x)
fcr_merge$PAR_umolm2s.x[ctd_fcr_na] <- fcr_merge$PAR_umolm2s.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$ORP_mV.x)
fcr_merge$ORP_mV.x[ctd_fcr_na] <- fcr_merge$ORP_mV.y[ctd_fcr_na]

ctd_fcr_na <- is.na(fcr_merge$pH.x)
fcr_merge$pH.x[ctd_fcr_na] <- fcr_merge$pH.y[ctd_fcr_na]

fcr_all <- fcr_merge %>% 
  select(time,Depth_m.x,Temp_C.x,DO_mgL.x,DO_pSat,Cond_uScm.x,Chla_ugL,Turb_NTU,pH.x,ORP_mV.x,PAR_umolm2s.x) %>% 
  dplyr::rename(depth=Depth_m.x,Temp_C=Temp_C.x,DO_mgL=DO_mgL.x,Cond_uScm=Cond_uScm.x,pH=pH.x,ORP_mV=ORP_mV.x,PAR_umolm2s=PAR_umolm2s.x)

fcr_date_list <- as.data.frame(unique(as.Date(fcr_all$time)))

## Average across date and depth
fcr_all <- fcr_all %>% group_by(time,depth) %>% summarize_all(funs(mean),na.rm=TRUE) %>% arrange(time,depth)

## Save merged CTD and YSI data for later analysis
write.csv(fcr_all,"./Data/merged_YSI_CTD.csv", row.names = FALSE)

layer <- fcr_all %>% 
  dplyr::rename(Date = time, Depth_m = depth)

### Formatting for LakeAnalyzer ----
df.final<-data.frame()

layer1<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.1)))
layer1$Depth_m <- 0.1
layer2<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.5)))
layer2$Depth_m <- 0.5
layer3<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1)))
layer3$Depth_m <- 1.0
layer4<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.5)))
layer4$Depth_m <- 1.5
layer5<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2)))
layer5$Depth_m <- 2.0
layer6<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.5)))
layer6$Depth_m <- 2.5
layer7<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3)))
layer7$Depth_m <- 3.0
layer8<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.5)))
layer8$Depth_m <- 3.5
layer9<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4)))
layer9$Depth_m <- 4.0
layer10<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.5)))
layer10$Depth_m <- 4.5
layer11<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5)))
layer11$Depth_m <- 5.0
layer12<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.5)))
layer12$Depth_m <- 5.5
layer13<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6)))
layer13$Depth_m <- 6.0
layer14<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.5)))
layer14$Depth_m <- 6.5
layer15<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7)))
layer15$Depth_m <- 7.0
layer16<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.5)))
layer16$Depth_m <- 7.5
layer17<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8)))
layer17$Depth_m <- 8.0
layer18<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.5)))
layer18$Depth_m <- 8.5
layer19<-layer %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9)))
layer19$Depth_m <- 9.0

df.final = rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7,layer8,layer9,layer10,layer11,layer12,layer13,
                 layer14,layer15,layer16,layer17,layer18,layer19)

fcr_layers <- arrange(df.final, Date)
fcr_layers$Depth_m <- round(as.numeric(fcr_layers$Depth_m), digits = .5)

fcr_layers_temp <- fcr_layers %>% select(Date,Depth_m,Temp_C) %>% group_by(Date,Depth_m) %>% summarise_each(funs(mean))

fcr_layers_temp <- fcr_layers_temp %>% 
  mutate(Temp_C = na.fill(na.approx(Temp_C,na.rm=FALSE),"extend"))

fcr_new <- fcr_layers_temp %>% spread(Depth_m,Temp_C)

names(fcr_new)[1] <- "dateTime"
names(fcr_new)[2] <- "wtr_0.1"
names(fcr_new)[3] <- "wtr_0.5"
names(fcr_new)[4] <- "wtr_1.0"
names(fcr_new)[5] <- "wtr_1.5"
names(fcr_new)[6] <- "wtr_2.0"
names(fcr_new)[7] <- "wtr_2.5"
names(fcr_new)[8] <- "wtr_3.0"
names(fcr_new)[9] <- "wtr_3.5"
names(fcr_new)[10] <- "wtr_4.0"
names(fcr_new)[11] <- "wtr_4.5"
names(fcr_new)[12] <- "wtr_5.0"
names(fcr_new)[13] <- "wtr_5.5"
names(fcr_new)[14] <- "wtr_6.0"
names(fcr_new)[15] <- "wtr_6.5"
names(fcr_new)[16] <- "wtr_7.0"
names(fcr_new)[17] <- "wtr_7.5"
names(fcr_new)[18] <- "wtr_8.0"
names(fcr_new)[19] <- "wtr_8.5"
names(fcr_new)[20] <- "wtr_9.0"

fcr_new$dateTime <- ymd(fcr_new$dateTime)

fcr_new = data.frame(fcr_new)

fcr_new <- fcr_new %>% 
  arrange(dateTime)

# Export out for LakeAnalyzer
write.table(fcr_new, "./Data/fcr.wtr", sep='\t', row.names=FALSE)

###############################################################################

## Load in rLakeAnalyzer
pacman::p_load(rLakeAnalyzer)

## Move fcr.wtr to rLakeAnalyzer folder and setpath
wtr.path <- "C:/Users/ahoun/Documents/R/win-library/4.0/rLakeAnalyzer/extdata/fcr.wtr"

fcr_temp = load.ts(wtr.path)

thermo_depth = ts.thermo.depth(fcr_temp,seasonal=TRUE)

# Replace NaNs with NA
thermo_depth <- thermo_depth %>% 
  mutate(thermo.depth = ifelse(thermo.depth == "NaN", NA, thermo.depth))

# Plot to look at results
ggplot(thermo_depth,mapping=aes(x=datetime,y=-thermo.depth)) +
  geom_line()+
  xlim(as.POSIXct("2015-01-01"), as.POSIXct("2021-12-31"))

write.csv(thermo_depth,"./Data/rev_FCR_results_LA.csv")

###############################################################################