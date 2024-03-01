### Script to look at environmental variables and conduct ARIMA modeling between 
### epi and hypo DOC concentrations and various limno. and met. parameters

### 18 Aug 2022, A. Hounshell

### Adding in various metrics of anoxia/oxygenation to hypo AR model
### 1 Sep 2022, A. Hounshell

### Updating following comments from CCC - DO mg/L -> DO%; remove GHG variables; ARIMA for DOC processing
### Adding in some additional visualizations/analyses for the Environmental data

###############################################################################

## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,zoo,scales,plyr,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

###############################################################################

## Load in thermocline data - obtained from CTD and YSI casts
## Calculated using LakeAnalyzer
thermo <- read.csv("./Data/rev_FCR_results_LA.csv") %>% 
  select(datetime,thermo.depth) %>% 
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = datetime)

## Estimate the location of the epi and hypo using thermocline depth
thermo <- thermo %>% 
  mutate(thermo.depth = round(thermo.depth,digits=1)) %>% 
  mutate(epi_bottomg_depth_m = ifelse(thermo.depth > 9.0, 9.0,
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

## Check average thermo depth from April-Oct from 2017-2021
thermo %>% 
  mutate(month = month(DateTime)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & month %in% c(4,5,6,7,8,9)) %>% 
  summarise_all(median,na.rm=TRUE)

## Plot timeseries of thermocline depth for SI
thermo %>% 
  drop_na(thermo.depth) %>% 
  ggplot(mapping=aes(x=DateTime,y=-thermo.depth))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = -0.1, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -1.6, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -3.8, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -5, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -6.2, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -8, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -9, linetype="dotted",color="darkgrey")+
  geom_line(size=1)+
  geom_point(size=2)+
  ylab("Thermocline depth (m)")+
  xlab("")+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

ggsave("./Fig_Output/SI_ThermoDepth.png",dpi=800,width=9,height=4)

###############################################################################

## Load in a format DOC concentration data from 2017-2021
chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

chem_50 <- chem %>% 
  filter(Site == 50) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5.0,6.2,8.0,9.0)) %>% 
  mutate(year = year(DateTime)) %>% 
  filter(DOC_mgL <= 15) %>% 
  drop_na(DOC_mgL)

## Calculate VW epi and hypo DOC concentrations

# Create vector of different volumes for each depth: based on L&O-L 2020 paper
vol_depths <- data.frame("Depth" = c(0.1,1.6,3.8,5.0,6.2,8.0,9.0), "Vol_m3" = c(138486.51,89053.28,59619.35,40197.90,13943.82,14038.52,1954.71))

# Select DOC concentrations and pivot longer
doc_mgL <- chem_50 %>% 
  select(DateTime,Depth_m,DOC_mgL) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = DOC_mgL, values_fil = NA, values_fn = mean, names_prefix = "DOC_")

doc_mgL <- left_join(doc_mgL, thermo, by="DateTime")

doc_mgL <- doc_mgL %>% 
  mutate(epi_DOC = ifelse(is.na(epi_bottomg_depth_m), ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(epi_bottomg_depth_m == 0.1, DOC_0.1,
                                 ifelse(epi_bottomg_depth_m == 1.6, (DOC_0.1*vol_depths$Vol_m3[1]+DOC_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                        ifelse(epi_bottomg_depth_m == 3.8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                               ifelse(epi_bottomg_depth_m == 5, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                      ifelse(epi_bottomg_depth_m == 6.2, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                             ifelse(epi_bottomg_depth_m == 8, ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DOC = ifelse(is.na(hypo_top_depth_m), ((DOC_0.1*vol_depths$Vol_m3[1])+(DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(hypo_top_depth_m == 1.6, ((DOC_1.6*vol_depths$Vol_m3[2])+(DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                  ifelse(hypo_top_depth_m == 3.8, ((DOC_3.8*vol_depths$Vol_m3[3])+(DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                         ifelse(hypo_top_depth_m == 5, ((DOC_5*vol_depths$Vol_m3[4])+(DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                ifelse(hypo_top_depth_m == 6.2, ((DOC_6.2*vol_depths$Vol_m3[5])+(DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                       ifelse(hypo_top_depth_m == 8, ((DOC_8*vol_depths$Vol_m3[6])+(DOC_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                              ifelse(hypo_top_depth_m == 9, DOC_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))
  
# Format for ARIMA modeling
final_doc_mgL <- doc_mgL %>% 
  select(DateTime,epi_DOC,hypo_DOC) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_DOC_mgL") %>% 
  mutate(Depth = ifelse(Depth == "epi_DOC", "Epi",
                        ifelse(Depth == "hypo_DOC", "Hypo", NA)))

## Calculate mean and SD for the hypo during the summer stratified period
final_doc_mgL %>% 
  mutate(month = month(DateTime)) %>% 
  filter(Depth == "Hypo" & month %in% c(6,7,8,9,10) & DateTime >= as.POSIXct("2017-01-01")) %>% 
  summarise(mean = mean(VW_DOC_mgL,na.rm=TRUE),
            sd = sd(VW_DOC_mgL,na.rm=TRUE))

###############################################################################
## Add in DOC processing - calculated from Eco_DOC.R
doc_processing <- read_csv("./Fig_Output/model_results.csv") %>% 
  select(DateTime, mean_doc_epi_process_mgL, mean_doc_hypo_process_mgL) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

doc_proc_mgL <- doc_processing %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "DOC_processing_mgL") %>% 
  mutate(Depth = ifelse(Depth == "mean_doc_epi_process_mgL","Epi",
                        ifelse(Depth == "mean_doc_hypo_process_mgL", "Hypo", NA)))

all_doc_mgL <- left_join(final_doc_mgL,doc_proc_mgL,by=c("DateTime","Depth"))

## Plot figure of DOC processing by year for SI: similar to figure for [DOC] by year
doc_proc_mgL <- doc_proc_mgL %>% 
  mutate(year = year(DateTime),
         month = month(DateTime))

doc_proc_mgL %>% 
  filter(month %in% c(6,7,8,9,10)) %>% 
  ggplot(mapping=aes(x=as.factor(year),y=DOC_processing_mgL,fill=Depth))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot(size=0.8,alpha=0.5)+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlab("")+
  ylab(expression(paste("DOC Internal Loading (mg L"^-1*")")))+
  theme_classic(base_size = 15)

ggsave("./Fig_Output/SI_DOC_Proc_Year_Boxplot.png",dpi=800,width=7,height=4)

###############################################################################

## Load in CTD + YSI data - temp, Sal, DO
## From merged spreadsheet in: LakeAnalyzer_thermo.R
casts <- read.csv("./Data/merged_YSI_CTD.csv") %>% 
  filter(depth >= 0) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = time)

## Calculate Epi and Hypo V.W. parameters using thermocline data
## Following Eco_DOC.R

### Create a dataframe for cast parameters at each sampling depth
depths <- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

#Initialize an empty matrix with the correct number of rows and columns 
temp <- matrix(data=NA, ncol=ncol(casts), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final <- matrix(data=NA, ncol=1, nrow=0)
dates <- unique(casts$DateTime)

#create a function to chose the matching depth closest to our focal depths
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j = dates[i]
  q <- subset(casts, casts$DateTime == j)
  
  layer1 <- q[q[, "depth"] == closest(q$depth,0.1),][1,]
  layer2 <- q[q[, "depth"] == closest(q$depth,1.6),][1,]
  layer3 <- q[q[, "depth"] == closest(q$depth,3.8),][1,]
  layer4 <- q[q[, "depth"] == closest(q$depth,5.0),][1,]
  layer5 <- q[q[, "depth"] == closest(q$depth,6.2),][1,]
  layer6 <- q[q[, "depth"] == closest(q$depth,8.0),][1,]
  layer7 <- q[q[, "depth"] == closest(q$depth,9.0),][1,]
  
  temp <- rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(casts))+1)] <- depths
  colnames(temp)[((ncol(casts))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
casts_depths <- as.data.frame(super_final) %>%
  select(-c(1,depth)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

# Separate by Env parameter of interest and pivot longer
# Temperature
temp_c <- casts_depths %>% 
  select(DateTime,new_depth,Temp_C) %>% 
  mutate(Temp_C = as.numeric(Temp_C)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Temp_C, values_fil = NA, values_fn = mean, names_prefix = "Temp_")

temp_c <- left_join(temp_c, thermo, by="DateTime")

temp_c <- temp_c %>% 
  mutate(epi_temp = ifelse(is.na(epi_bottomg_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Temp_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Temp_0.1*vol_depths$Vol_m3[1]+Temp_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_temp = ifelse(is.na(hypo_top_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Temp_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# temp_c <- temp_c[!(temp_c$DateTime = as.POSIXct("2021-08-20") | temp_c$DateTime == as.POSIXct("2020-07-08") | temp_c$DateTime == as.POSIXct("2019-05-30") | temp_c$DateTime == as.POSIXct("2019-04-29") | temp_c$DateTime == as.POSIXct("2017-09-17")),]

temp_c <- temp_c[-c(267,348,353,418,484),]

## Plot data by epi and hypo
temp_plot <- temp_c %>%  
  drop_na(epi_temp,hypo_temp) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Temp~(C^o)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_temp_c <- temp_c %>% 
  select(DateTime,epi_temp,hypo_temp) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Temp_C") %>% 
  mutate(Depth = ifelse(Depth == "epi_temp", "Epi",
                        ifelse(Depth == "hypo_temp", "Hypo", NA)))
## Dissolved oxygen
do_pSat <- casts_depths %>% 
  select(DateTime,new_depth,DO_pSat) %>% 
  mutate(DO_pSat = as.numeric(DO_pSat)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = DO_pSat, values_fil = NA, values_fn = mean, names_prefix = "DO_")

do_pSat <- left_join(do_pSat, thermo, by="DateTime")

do_pSat <- do_pSat %>% 
  mutate(epi_DO = ifelse(is.na(epi_bottomg_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, DO_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (DO_0.1*vol_depths$Vol_m3[1]+DO_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DO = ifelse(is.na(hypo_top_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, DO_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2021-08-20;479
# 2020-07-08; 413

do_mgL <- do_mgL[-c(413,479),]

do_plot <- do_pSat %>%  
  drop_na(epi_DO,hypo_DO) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab("V.W. DO %Sat")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

## Save %DO Plot for SI
ggsave("./Fig_Output/SI_PercDO.png",dpi=800,width=9,height=4)

# Format for ARIMA modeling
final_do_pSat <- do_pSat %>% 
  select(DateTime,epi_DO,hypo_DO) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_DO_mgL") %>% 
  mutate(Depth = ifelse(Depth == "epi_DO", "Epi",
                        ifelse(Depth == "hypo_DO", "Hypo", NA)))

## Calculate DO mg/L to determine days since anoxia
do_mgL <- casts_depths %>% 
  select(DateTime,new_depth,DO_mgL) %>% 
  mutate(DO_mgL = as.numeric(DO_mgL)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = DO_mgL, values_fil = NA, values_fn = mean, names_prefix = "DO_")

do_mgL <- left_join(do_mgL, thermo, by="DateTime")

do_mgL <- do_mgL %>% 
  mutate(epi_DO = ifelse(is.na(epi_bottomg_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                         ifelse(epi_bottomg_depth_m == 0.1, DO_0.1,
                                ifelse(epi_bottomg_depth_m == 1.6, (DO_0.1*vol_depths$Vol_m3[1]+DO_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                       ifelse(epi_bottomg_depth_m == 3.8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                              ifelse(epi_bottomg_depth_m == 5, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                     ifelse(epi_bottomg_depth_m == 6.2, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                            ifelse(epi_bottomg_depth_m == 8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DO = ifelse(is.na(hypo_top_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(hypo_top_depth_m == 1.6, ((DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                 ifelse(hypo_top_depth_m == 3.8, ((DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                        ifelse(hypo_top_depth_m == 5, ((DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                               ifelse(hypo_top_depth_m == 6.2, ((DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                      ifelse(hypo_top_depth_m == 8, ((DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                             ifelse(hypo_top_depth_m == 9, DO_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2021-08-20;479
# 2020-07-08; 413

do_mgL <- do_mgL[-c(413,479),]

# Format for ARIMA modeling
final_do_mgL <- do_mgL %>% 
  select(DateTime,epi_DO,hypo_DO) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_DO_mgL") %>% 
  mutate(Depth = ifelse(Depth == "epi_DO", "Epi",
                        ifelse(Depth == "hypo_DO", "Hypo", NA)))

# Calculate days since anoxia and oxygenation status for Hypo
hypo_do_mgL <- final_do_mgL %>% 
  filter(Depth == "Hypo") %>% 
  mutate(anoxia = ifelse(VW_DO_mgL < 1.0, 1, 0)) %>% 
  mutate(anoxia_time_d = 0) 

for (i in 1:length(hypo_do_mgL$DateTime)){
  if (hypo_do_mgL$anoxia[i] == 1){
    a = seq(from = hypo_do_mgL$DateTime[i-1], to = hypo_do_mgL$DateTime[i], by = 'day')
    hypo_do_mgL$anoxia_time_d[i] = hypo_do_mgL$anoxia_time_d[i-1]+length(a)
  } else {
    hypo_do_mgL$anoxia_time_d[i] = 0
  }
}

# Plot
ggplot(hypo_do_mgL,mapping=aes(x=DateTime,y=anoxia_time_d))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=0.75)+
  geom_point(size=2)+
  xlab("")+
  ylab("Time since anoxia (d)")+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

ggsave("./Fig_Output/SI_DaysAnoxia.jpg",width=7,height=5,units="in",dpi=320)

## Load in model results and plot by oxic vs. anoxic waters in the hypolimnion
doc_model <- read.csv("./Fig_Output/model_results.csv") %>% 
  select(-X) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

doc_model_oxy <- left_join(hypo_do_mgL,doc_model,by="DateTime") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(4,5,6,7,8,9,10,11)) %>% 
  drop_na(anoxia)

## Plot by hypo internal loading under oxic vs. anoxic waters
wilcox.test(mean_doc_hypo_process_g~anoxia,doc_model_oxy)

doc_hypo <- doc_model_oxy %>% 
  ggplot(mapping=aes(x=as.character(anoxia),y=mean_doc_hypo_process_g/1000))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  ylab(expression(paste("Hypo Internal DOC (kg)")))+
  xlab("")+
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("Oxic", "Anoxic"))+
  theme_classic(base_size = 15)

doc_depth <- doc_model %>% 
  select(DateTime,mean_doc_epi_process_g,mean_doc_hypo_process_g) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(4,5,6,7,8,9,10,11)) %>% 
  pivot_longer(!c(DateTime,month),names_to="Loc",values_to="mean_doc_process_g") %>% 
  ggplot(mapping=aes(x=Loc,y=mean_doc_process_g/1000))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_boxplot()+
  ylab(expression(paste("Internal DOC (kg)")))+
  xlab("")+
  scale_x_discrete(breaks=c("mean_doc_epi_process_g","mean_doc_hypo_process_g"),
                   labels=c("Epi", "Hypo"))+
  theme_classic(base_size = 15)

ggarrange(doc_hypo,doc_depth,nrow=1,ncol=2,labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/SI_Internal_proc_comps.jpg",width=7,height=4,units="in",dpi=320)

## Contribution of Epi internal DOC loading as compared to total loading
doc_model_select <- doc_model %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  mutate(month = month(DateTime)) %>% 
  filter(month %in% c(4,5,6,7,8,9,10,11))

mean(doc_model_select$mean_doc_epi_process_g)/(mean(doc_model_select$mean_doc_epi_process_g)+mean(doc_model_select$mean_doc_hypo_outflow_g)+(mean(doc_model_select$mean_doc_inflow_g)*0.74))*100

mean(doc_model_select$mean_doc_hypo_process_g)/(mean(doc_model_select$mean_doc_hypo_process_g)+(mean(doc_model_select$mean_doc_inflow_g)*0.26))*100

## Conductivity
cond_uScm <- casts_depths %>% 
  select(DateTime,new_depth,Cond_uScm) %>% 
  mutate(Cond_uScm = as.numeric(Cond_uScm)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Cond_uScm, values_fil = NA, values_fn = mean, names_prefix = "Cond_")

cond_uScm <- left_join(cond_uScm, thermo, by="DateTime")

cond_uScm <- cond_uScm %>% 
  mutate(epi_Cond = ifelse(is.na(epi_bottomg_depth_m), ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                         ifelse(epi_bottomg_depth_m == 0.1, Cond_0.1,
                                ifelse(epi_bottomg_depth_m == 1.6, (Cond_0.1*vol_depths$Vol_m3[1]+Cond_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                       ifelse(epi_bottomg_depth_m == 3.8, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                              ifelse(epi_bottomg_depth_m == 5, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                     ifelse(epi_bottomg_depth_m == 6.2, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                            ifelse(epi_bottomg_depth_m == 8, ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Cond = ifelse(is.na(hypo_top_depth_m), ((Cond_0.1*vol_depths$Vol_m3[1])+(Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(hypo_top_depth_m == 1.6, ((Cond_1.6*vol_depths$Vol_m3[2])+(Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                 ifelse(hypo_top_depth_m == 3.8, ((Cond_3.8*vol_depths$Vol_m3[3])+(Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                        ifelse(hypo_top_depth_m == 5, ((Cond_5.0*vol_depths$Vol_m3[4])+(Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                               ifelse(hypo_top_depth_m == 6.2, ((Cond_6.2*vol_depths$Vol_m3[5])+(Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                      ifelse(hypo_top_depth_m == 8, ((Cond_8.0*vol_depths$Vol_m3[6])+(Cond_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                             ifelse(hypo_top_depth_m == 9, Cond_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2017-03-12; 208
cond_uScm <- cond_uScm[-c(208),]

cond_uScm %>%  
  drop_na(epi_Cond,hypo_Cond) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Cond,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Cond,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Cond,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Cond,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Cond.~(~mu*s~cm^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_cond_uScm <- cond_uScm %>% 
  select(DateTime,epi_Cond,hypo_Cond) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Cond_uScm") %>% 
  mutate(Depth = ifelse(Depth == "epi_Cond", "Epi",
                        ifelse(Depth == "hypo_Cond", "Hypo", NA)))

# Turbidity
turb_ntu <- casts_depths %>% 
  select(DateTime,new_depth,Turb_NTU) %>% 
  mutate(Turb_NTU = as.numeric(Turb_NTU)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Turb_NTU, values_fil = NA, values_fn = mean, names_prefix = "Turb_")

turb_ntu <- left_join(turb_ntu, thermo, by="DateTime")

turb_ntu <- turb_ntu %>% 
  mutate(epi_Turb = ifelse(is.na(epi_bottomg_depth_m), ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Turb_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Turb_0.1*vol_depths$Vol_m3[1]+Turb_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Turb = ifelse(is.na(hypo_top_depth_m), ((Turb_0.1*vol_depths$Vol_m3[1])+(Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Turb_1.6*vol_depths$Vol_m3[2])+(Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Turb_3.8*vol_depths$Vol_m3[3])+(Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Turb_5.0*vol_depths$Vol_m3[4])+(Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Turb_6.2*vol_depths$Vol_m3[5])+(Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Turb_8.0*vol_depths$Vol_m3[6])+(Turb_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Turb_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Some light QA/QC'ing
# 2019-07-03; 208
turb_ntu <- turb_ntu[-c(270),]

turb_ntu %>%  
  drop_na(epi_Turb,hypo_Turb) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Turb,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Turb,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Turb,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Turb,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylim(0,50)+
  ylab(expression(V.W.~Turb.~(~NTU)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_turb_ntu <- turb_ntu %>% 
  select(DateTime,epi_Turb,hypo_Turb) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Turb_NTU") %>% 
  mutate(Depth = ifelse(Depth == "epi_Turb", "Epi",
                        ifelse(Depth == "hypo_Turb", "Hypo", NA)))

###############################################################################

## Load in Flora data - Chla and community analysis
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/5/a24f4dbe9f0d856688f257547069d0a3" 
#infile1 <- paste0(getwd(),"/Data/FluoroProbe.csv")
#download.file(inUrl1,infile1,method="curl")
flora <- read.csv("./Data/FluoroProbe_2014_2021.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2015-01-01")) %>% 
  filter(Site == 50)

## Create dataframe of Flora data that is closest to pre-defined sampling depths; then calculate epi and hypo VW Chla concentrations
#Initialize an empty matrix with the correct number of rows and columns 
temp <- matrix(data=NA, ncol=ncol(flora), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final <- matrix(data=NA, ncol=1, nrow=0)
dates <- unique(flora$DateTime)

library(plyr) #only use plyr for this for loop, then detach!

#For loop to retrieve CTD depth with the closest function and fill in matrix
for (i in 1:length(dates)){
  j = dates[i]
  q <- subset(flora, flora$DateTime == j)
  
  layer1 <- q[q[, "Depth_m"] == closest(q$Depth_m,0.1),][1,]
  layer2 <- q[q[, "Depth_m"] == closest(q$Depth_m,1.6),][1,]
  layer3 <- q[q[, "Depth_m"] == closest(q$Depth_m,3.8),][1,]
  layer4 <- q[q[, "Depth_m"] == closest(q$Depth_m,5.0),][1,]
  layer5 <- q[q[, "Depth_m"] == closest(q$Depth_m,6.2),][1,]
  layer6 <- q[q[, "Depth_m"] == closest(q$Depth_m,8.0),][1,]
  layer7 <- q[q[, "Depth_m"] == closest(q$Depth_m,9.0),][1,]
  
  temp <- rbind(layer1,layer2,layer3,layer4,layer5,layer6,layer7)
  temp[,((ncol(flora))+1)] <- depths
  colnames(temp)[((ncol(flora))+1)]<-"new_depth"
  final <- temp
  final <- data.frame(final)
  super_final <- rbind.fill.matrix(super_final,final)
}

detach(package:plyr)#to prevent issues with dplyr vs plyr not playing well together!

#now need to clean up the data frame and make all factors numeric
flora_depths <- as.data.frame(super_final) %>%
  select(-c(1,Depth_m)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Calculate VW Epi and Hypo concentrations
chla_ugL <- flora_depths %>% 
  select(DateTime,new_depth,TotalConc_ugL) %>% 
  mutate(TotalConc_ugL = as.numeric(TotalConc_ugL)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = TotalConc_ugL, values_fil = NA, values_fn = mean, names_prefix = "Chla_")

chla_ugL <- left_join(chla_ugL, thermo, by="DateTime")

chla_ugL <- chla_ugL %>% 
  mutate(epi_Chla = ifelse(is.na(epi_bottomg_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottomg_depth_m == 0.1, Chla_0.1,
                                  ifelse(epi_bottomg_depth_m == 1.6, (Chla_0.1*vol_depths$Vol_m3[1]+Chla_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottomg_depth_m == 3.8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottomg_depth_m == 5, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottomg_depth_m == 6.2, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottomg_depth_m == 8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Chla = ifelse(is.na(hypo_top_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Chla_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Plot
chla_plot <- chla_ugL %>%  
  drop_na(epi_Chla,hypo_Chla) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(V.W.~Chla~(~mu*g~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_chla_ugL <- chla_ugL %>% 
  select(DateTime,epi_Chla,hypo_Chla) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Chla_ugL") %>% 
  mutate(Depth = ifelse(Depth == "epi_Chla", "Epi",
                        ifelse(Depth == "hypo_Chla", "Hypo", NA)))

###############################################################################

## Load in inflow
inflow <- read.csv("./Data/Inflow_2013_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  select(Reservoir:VT_Temp_C)

# Find relationship between WVWA and VT flow
# Use to fill in any missing WVWA values with VT flow (where possible)
flow_lm <- lm(WVWA_Flow_cms ~ VT_Flow_cms, data = inflow)

inflow <- inflow %>% 
  mutate(WVWA_Flow_cms = ifelse(is.na(WVWA_Flow_cms), 0.8917528*VT_Flow_cms-0.0009302, WVWA_Flow_cms)) %>% 
  mutate(flow_diff = abs(VT_Flow_cms - WVWA_Flow_cms))

# Average inflow by day
inflow_daily <- inflow %>% 
  group_by(DateTime) %>% 
  summarize_at(vars("WVWA_Flow_cms"),funs(mean(.,na.rm=TRUE),sd)) %>% 
  filter(DateTime >= as.POSIXct("2015-01-01"))

# Plot daily inflow for the study period
inflow_plot <- inflow_daily %>% 
  na.omit(mean) %>% 
  ggplot(mapping=aes(x=DateTime,y=mean))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  geom_point(size=2)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Inflow~(~m^3~s^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_inflow_m3s <- inflow_daily %>% 
  select(DateTime,mean) %>% 
  dplyr::rename(Inflow_m3s = mean)

# Calculate residence time assuming full pond
wtr_d <- final_inflow_m3s %>% 
  mutate(wtr_d = (310000)/Inflow_m3s/60/60/24)

wtr_d %>% 
  na.omit(wtr_d) %>% 
  summarise(med = median(wtr_d),
            sd = sd(wtr_d)) %>% 
  mutate(se = sd/length(wtr_d$wtr_d))

## Plot and calculate mean and SD
wtr_d %>% 
  na.omit(wtr_d) %>% 
  ggplot(mapping=aes(x=DateTime,y=wtr_d))+
  geom_hline(yintercept = 139, linetype="dashed")+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  geom_point(size=2)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Water~Residence~Time~(d)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

###############################################################################

## Plot Temp, DO, Chla, and Inflow for main MS
ggarrange(temp_plot,do_plot,chla_plot,inflow_plot,ncol=1,nrow=4,common.legend = TRUE, labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/Fig4_EnvParameters.jpg",width=10,height=12,units="in",dpi=320)
  
###############################################################################

## Load in met data - shortwave radiation and rainfall
# Downloaded on 23 Aug 2022
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb"
#infile1 <- paste0(getwd(),"/Data/FCR_Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

## Load in Met data
met <- read.csv("./Data/FCR_Met_final_2015_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01 00:00:00"))

## Average to daily
met_daily <- met %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  dplyr::summarise(rain_tot_mm = sum(Rain_Total_mm, na.rm = TRUE),
            ShortwaveRadiationUp_Average_W_m2 = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE))

## Then plot total rainfall and shortwave radiation
rain_plot <- met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=rain_tot_mm))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Total~Rainfall~(mm)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

sw_plot <- met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=ShortwaveRadiationUp_Average_W_m2))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(S.W.~Radiation~(W~m^2)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(rain_plot,sw_plot,ncol=1,nrow=2,common.legend = TRUE, labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

ggsave("./Fig_Output/SI_MetParameters.jpg",width=10,height=7,units="in",dpi=320)

###############################################################################

## Format and add thermocline depth as a potential predictor variable, too
final_thermo <- thermo %>% 
  select(DateTime, thermo.depth)

###############################################################################

###############################################################################

## Organize data for ARIMA modeling - following EddyFlux, 3_Rev_EnvAnalysis.R
# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, and SW Radiation for Epi and Hypo
arima_epi <- plyr::join_all(list(all_doc_mgL,final_temp_c,final_do_pSat,final_chla_ugL),by=c("DateTime","Depth"),type="left") 

# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, SW Radiation, CO2, and CH4 for Epi and Hypo
arima_hypo <- plyr::join_all(list(all_doc_mgL,final_temp_c,final_do_pSat,final_chla_ugL),by=c("DateTime","Depth"),type="left")

## Select time points where we have DOC concentrations
arima_epi <- plyr::join_all(list(arima_epi,met_daily,final_inflow_m3s),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Epi")

arima_hypo <- plyr::join_all(list(arima_hypo,met_daily,final_inflow_m3s),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Hypo")

## Add in days since anoxia for Hypo
arima_hypo <- left_join(arima_hypo,hypo_do_mgL,by=c("DateTime","Depth")) %>% 
  dplyr::select(DateTime,Depth,VW_DOC_mgL,DOC_processing_mgL,VW_Temp_C,VW_DO_mgL.x,VW_Chla_ugL,rain_tot_mm,
                ShortwaveRadiationUp_Average_W_m2,Inflow_m3s,anoxia_time_d) %>% 
  dplyr::rename(VW_DO_pSat = VW_DO_mgL.x)

## Add in thermocline depth information
arima_epi <- left_join(arima_epi,final_thermo,by="DateTime")

arima_hypo <- left_join(arima_hypo,final_thermo,by="DateTime")

## Calculate stats for env parameters - limited to summer stratified period (May-Oct)
epi_stats <- arima_epi %>% 
  mutate(month = month(DateTime)) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & month %in% c(5,6,7,8,9,10)) %>%
  select(VW_Temp_C,VW_DO_mgL,VW_Chla_ugL,rain_tot_mm,ShortwaveRadiationUp_Average_W_m2,Inflow_m3s) %>% 
  summarise_all(list(min,max,median,mean,sd),na.rm=TRUE)

hypo_stats <- arima_hypo %>% 
  mutate(month = month(DateTime)) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & month %in% c(5,6,7,8,9,10)) %>% 
  select(VW_Temp_C,VW_DO_mgL,VW_Chla_ugL,rain_tot_mm,ShortwaveRadiationUp_Average_W_m2,Inflow_m3s) %>% 
  summarise_all(list(min,max,median,mean,sd),na.rm=TRUE)

###############################################################################

## Plot Hypo DOC concentrations under oxic vs. anoxic conditions
arima_hypo %>% 
  mutate(anoxia = ifelse(VW_DO_mgL >= 1.0, "Oxic",
                         ifelse(VW_DO_mgL < 1.0, "Anoxic", NA))) %>% 
  filter(anoxia == "Oxic" | anoxia == "Anoxic") %>% 
  ggplot(mapping=aes(x=anoxia,y=VW_DOC_mgL,color=anoxia))+
  geom_boxplot()+
  theme_classic(base_size = 15)

###############################################################################

## Check correlations among environmental variables - what needs to be removed?
# Epi
epi_cor = as.data.frame(cor(arima_epi[,3:11],use = "complete.obs"),method=c("pearson"))
write_csv(epi_cor, "./Fig_Output/22Dec22_epi_cor.csv")

chart.Correlation(arima_epi[,3:11],histogram = TRUE,method=c("pearson"))
# No correlations!

# Hypo
hypo_cor = as.data.frame(cor(arima_hypo[,3:12],use = "complete.obs"),method=c("pearson"))
write_csv(hypo_cor, "./Fig_Output/22Dec22_hypo_cor.csv")

chart.Correlation(arima_hypo[,3:12],histogram = TRUE,method=c("pearson"))
# No correlations!

###############################################################################

## Check for skewness following MEL script!
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# Epi
for (i in 3:11){
  print(colnames(arima_epi)[i])
  var <- arima_epi[,i]
  hist(as.matrix(var), main = colnames(arima_epi)[i])
  print(skewness(arima_epi[,i], na.rm = TRUE))
  print(skewness(log(arima_epi[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(arima_epi[,i]), na.rm = TRUE))
  print(skewness(arima_epi[,i]^2),na.rm=TRUE)
  var <- log(arima_epi[,i])
  hist(as.matrix(var), main = c("Log",colnames(arima_epi)[i]))
  var <- Math.cbrt(arima_epi[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(arima_epi)[i]))
  var <- (arima_epi[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(arima_epi)[i]))
}
# Nothing: DOC_processing, Temp, DO, SW Radiation, Rain, Thermo
# Log: DOC, Chla, Inflow

# Transform and scale data
arima_epi_scale <- arima_epi %>% 
  mutate(VW_DOC_mgL = log(VW_DOC_mgL),
         VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s))

arima_epi_scale[,3:11] <- scale(arima_epi_scale[,3:11])

# Hypo
for (i in 3:12){
  print(colnames(arima_hypo)[i])
  var <- arima_hypo[,i]
  hist(as.matrix(var), main = colnames(arima_hypo)[i])
  print(skewness(arima_hypo[,i], na.rm = TRUE))
  print(skewness(log(arima_hypo[,i]+0.0001), na.rm = TRUE))
  print(skewness(Math.cbrt(arima_hypo[,i]), na.rm = TRUE))
  print(skewness(arima_hypo[,i]^2),na.rm=TRUE)
  var <- log(arima_hypo[,i])
  hist(as.matrix(var), main = c("Log",colnames(arima_hypo)[i]))
  var <- Math.cbrt(arima_hypo[,i])
  hist(as.matrix(var), main = c("cube_rt",colnames(arima_hypo)[i]))
  var <- (arima_hypo[,i]^2)
  hist(as.matrix(var), main = c("sq",colnames(arima_hypo)[i]))
}
# Nothing: DOC, Temp, DO, Rain, SW Radiation, Thermo, Anoxia_time
# Log: Chla, Inflow
# Cube root: DOC_processing

# Transform and scale data
arima_hypo_scale <- arima_hypo %>% 
  mutate(DOC_processing_mgL = Math.cbrt(DOC_processing_mgL),
         VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s))

arima_hypo_scale[,3:12] <- scale(arima_hypo_scale[,3:12])

###############################################################################

## ARIMA modeling - DOC concentrations!
# Following MEL code : )
# THINGS TO CHANGE: 'cols' (change to the environmental variables); best fit!

###############################################################################

# Epi
colnames(arima_epi_scale)

cols <- c(5:11) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_epi_scale[,3] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "epi"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "epi"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for epi",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "epi"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_epi_scale)[combn(cols,4)[,16]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,16] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_epi_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_epi_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_epi_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_epi_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################

# Hypo
colnames(arima_hypo_scale)

cols <- c(5:7,10:12) # UPDATE THIS TO THE ENV. VARIABLES - exclude Rainfall and SW radiation (primarily surface processes)
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,3] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "hypo"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "hypo"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for hypo",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "hypo"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_hypo_scale)[combn(cols,5)[,2]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,5)[,2] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_hypo_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################

## Epi DOC Processing
colnames(arima_epi_scale)

cols <- c(5:11) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,4] # UPDATE THIS TO DOC PROCESSING

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "epi"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "epi"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for epi",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "hypo"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_hypo_scale)[combn(cols,3)[,6]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,3)[,6] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_hypo_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################
## Hypo DOC processing
colnames(arima_hypo_scale)

cols <- c(5:7,10:12) # UPDATE THIS TO THE ENV. VARIABLES - exclude Rainfall and SW radiation (primarily surface processes)
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,4] # UPDATE THIS TO DOC PROCESSING

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "hypo"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "hypo"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for hypo",sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "hypo"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(arima_hypo_scale)[combn(cols,4)[,13]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,13] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_hypo_scale[,1]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,1])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

###############################################################################

## Random plots

## Plot [DOC] at 9 m under anoxic vs. oxic conditions
# First combine DOC and DO concentrations at 9m
casts_depths <- casts_depths %>% 
  rename(Depth_m = new_depth) 

casts_depths <- casts_depths %>% 
  mutate(Depth_m = as.numeric(Depth_m),
         DO_mgL= as.numeric(DO_mgL))

hypo_9m <- left_join(chem_50,casts_depths,by=c("DateTime","Depth_m")) %>% 
  select(DateTime,Depth_m,DOC_mgL,DO_mgL) %>% 
  filter(Depth_m == 9) %>% 
  mutate(oxy = ifelse(DO_mgL < 1.0, "Anoxic",
                      ifelse(DO_mgL >= 1.0, "Oxic", NA))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  drop_na() %>% 
  mutate(year = year(DateTime))

ggplot(hypo_9m,mapping=aes(x=oxy,y=DOC_mgL))+
  geom_boxplot()

###############################################################################