###############################################################################

### Script to look at environmental variables and conduct ARIMA modeling between 
### epi and hypo DOC concentrations and various limno. and met. parameters
### 18 Aug 2022, A. Hounshell

### Adding in various metrics of anoxia/oxygenation to hypo AR model
### 1 Sep 2022, A. Hounshell

### Updating following comments from CCC - DO mg/L -> DO%; remove GHG variables; ARIMA for DOC processing
### Adding in some additional visualizations/analyses for the Environmental data

### Cleaned-up for code review and final verison of MS
### 12 Mar 2024, A. Hounshell

### Last modified: 12 Apr. 2024, A. Hounshell
## Use bathy data from EDI
## Include Epi phyto biomass in Hypo ARIMA models (for 'sinking' phytos!)
## Updates following DWH code-review (thanks, Dexter!!)

### Modified: 14 Feb. 2025, A. Hounshell
## Updating to include model output which incorporates Falling Creek (FC) inputs
## Include total nutrients (TN, TP) at deepest point as explantory variables
## in ARIMA modeling
## Use PACF to ID important AR lags

## Modified: 01 August 2025, A. Hounshell
## Updating to include updated discharge uncertainty estimates for model
## input

###############################################################################
## Clear workspace
rm(list = ls())

## Set working directory
wd <- getwd()
setwd(wd)

## Load in libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,zoo,scales,plyr,ggridges,
               lubridate,lognorm,forecast,utils,igraph,RColorBrewer,PerformanceAnalytics)

###############################################################################
## Load in Epi and Hypo V.W. DOC concentrations - from Eco_DOC_rlnorm_FC.R
doc_mgL <- read.csv("./Data/EpiHypo_Weir_FC_DOC.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(Depth = Loc)

###############################################################################
## Add in DOC processing - calculated from Eco_DOC_rlnorm_FC.R
doc_processing <- read_csv("./Data/01Aug25_final_doc_inputs_fc.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

doc_proc_mgL <- doc_processing %>% 
  select(DateTime, mean_doc_epi_process_mgL, mean_doc_hypo_process_mgL) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "DOC_processing_mgL") %>% 
  mutate(Depth = ifelse(Depth == "mean_doc_epi_process_mgL","Epi",
                        ifelse(Depth == "mean_doc_hypo_process_mgL", "Hypo", NA))) %>% 
  mutate(year = year(DateTime),
         month = month(DateTime))

all_doc_mgL <- left_join(doc_mgL,doc_proc_mgL,by=c("DateTime","Depth"))

doc_proc_g <- doc_processing %>% 
  select(DateTime, mean_doc_epi_process_g, mean_doc_hypo_process_g) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "DOC_processing_g") %>% 
  mutate(Depth = ifelse(Depth == "mean_doc_epi_process_g","Epi",
                        ifelse(Depth == "mean_doc_hypo_process_g", "Hypo", NA))) %>% 
  mutate(year = year(DateTime),
         month = month(DateTime))

## Plot distributions of doc_proc_g - for summer stratified period ONLY (May 1-Nov. 15)
summer_epi_processing <- doc_proc_g %>% 
  filter(Depth == "Epi") %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320)

## Test statistical significance of distributions from zero
t.test(summer_epi_processing$DOC_processing_g[summer_epi_processing$year==2017],mu=0)
t.test(summer_epi_processing$DOC_processing_g[summer_epi_processing$year==2018],mu=0) # Sig. dif from zero (p<0.05)
t.test(summer_epi_processing$DOC_processing_g[summer_epi_processing$year==2019],mu=0)
t.test(summer_epi_processing$DOC_processing_g[summer_epi_processing$year==2020],mu=0) 
t.test(summer_epi_processing$DOC_processing_g[summer_epi_processing$year==2021],mu=0) # Sig. diff rom zero (p<0.05)
t.test(summer_epi_processing$DOC_processing_g,mu=0) # Sig. diff rom zero (p<0.05)

epi_distribution <- summer_epi_processing %>% 
  ggplot(mapping=aes(x=DOC_processing_g/1000,y=as.factor(year),fill=as.factor(year)))+
  ggtitle("Epilimnion")+
  stat_density_ridges(quantile_lines=TRUE,quantiles=2)+
  geom_vline(mapping=aes(xintercept=0),linetype="dashed")+
  xlab(expression(paste("Internal DOC (kg d"^-1*")")))+
  ylab("")+
  xlim(-242,255)+
  scale_y_discrete(labels=c("2017","*2018","2019","2020","*2021"))+
  scale_fill_manual(values=c("#E7804B","#F0B670","#91B374","#7EBDC2", "#9B9B9B"))+
  theme_ridges()+
  theme(legend.position = "none")

summer_hypo_processing <- doc_proc_g %>% 
  filter(Depth == "Hypo") %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320)

## Test statistical significance of distributions from zero
t.test(summer_hypo_processing$DOC_processing_g[summer_hypo_processing$year==2017],mu=0)
t.test(summer_hypo_processing$DOC_processing_g[summer_hypo_processing$year==2018],mu=0)
t.test(summer_hypo_processing$DOC_processing_g[summer_hypo_processing$year==2019],mu=0)
t.test(summer_hypo_processing$DOC_processing_g[summer_hypo_processing$year==2020],mu=0)
t.test(summer_hypo_processing$DOC_processing_g[summer_hypo_processing$year==2021],mu=0)
t.test(summer_hypo_processing$DOC_processing_g,mu=0)

hypo_distribution <-  summer_hypo_processing %>% 
  ggplot(mapping=aes(x=DOC_processing_g/1000,y=as.factor(year),fill=as.factor(year)))+
  ggtitle("Hypolimnion")+
  stat_density_ridges(quantile_lines=TRUE,quantiles=2)+
  geom_vline(mapping=aes(xintercept=0),linetype="dashed")+
  xlab(expression(paste("Internal DOC (kg d"^-1*")")))+
  ylab("")+
  scale_fill_manual(values=c("#E7804B","#F0B670","#91B374","#7EBDC2", "#9B9B9B"))+
  theme_ridges()+
  theme(legend.position = "none")

### Load in model summary from 2_Eco_DOC_rlnorm_FC.R
all_summary <- read.csv("./Data/Table_S5_model_summary_fc.csv")

## Plot
func_order <- c("mean_doc_entr_g","mean_doc_inflow_g_comb", "mean_doc_hypo_outflow_g", 
                "mean_doc_epi_outflow_g","mean_doc_epi_process_g","mean_doc_hypo_process_g")

model_summary <- all_summary %>% 
  filter(! func %in% c("mean_doc_epi_process_mgL","mean_doc_hypo_process_mgL","mean_doc_dt_epi_g","mean_doc_dt_hypo_g","mean_doc_inflow_g","mean_doc_fc_inflow_g")) %>% 
  ggplot(mapping=aes(x=factor(func,func_order),y=mean/1000,fill=year))+
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_hline(mapping=aes(yintercept=0))+
  ylab(expression(paste("Mean DOC (kg d"^-1*")")))+
  scale_fill_manual(values=c("#E7804B","#F0B670","#91B374","#7EBDC2", "#9B9B9B", "#393E41"))+
  scale_x_discrete("Model Term",labels=c("Entr.","Inflow",expression(atop("Hypo.", "Outflow")),expression(atop("Epi.", "Outflow")),
                                              expression(atop("Epi.", "Internal")),expression(atop("Hypo.", "Internal"))))+
  theme_bw(base_size = 15) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="top")

## Calculate Percent contribution of internal vs. external loading
contributions <- all_summary %>% 
  select(func,year,mean) %>% 
  pivot_wider(names_from = func,values_from=mean) %>% 
  mutate(epi_internal = abs(mean_doc_epi_process_g)/(mean_doc_inflow_g_comb+abs(mean_doc_epi_process_g)+mean_doc_hypo_process_g)*100,
         hypo_internal = mean_doc_hypo_process_g/(mean_doc_inflow_g_comb+abs(mean_doc_epi_process_g)+mean_doc_hypo_process_g)*100,
         inflow = mean_doc_inflow_g_comb/(mean_doc_inflow_g_comb+abs(mean_doc_epi_process_g)+mean_doc_hypo_process_g)*100) %>% 
  select(year,epi_internal,hypo_internal,inflow) %>% 
  pivot_longer(cols = epi_internal:inflow,
               names_to = "type")

write.csv(contributions, "./Data/Table_S6_model_contributions_fc.csv",row.names=FALSE)

## Plot overall contributions
type_order <- c('hypo_internal','epi_internal','inflow')

contribution_fig <- ggplot(contributions,mapping=aes(x=year,y=value,fill=factor(type,type_order)))+
  geom_bar(stat="identity",position="stack",color="black",alpha=0.8)+
  ylab("Mean Contribution (%)")+
  scale_fill_manual(values=c("#393E41","#7EBDC2","#F0B670"),labels=c("Hypo. Internal","Epi. Internal", "Inflow"))+
  scale_x_discrete("")+
  theme_bw(base_size = 15) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="top")

ggarrange(ggarrange(model_summary,contribution_fig,ncol=2,labels=c("A.","B."),font.label=list(face="plain",size=15)),
          epi_distribution,hypo_distribution,
          nrow=3,ncol=1,labels = c("","C.", "D."),font.label=list(face="plain",size=15),heights = c(1, 0.7,0.7))

ggsave("./Figs/Fig5_Model_Summary_FC.jpg",width=9,height=12,units="in",dpi=320)

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

## Check average thermo depth from April-Oct from 2017-2021
thermo %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  summarise_all(median,na.rm=TRUE)

## Check by year, too!
thermo %>% 
  mutate(year = year(DateTime),
         doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  group_by(year) %>% 
  summarise_all(median,na.rm=TRUE)

###############################################################################
## Load in CTD + YSI data - temp, Sal, DO
## From merged spreadsheet in: LakeAnalyzer_thermo.R
casts <- read.csv("./Data/merged_YSI_CTD.csv") %>% 
  filter(depth >= 0) %>% 
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>% 
  dplyr::rename(DateTime = time)

## Calculate Epi and Hypo V.W. parameters using thermocline data
## Following Eco_DOC_rlnorm.R

### Create a dataframe for cast parameters at each sampling depth
depths <- c(0.1, 1.6, 3.8, 5.0, 6.2, 8.0, 9.0) 

# Create vector of different volumes for each depth: from Eco_DOC_rlnorm.R
# Table SI.3
vol_depths <- read.csv("./Data/vol_depths.csv")

# Initialize an empty matrix with the correct number of rows and columns 
temp <- matrix(data=NA, ncol=ncol(casts), nrow=length(depths)) #of cols in CTD data, and then nrows = # of layers produced
super_final <- matrix(data=NA, ncol=1, nrow=0)
dates <- unique(casts$DateTime)

# Create a function to choose the matching depth closest to our focal depths
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

library(plyr) #only use plyr for this for loop, then detach!

# For loop to retrieve CTD depth with the closest function and fill in matrix
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

detach(package:plyr) # to prevent issues with dplyr vs plyr not playing well together!

# Now need to clean up the data frame and make all factors numeric
casts_depths <- as.data.frame(super_final) %>%
  select(-c(1,depth)) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  drop_na(Temp_C)

# Separate by Env parameter of interest and pivot longer
# Temperature
temp_c <- casts_depths %>% 
  select(DateTime,new_depth,Temp_C) %>% 
  mutate(Temp_C = as.numeric(Temp_C)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = new_depth, values_from = Temp_C, values_fill = NA, values_fn = mean, names_prefix = "Temp_")

temp_c <- left_join(temp_c, thermo, by="DateTime")

temp_c <- temp_c %>% 
  mutate(epi_temp = ifelse(is.na(epi_bottom_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottom_depth_m == 0.1, Temp_0.1,
                                  ifelse(epi_bottom_depth_m == 1.6, (Temp_0.1*vol_depths$Vol_m3[1]+Temp_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottom_depth_m == 3.8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottom_depth_m == 5, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottom_depth_m == 6.2, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottom_depth_m == 8, ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_temp = ifelse(is.na(hypo_top_depth_m), ((Temp_0.1*vol_depths$Vol_m3[1])+(Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Temp_1.6*vol_depths$Vol_m3[2])+(Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Temp_3.8*vol_depths$Vol_m3[3])+(Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Temp_5.0*vol_depths$Vol_m3[4])+(Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Temp_6.2*vol_depths$Vol_m3[5])+(Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Temp_8.0*vol_depths$Vol_m3[6])+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Temp_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Fill in some 2021 data where we're missing 'in between' depths
temp_c <- temp_c %>% 
  mutate(epi_temp = ifelse(is.na(epi_temp) & epi_bottom_depth_m == 5, ((Temp_0.1*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(Temp_5.0*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_temp),
         hypo_temp = ifelse(is.na(hypo_temp) & hypo_top_depth_m == 6.2, ((Temp_8.0*sum(vol_depths$Vol_m3[5]+vol_depths$Vol_m3[6]))+(Temp_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),hypo_temp))

## Some light QA/QC'ing - weird temps
# temp_c <- temp_c[!(temp_c$DateTime = as.POSIXct("2021-08-20") | temp_c$DateTime == as.POSIXct("2020-07-08") | temp_c$DateTime == as.POSIXct("2019-05-30") | temp_c$DateTime == as.POSIXct("2019-04-29") | temp_c$DateTime == as.POSIXct("2017-09-17")),]
temp_c <- temp_c[-c(267,348,353,418,484),]

## Plot data by epi and hypo
temp_plot <- temp_c %>%  
  drop_na(epi_temp,hypo_temp) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_temp,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_temp,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(VW~Temp~(C^o)))+
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
  pivot_wider(names_from = new_depth, values_from = DO_pSat, values_fill = NA, values_fn = mean, names_prefix = "DO_")

do_pSat <- left_join(do_pSat, thermo, by="DateTime")

do_pSat <- do_pSat %>% 
  mutate(epi_DO = ifelse(is.na(epi_bottom_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottom_depth_m == 0.1, DO_0.1,
                                  ifelse(epi_bottom_depth_m == 1.6, (DO_0.1*vol_depths$Vol_m3[1]+DO_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottom_depth_m == 3.8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottom_depth_m == 5, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottom_depth_m == 6.2, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottom_depth_m == 8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DO = ifelse(is.na(hypo_top_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, DO_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Fill in some 2021 data where we're missing 'in between' depths
do_pSat <- do_pSat %>% 
  mutate(epi_DO = ifelse(is.na(epi_DO) & epi_bottom_depth_m == 5, ((DO_0.1*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(DO_5.0*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_DO),
         hypo_DO = ifelse(is.na(hypo_DO) & hypo_top_depth_m == 6.2, ((DO_8.0*sum(vol_depths$Vol_m3[5]+vol_depths$Vol_m3[6]))+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),hypo_DO))

## Some light QA/QC'ing
# 2021-08-20; 479
# 2020-07-08; 413
do_pSat <- do_pSat[-c(413,479),]

do_plot <- do_pSat %>%  
  drop_na(epi_DO,hypo_DO) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab("VW DO %Sat")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

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
  pivot_wider(names_from = new_depth, values_from = DO_mgL, values_fill = NA, values_fn = mean, names_prefix = "DO_")

do_mgL <- left_join(do_mgL, thermo, by="DateTime")

do_mgL <- do_mgL %>% 
  mutate(epi_DO = ifelse(is.na(epi_bottom_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                         ifelse(epi_bottom_depth_m == 0.1, DO_0.1,
                                ifelse(epi_bottom_depth_m == 1.6, (DO_0.1*vol_depths$Vol_m3[1]+DO_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                       ifelse(epi_bottom_depth_m == 3.8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                              ifelse(epi_bottom_depth_m == 5, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                     ifelse(epi_bottom_depth_m == 6.2, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                            ifelse(epi_bottom_depth_m == 8, ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_DO = ifelse(is.na(hypo_top_depth_m), ((DO_0.1*vol_depths$Vol_m3[1])+(DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(hypo_top_depth_m == 1.6, ((DO_1.6*vol_depths$Vol_m3[2])+(DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                 ifelse(hypo_top_depth_m == 3.8, ((DO_3.8*vol_depths$Vol_m3[3])+(DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                        ifelse(hypo_top_depth_m == 5, ((DO_5.0*vol_depths$Vol_m3[4])+(DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                               ifelse(hypo_top_depth_m == 6.2, ((DO_6.2*vol_depths$Vol_m3[5])+(DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                      ifelse(hypo_top_depth_m == 8, ((DO_8.0*vol_depths$Vol_m3[6])+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                             ifelse(hypo_top_depth_m == 9, DO_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Fill in some 2021 data where we're missing 'in between' depths
do_mgL <- do_mgL %>% 
  mutate(epi_DO = ifelse(is.na(epi_DO) & epi_bottom_depth_m == 5, ((DO_0.1*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(DO_5.0*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_DO),
         hypo_DO = ifelse(is.na(hypo_DO) & hypo_top_depth_m == 6.2, ((DO_8.0*sum(vol_depths$Vol_m3[5]+vol_depths$Vol_m3[6]))+(DO_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),hypo_DO))

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
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(size=0.75)+
  geom_point(size=2)+
  xlab("")+
  ylab("Time since anoxia (d)")+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

ggsave("./Figs/Fig_S11_DaysAnoxia.jpg",width=7,height=4,units="in",dpi=320)

## Plot model results by oxic vs. anoxic waters in the hypolimnion
doc_model_oxy <- left_join(hypo_do_mgL,doc_processing,by="DateTime") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01")) %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  drop_na(anoxia)

## Add in VW Hypo DOC (mg/L)
doc_model_oxy <- left_join(doc_model_oxy,doc_mgL,by=c("DateTime","Depth"))

## Plot by hypo internal loading under oxic vs. anoxic waters
wilcox.test(mean_doc_hypo_process_g~anoxia,doc_model_oxy)

## Plot distributions - oxic vs. anoxic
## Test statistical significance of distributions from zero
t.test(doc_model_oxy$mean_doc_hypo_process_g[doc_model_oxy$anoxia==1],mu=0)
t.test(doc_model_oxy$mean_doc_hypo_process_g[doc_model_oxy$anoxia==0],mu=0)

## Test that the two distributions are different
t.test(doc_model_oxy$mean_doc_hypo_process_g[doc_model_oxy$anoxia==1],doc_model_oxy$mean_doc_hypo_process_g[doc_model_oxy$anoxia==0])

hypo_oxy_proc <- doc_model_oxy %>% 
  ggplot(mapping=aes(x=mean_doc_hypo_process_g/1000,y=as.factor(anoxia),fill=as.factor(anoxia)))+
  stat_density_ridges(quantile_lines=TRUE,quantiles=2)+
  geom_vline(mapping=aes(xintercept=0),linetype="dashed")+
  scale_fill_manual(values=c("#7EBDC2", "#9B9B9B"))+
  scale_y_discrete(breaks=c("0","1"),
                   labels=c("Oxic", "Anoxic"))+
  xlab(expression(paste("Hypo. Internal DOC (kg d"^-1*")")))+
  ylab("")+
  theme_ridges()+
  theme(legend.position = "none")

## Find mean under oxic conditions
doc_model_oxy %>% 
  group_by(anoxia) %>% 
  summarise_at(vars(DOC_mgL), funs(mean(., na.rm=TRUE)))

## Test that the two distributions are different
t.test(doc_model_oxy$DOC_mgL[doc_model_oxy$anoxia==1],doc_model_oxy$DOC_mgL[doc_model_oxy$anoxia==0])

hypo_oxy_conc <- doc_model_oxy %>% 
  ggplot(mapping=aes(x=DOC_mgL,y=as.factor(anoxia),fill=as.factor(anoxia)))+
  stat_density_ridges(quantile_lines=TRUE,quantiles=2)+
  geom_vline(mapping=aes(xintercept=2.64),linetype="dashed")+ #oxic
  annotate("text", x = 2.8, y=1.2, label = "Oxic",angle = 90,size=5)+
  geom_vline(mapping=aes(xintercept=2.31),linetype="dashed")+ #anoxic
  annotate("text", x = 2.4, y=3.1, label = "Anoxic",angle = 90,size=5)+
  scale_fill_manual(values=c("#7EBDC2", "#9B9B9B"))+
  scale_y_discrete(breaks=c("0","1"),
                   labels=c("Oxic", "Anoxic"))+
  xlab(expression(paste("Hypo. DOC (mg L"^-1*")")))+
  ylab("")+
  theme_ridges()+
  theme(legend.position = "none")

ggarrange(hypo_oxy_proc,hypo_oxy_conc,ncol=1,nrow=2, labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/FigS10_Hypo_Oxy_FC.jpg",width=9,height=9,units="in",dpi=320)

###############################################################################
## Load in Flora data - Chla and community analysis
## Data last downloaded on 12 Mar 2024
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/6/6b3151c0fdd913e02641363c2b00ae57" 
#infile1 <- paste0(getwd(),"/Data/FluoroProbe.csv")
#download.file(inUrl1,infile1,method="curl")

flora <- read.csv("./Data/FluoroProbe.csv", header=T) %>%
  select(Reservoir:Depth_m,TotalConc_ugL) %>%
  dplyr::filter(Reservoir=="FCR") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime > as.POSIXct("2017-01-01")) %>% 
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
  pivot_wider(names_from = new_depth, values_from = TotalConc_ugL, values_fill = NA, values_fn = mean, names_prefix = "Chla_")

chla_ugL <- left_join(chla_ugL, thermo, by="DateTime")

chla_ugL <- chla_ugL %>% 
  drop_na(thermo.depth) %>% 
  mutate(epi_Chla = ifelse(is.na(epi_bottom_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottom_depth_m == 0.1, Chla_0.1,
                                  ifelse(epi_bottom_depth_m == 1.6, (Chla_0.1*vol_depths$Vol_m3[1]+Chla_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottom_depth_m == 3.8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottom_depth_m == 5, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottom_depth_m == 6.2, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottom_depth_m == 8, ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_Chla = ifelse(is.na(hypo_top_depth_m), ((Chla_0.1*vol_depths$Vol_m3[1])+(Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((Chla_1.6*vol_depths$Vol_m3[2])+(Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((Chla_3.8*vol_depths$Vol_m3[3])+(Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((Chla_5.0*vol_depths$Vol_m3[4])+(Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((Chla_6.2*vol_depths$Vol_m3[5])+(Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((Chla_8.0*vol_depths$Vol_m3[6])+(Chla_9.0*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, Chla_9.0, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Save formatted data to compare with EEMs
write_csv(chla_ugL, "./Data/Formated_Chla_ugL.csv")

## Plot
chla_plot <- chla_ugL %>%  
  drop_na(epi_Chla,hypo_Chla) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Chla,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(VW~Phyto~(~mu*g~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

# Format for ARIMA modeling
final_chla_ugL <- chla_ugL %>% 
  select(DateTime,epi_Chla,hypo_Chla) %>% 
  pivot_longer(!DateTime, names_to = "Depth", values_to = "VW_Chla_ugL") %>% 
  mutate(Depth = ifelse(Depth == "epi_Chla", "Epi",
                        ifelse(Depth == "hypo_Chla", "Hypo", NA)))

###############################################################################
## Explore adding in nutrient data (N and P) - from the FCR dam
chem <- read.csv("./Data/chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR",Site=="50") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))%>% 
  filter(DateTime >= as.POSIXct("2017-01-01") & DateTime < as.POSIXct("2022-01-01")) %>% 
  mutate(DIN_ugL = NH4_ugL+NO3NO2_ugL) %>% 
  select(DateTime,Depth_m,TN_ugL,TP_ugL,DIN_ugL,SRP_ugL)

## Calculate epi and hypo TP, TN, DIN, and SRP numbers
## TN
TN_ugL <- chem %>% 
  select(DateTime,Depth_m,TN_ugL) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5,6.2,8,9)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = TN_ugL, values_fill = NA, values_fn = mean, names_prefix = "TN_") %>% 
  filter(!if_all(c(TN_0.1, TN_1.6), is.na)) %>% 
  filter(!if_all(c(TN_3.8,TN_5,TN_6.2), is.na)) %>% 
  filter(!if_all(c(TN_8,TN_9),is.na))

TN_ugL <- left_join(TN_ugL, thermo, by="DateTime")

TN_ugL <- TN_ugL %>% 
  drop_na(thermo.depth) %>% 
  mutate(epi_TN = ifelse(is.na(epi_bottom_depth_m), ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                           ifelse(epi_bottom_depth_m == 0.1, TN_0.1,
                                  ifelse(epi_bottom_depth_m == 1.6, (TN_0.1*vol_depths$Vol_m3[1]+TN_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                         ifelse(epi_bottom_depth_m == 3.8, ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                                ifelse(epi_bottom_depth_m == 5, ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                       ifelse(epi_bottom_depth_m == 6.2, ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                              ifelse(epi_bottom_depth_m == 8, ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_TN = ifelse(is.na(hypo_top_depth_m), ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                            ifelse(hypo_top_depth_m == 1.6, ((TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                   ifelse(hypo_top_depth_m == 3.8, ((TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                          ifelse(hypo_top_depth_m == 5, ((TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                                 ifelse(hypo_top_depth_m == 6.2, ((TN_6.2*vol_depths$Vol_m3[5])+(TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                        ifelse(hypo_top_depth_m == 8, ((TN_8*vol_depths$Vol_m3[6])+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                               ifelse(hypo_top_depth_m == 9, TN_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Fill in some 2021 data where we're missing 'in between' depths
TN_ugL <- TN_ugL %>% 
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 1.6, TN_1.6, epi_TN),
         hypo_TN = ifelse(is.na(hypo_TN) & hypo_top_depth_m == 3.8, ((TN_5*(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4]+vol_depths$Vol_m3[5]))+(TN_9*(vol_depths$Vol_m3[6]+vol_depths$Vol_m3[7])))/sum(vol_depths$Vol_m3[3:7]),hypo_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 1.6, TN_0.1, epi_TN),
         epi_hypo = ifelse(is.na(hypo_TN) & hypo_top_depth_m == 3.8, ((TN_3.8*vol_depths$Vol_m3[3])+(TN_5*vol_depths$Vol_m3[4])+(TN_6.2*vol_depths$Vol_m3[5])+(TN_8*(vol_depths$Vol_m3[6]+vol_depths$Vol_m3[7])))/sum(vol_depths$Vol_m3[3:7]),hypo_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 3.8, ((TN_0.1*(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TN_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),epi_TN),
         hypo_TN = ifelse(is.na(hypo_TN) & hypo_top_depth_m == 5, ((TN_5*(vol_depths$Vol_m3[4]+vol_depths$Vol_m3[5]))+(TN_9*(vol_depths$Vol_m3[7]+vol_depths$Vol_m3[6])))/sum(vol_depths$Vol_m3[4:7]),hypo_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 3.8, TN_0.1 ,epi_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 3.8, TN_1.6 ,epi_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 5, ((TN_0.1*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TN_5*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_TN),
         hypo_TN = ifelse(is.na(hypo_TN) & hypo_top_depth_m == 6.2, ((TN_8*sum(vol_depths$Vol_m3[5]+vol_depths$Vol_m3[6]))+(TN_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),hypo_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 5, ((TN_1.6*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TN_5*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_TN),
         hypo_TN = ifelse(is.na(hypo_TN) & hypo_top_depth_m == 6.2, TN_9, hypo_TN)) %>% 
  
  mutate(epi_TN = ifelse(is.na(epi_TN) & epi_bottom_depth_m == 5, ((TN_0.1*vol_depths$Vol_m3[1])+(TN_1.6*vol_depths$Vol_m3[2])+(TN_3.8*(vol_depths$Vol_m3[3])+vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),epi_TN))

## TP
TP_ugL <- chem %>% 
  select(DateTime,Depth_m,TP_ugL) %>% 
  filter(Depth_m %in% c(0.1,1.6,3.8,5,6.2,8,9)) %>% 
  drop_na() %>% 
  pivot_wider(names_from = Depth_m, values_from = TP_ugL, values_fill = NA, values_fn = mean, names_prefix = "TP_") %>% 
  filter(!if_all(c(TP_0.1, TP_1.6), is.na)) %>% 
  filter(!if_all(c(TP_3.8,TP_5,TP_6.2), is.na)) %>% 
  filter(!if_all(c(TP_8,TP_9),is.na))

TP_ugL <- left_join(TP_ugL, thermo, by="DateTime")

TP_ugL <- TP_ugL %>% 
  drop_na(thermo.depth) %>% 
  mutate(epi_TP = ifelse(is.na(epi_bottom_depth_m), ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                         ifelse(epi_bottom_depth_m == 0.1, TP_0.1,
                                ifelse(epi_bottom_depth_m == 1.6, (TP_0.1*vol_depths$Vol_m3[1]+TP_1.6*vol_depths$Vol_m3[2])/sum(vol_depths$Vol_m3[1:2]),
                                       ifelse(epi_bottom_depth_m == 3.8, ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),
                                              ifelse(epi_bottom_depth_m == 5, ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),
                                                     ifelse(epi_bottom_depth_m == 6.2, ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5]))/sum(vol_depths$Vol_m3[1:5]),
                                                            ifelse(epi_bottom_depth_m == 8, ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6]))/sum(vol_depths$Vol_m3[1:6]), NA)))))))) %>% 
  mutate(hypo_TP = ifelse(is.na(hypo_top_depth_m), ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[1:7]),
                          ifelse(hypo_top_depth_m == 1.6, ((TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[2:7]),
                                 ifelse(hypo_top_depth_m == 3.8, ((TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[3:7]),
                                        ifelse(hypo_top_depth_m == 5, ((TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[4:7]),
                                               ifelse(hypo_top_depth_m == 6.2, ((TP_6.2*vol_depths$Vol_m3[5])+(TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),
                                                      ifelse(hypo_top_depth_m == 8, ((TP_8*vol_depths$Vol_m3[6])+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[6:7]),
                                                             ifelse(hypo_top_depth_m == 9, TP_9, NA)))))))) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

## Fill in some 2021 data where we're missing 'in between' depths
TP_ugL <- TP_ugL %>% 
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 1.6, TP_1.6, epi_TP),
         hypo_TP = ifelse(is.na(hypo_TP) & hypo_top_depth_m == 3.8, ((TP_5*(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4]+vol_depths$Vol_m3[5]))+(TP_9*(vol_depths$Vol_m3[6]+vol_depths$Vol_m3[7])))/sum(vol_depths$Vol_m3[3:7]),hypo_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 1.6, TP_0.1, epi_TP),
         epi_hypo = ifelse(is.na(hypo_TP) & hypo_top_depth_m == 3.8, ((TP_3.8*vol_depths$Vol_m3[3])+(TP_5*vol_depths$Vol_m3[4])+(TP_6.2*vol_depths$Vol_m3[5])+(TP_8*(vol_depths$Vol_m3[6]+vol_depths$Vol_m3[7])))/sum(vol_depths$Vol_m3[3:7]),hypo_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 3.8, ((TP_0.1*(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TP_3.8*vol_depths$Vol_m3[3]))/sum(vol_depths$Vol_m3[1:3]),epi_TP),
         hypo_TP = ifelse(is.na(hypo_TP) & hypo_top_depth_m == 5, ((TP_5*(vol_depths$Vol_m3[4]+vol_depths$Vol_m3[5]))+(TP_9*(vol_depths$Vol_m3[7]+vol_depths$Vol_m3[6])))/sum(vol_depths$Vol_m3[4:7]),hypo_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 3.8, TP_0.1 ,epi_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 3.8, TP_1.6 ,epi_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 5, ((TP_0.1*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TP_5*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_TP),
         hypo_TP = ifelse(is.na(hypo_TP) & hypo_top_depth_m == 6.2, ((TP_8*sum(vol_depths$Vol_m3[5]+vol_depths$Vol_m3[6]))+(TP_9*vol_depths$Vol_m3[7]))/sum(vol_depths$Vol_m3[5:7]),hypo_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 5, ((TP_1.6*sum(vol_depths$Vol_m3[1]+vol_depths$Vol_m3[2]))+(TP_5*sum(vol_depths$Vol_m3[3]+vol_depths$Vol_m3[4])))/sum(vol_depths$Vol_m3[1:4]),epi_TP),
         hypo_TP = ifelse(is.na(hypo_TP) & hypo_top_depth_m == 6.2, TP_9, hypo_TP)) %>% 
  
  mutate(epi_TP = ifelse(is.na(epi_TP) & epi_bottom_depth_m == 5, ((TP_0.1*vol_depths$Vol_m3[1])+(TP_1.6*vol_depths$Vol_m3[2])+(TP_3.8*(vol_depths$Vol_m3[3])+vol_depths$Vol_m3[4]))/sum(vol_depths$Vol_m3[1:4]),epi_TP))

## Combine all Chem data
chem_epi_hypo <- plyr::join_all(list(TN_ugL,TP_ugL),by=c("DateTime","thermo.depth","epi_bottom_depth_m","hypo_top_depth_m"),type="left") %>% 
  select(DateTime,epi_TN,hypo_TN,epi_TP,hypo_TP)

## Plot
tn_plot <- chem_epi_hypo %>%  
  drop_na(epi_TN,hypo_TN) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_TN,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_TN,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_TN,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_TN,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(VW~TN~(~mu*g~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

TP_plot <- chem_epi_hypo %>%  
  drop_na(epi_TP,hypo_TP) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_TP,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_TP,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_TP,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_TP,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(VW~TP~(~mu*g~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(tn_plot,TP_plot,ncol=1,nrow=2, labels = c("A.", "B."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig_S8_FCR50_TotalNuts.jpg",width=10,height=6,units="in",dpi=320)

###############################################################################
## Load in daily inflow - from Eco_DOC_rlnorm_FC.R
## Now includes both TB and FC input
inflow_daily <- read.csv("./Data/inflow_daily.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

inflow_daily_fc <- read.csv("./Data/inflow_daily_fc.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

# Plot daily inflow for the study period
inflow_plot <- 
  ggplot(mapping=aes(x=DateTime,y=mean_flow_cms))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(inflow_daily %>% na.omit(mean_flow_cms),mapping=aes(x=DateTime,y=mean_flow_cms,color="Weir"),size=1)+
  geom_point(inflow_daily %>% na.omit(mean_flow_cms),mapping=aes(x=DateTime,y=mean_flow_cms,color="Weir"),size=2)+
  geom_line(inflow_daily_fc %>% drop_na(est_flow_cms),mapping=aes(x=DateTime,y=est_flow_cms,color="FC"),size=1)+
  scale_color_manual(breaks=c('Weir','FC'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Weir','FC'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Inflow~(~m^3~s^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

## Plot seasonal discharge by year
func_order <- c("Pre-Stratification","Stratified Period","Post-Stratification")

season_inflow <- inflow_daily %>% 
  mutate(doy = yday(DateTime),
         year = year(DateTime)) %>% 
  mutate(season = ifelse(doy>=1 & doy<122, "Pre-Stratification", 
                         ifelse(doy>=122 & doy<=320, "Stratified Period","Post-Stratification"))) %>% 
  ggplot()+
  geom_boxplot(aes(x=as.factor(year),y=mean_flow_cms,fill=factor(season,func_order)))+
  scale_fill_manual(values=c("#9B9B9B","#F0B670","#7EBDC2"))+
  ylab(expression(paste("Primary Inflow (m"^3*" s"^-1*")")))+
  xlab("")+
  theme_bw(base_size = 15) +
  theme(legend.title=element_blank())+
  theme(legend.position="top")

## Calculate mean_flow_cms stratified flow from TB and FC
## For R2R
inflow_daily %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<320) %>% 
  select(mean_flow_cms) %>% 
  summarize_all(mean_flow_cms,na.rm=TRUE)

inflow_daily %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<320) %>% 
  select(mean_flow_cms) %>% 
  summarize_all(sd,na.rm=TRUE)

inflow_daily_fc %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<320) %>% 
  select(est_flow_cms) %>% 
  summarize_all(mean_flow_cms,na.rm=TRUE)

inflow_daily_fc %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<320) %>% 
  select(est_flow_cms) %>% 
  summarize_all(sd,na.rm=TRUE)

# Format for ARIMA modeling - sum Weir inflow and FC inflow = daily inflow for ARIMA modeling
final_inflow_m3s <- left_join(inflow_daily,inflow_daily_fc,by="DateTime") %>% 
  mutate(Inflow_m3s = mean_flow_cms+est_flow_cms) %>% 
  select(DateTime,Inflow_m3s)

# Estimate percentage of of TB vs. FC flow
left_join(inflow_daily,inflow_daily_fc,by="DateTime") %>% 
  mutate(Inflow_m3s = mean_flow_cms+est_flow_cms) %>% 
  mutate(perc_tb = mean_flow_cms/Inflow_m3s*100,
         perc_fc = est_flow_cms/Inflow_m3s*100) %>% 
  summarise(mean_flow_cms_tb = mean_flow_cms(perc_tb,na.rm=TRUE),
            mean_flow_cms_fc = mean_flow_cms(perc_fc,na.rm=TRUE),
            sd_tb = sd(perc_tb,na.rm=TRUE),
            sd_fc = sd(perc_fc,na.rm=TRUE))

# Calculate residence time assuming full pond
wtr_d <- final_inflow_m3s %>% 
  mutate(wtr_d = (310000)/Inflow_m3s/60/60/24,
         month = month(DateTime),
         year = year(DateTime))

wtr_d %>% 
  na.omit(wtr_d) %>% 
  summarise(med = median(wtr_d),
            sd = sd(wtr_d)) %>% 
  mutate(se = sd/length(wtr_d$wtr_d))

## Plot WRT annually for each summer stratified period
annual_wrt <- ggplot(wtr_d,mapping=aes(x=as.character(year),y=log10(wtr_d)))+
  geom_boxplot(size=0.8,alpha=0.5)+
  xlab("")+
  ylab("logt10(WRT, d)")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(season_inflow, annual_wrt, ncol=1, nrow=2, labels=c("A.","B."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig_S9_SeasonalInflow_WRT.png",dpi=800,width=8,height=8)

###############################################################################
## Plot Temp, DO, Chla, and Total Inflow for main MS
ggarrange(temp_plot,do_plot,chla_plot,inflow_plot,ncol=1,nrow=4, labels = c("A.", "B.", "C.", "D."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig6_EnvParameters_FC.jpg",width=10,height=12,units="in",dpi=320)
  
###############################################################################
## Load in met data - shortwave radiation and rainfall
# Downloaded on 12 Mar 2024
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/389/6/a5524c686e2154ec0fd0459d46a7d1eb"
#infile1 <- paste0(getwd(),"/Data/FCR_Met_final_2015_2021.csv")
#download.file(inUrl1,infile1,method="curl")

## Load in Met data
met <- read.csv("./Data/FCR_Met_final_2015_2021.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) %>% 
  filter(DateTime >= as.POSIXct("2017-01-01 00:00:00"))

## Average to daily
met_daily <- met %>% 
  mutate(DateTime = format(as.POSIXct(DateTime, "%Y-%m-%d"),"%Y-%m-%d", tz="EST")) %>% 
  mutate(DateTime = as.POSIXct(DateTime, "%Y-%m-%d", tz = "EST")) %>% 
  group_by(DateTime) %>% 
  dplyr::summarise(rain_tot_mm = sum(Rain_Total_mm, na.rm = TRUE),
            ShortwaveRadiationUp_Average_W_m2 = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE),
            Temp_Average_C = mean(AirTemp_Average_C,na.rm=TRUE),
            BP_Average_kPa = mean(BP_Average_kPa,na.rm=TRUE),
            InfraredRadiationDown_Average_W_m2 = mean(InfraredRadiationDown_Average_W_m2,na.rm=TRUE),
            InfraredRadiationUp_Average_W_m2 = mean(InfraredRadiationUp_Average_W_m2,na.rm=TRUE),
            Albedo_Average_W_m2 = mean(Albedo_Average_W_m2,na.rm=TRUE)) %>% 
  filter(ShortwaveRadiationUp_Average_W_m2 < 450) # Remove 2018-04-05 - maintenance and only measured afternoon values!

## Then plot total rainfall and shortwave radiation with thermocline depth and DO
## Plot timeseries of thermocline depth for SI
vw_do_plot <- do_mgL %>%  
  drop_na(epi_DO,hypo_DO) %>% 
  ggplot()+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_DO,color="Epi"),size=2)+
  geom_line(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_DO,color="Hypo"),size=2)+
  scale_color_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  scale_fill_manual(breaks=c('Epi','Hypo'),values=c("#7EBDC2","#393E41"))+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(VW~DO~(mg~L^-1)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

thermo_plot <- thermo %>% 
  drop_na(thermo.depth) %>% 
  ggplot(mapping=aes(x=DateTime,y=-thermo.depth))+
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_hline(yintercept = -0.1, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -1.6, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -3.8, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -5, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -6.2, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -8, linetype="dotted",color="darkgrey")+
  geom_hline(yintercept = -9, linetype="dotted",color="darkgrey")+
  geom_line(size=1)+
  geom_point(size=2)+
  ylab("Thermo. depth (m)")+
  xlab("")+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  theme_classic(base_size = 15)

rain_plot <- met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=rain_tot_mm))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(Total~Rainfall~(mm)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

sw_plot <- met_daily %>% 
  ggplot(mapping=aes(x=DateTime,y=ShortwaveRadiationUp_Average_W_m2))+ 
  geom_vline(xintercept = as.POSIXct("2017-10-25"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2017-05-01"), xmax = as.POSIXct("2017-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2018-10-21"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2018-05-01"), xmax = as.POSIXct("2018-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2019-05-01"), xmax = as.POSIXct("2019-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2020-11-01"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2020-05-01"), xmax = as.POSIXct("2020-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_vline(xintercept = as.POSIXct("2021-11-03"),linetype="dashed",color="darkgrey")+
  annotate("rect", xmin = as.POSIXct("2021-05-01"), xmax = as.POSIXct("2021-11-15"), ymin = -Inf, ymax = Inf,alpha = .3,fill = "darkgrey")+
  geom_line(size=1)+
  xlim(as.POSIXct("2017-01-01"),as.POSIXct("2021-12-31"))+
  xlab("") + 
  ylab(expression(SW~Rad~(W~m^2)))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(vw_do_plot,thermo_plot,rain_plot,sw_plot,ncol=1,nrow=4,common.legend = TRUE, labels = c("A.","B.","C.","D."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig_S7_WaterCol_MetParameters.jpg",width=10,height=12,units="in",dpi=320)

###############################################################################
## Format and add thermocline depth as a potential predictor variable, too
final_thermo <- thermo %>% 
  select(DateTime, thermo.depth)

###############################################################################
###############################################################################
## Organize data for ARIMA modeling - following EddyFlux, 3_Rev_EnvAnalysis.R
# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, and SW Radiation for Epi and Hypo
arima_epi <- plyr::join_all(list(all_doc_mgL,final_temp_c,final_do_pSat,final_chla_ugL),by=c("DateTime","Depth"),type="left") %>% 
  select(-month,-year.x,-year.y)

# Include: DOC data, Temp, DO, Flora, Inflow, Rainfall, SW Radiation, CO2, and CH4 for Epi and Hypo
arima_hypo <- plyr::join_all(list(all_doc_mgL,final_temp_c,final_do_pSat,final_chla_ugL),by=c("DateTime","Depth"),type="left") %>% 
  select(-month,-year.x,-year.y)

## Select time points where we have DOC concentrations
arima_epi <- plyr::join_all(list(arima_epi,met_daily,wtr_d),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Epi") %>% 
  select(-month,-year)

arima_hypo <- plyr::join_all(list(arima_hypo,met_daily,wtr_d),by="DateTime",type="left") %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"),Depth == "Hypo") %>% 
  select(-month,-year)

## Add in days since anoxia for Hypo
arima_hypo <- left_join(arima_hypo,hypo_do_mgL,by=c("DateTime","Depth")) %>% 
  dplyr::select(DateTime,Depth,DOC_mgL,DOC_processing_mgL,VW_Temp_C,VW_DO_mgL.x,VW_Chla_ugL,rain_tot_mm,
                Inflow_m3s,anoxia_time_d) %>% 
  dplyr::rename(VW_DO_pSat = VW_DO_mgL.x)

## Add in thermocline depth information
arima_epi <- left_join(arima_epi,final_thermo,by="DateTime")

arima_epi <- arima_epi %>% 
  dplyr::rename(VW_DO_pSat = VW_DO_mgL)

arima_hypo <- left_join(arima_hypo,final_thermo,by="DateTime")

## Add Epi Chla to Hypo (to represent sinking phytos)
epi_chla <- arima_epi %>% 
  select(DateTime,VW_Chla_ugL) %>% 
  rename(Epi_Chla_ugL = VW_Chla_ugL)

arima_hypo <- left_join(arima_hypo,epi_chla,by="DateTime")

## Add nutrient data, too
arima_epi <- left_join(arima_epi,chem_epi_hypo,by="DateTime") %>% 
  select(-hypo_TN,-hypo_TP)

arima_hypo <- left_join(arima_hypo,chem_epi_hypo,by="DateTime") %>% 
  select(-epi_TN,-epi_TP)

## Constrain to stratified time periods
arima_epi <- arima_epi %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  select(DateTime,DOC_mgL,DOC_processing_mgL,VW_Temp_C,VW_DO_pSat,VW_Chla_ugL,rain_tot_mm,ShortwaveRadiationUp_Average_W_m2,Inflow_m3s,thermo.depth,epi_TN,epi_TP)

arima_hypo <- arima_hypo %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  select(DateTime,DOC_mgL,DOC_processing_mgL,VW_Temp_C,VW_DO_pSat,VW_Chla_ugL,rain_tot_mm,Inflow_m3s,anoxia_time_d,thermo.depth,Epi_Chla_ugL,hypo_TN,hypo_TP)

## Calculate stats for env parameters - limited to summer stratified period (May-Oct)
epi_stats <- arima_epi %>% 
  mutate(month = month(DateTime)) %>% 
  select(VW_Temp_C,VW_DO_pSat,VW_Chla_ugL,epi_TN,epi_TP) %>% 
  summarise_all(list(min,max,median,mean,sd),na.rm=TRUE)%>% 
  pivot_longer(cols = contains("fn"),
               names_to = "func") %>% 
  mutate(stat = ifelse(grepl("fn1",func),"min",
                       ifelse(grepl("fn2",func),"max",
                              ifelse(grepl("fn3",func),"median",
                                     ifelse(grepl("fn4",func),"mean",
                                            ifelse(grepl("fn5",func),"sd",NA))))))

epi_stats$func <- substr(epi_stats$func,1,nchar(epi_stats$func)-4)

epi_stats <- epi_stats %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(year = "all")

hypo_stats <- arima_hypo %>% 
  mutate(month = month(DateTime)) %>% 
  select(VW_Temp_C,VW_DO_pSat,VW_Chla_ugL,hypo_TN,hypo_TP,rain_tot_mm,ShortwaveRadiationUp_Average_W_m2,Inflow_m3s,wtr_d) %>% 
  summarise_all(list(min,max,median,mean,sd),na.rm=TRUE) %>% 
  pivot_longer(cols = contains("fn"),
               names_to = "func") %>% 
  mutate(stat = ifelse(grepl("fn1",func),"min",
                       ifelse(grepl("fn2",func),"max",
                              ifelse(grepl("fn3",func),"median",
                                     ifelse(grepl("fn4",func),"mean",
                                            ifelse(grepl("fn5",func),"sd",NA))))))

hypo_stats$func <- substr(hypo_stats$func,1,nchar(hypo_stats$func)-4)

hypo_stats <- hypo_stats %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(year = "all")

## Calculate DO_mgL stats, too
do_mgL %>%  mutate(month = month(DateTime)) %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320) %>% 
  select(epi_DO,hypo_DO) %>% 
  summarise_all(list(min,max,median,mean,sd),na.rm=TRUE)

###############################################################################
## Check correlations among environmental variables - what needs to be removed?
# Epi
epi_cor = as.data.frame(cor(arima_epi[,2:12],use = "complete.obs"),method=c("pearson"))
write_csv(epi_cor, "./Figs/epi_cor_fc.csv")

chart.Correlation(arima_epi[,4:12],histogram = TRUE,method=c("pearson"))
# Remove Epi TN
arima_epi <- arima_epi %>% 
  select(-epi_TN)

# Hypo
hypo_cor = as.data.frame(cor(arima_hypo[,2:13],use = "complete.obs"),method=c("pearson"))
write_csv(hypo_cor, "./Figs/hypo_cor_fc.csv")

chart.Correlation(arima_hypo[,2:13],histogram = TRUE,method=c("pearson"))

###############################################################################
## Check for skewness following MEL script!
Math.cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

# Epi
for (i in 2:12){
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
# Nothing: DO, SW Radiation, Thermo, rain, TP
# Log: DOC, Chla, Inflow
# Cube root: DOC processing
# sqrt: Temp

# Transform and scale data
arima_epi_scale <- arima_epi %>% 
  mutate(DOC_mgL = log(DOC_mgL),
         VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s),
         DOC_processing_mgL = Math.cbrt(DOC_processing_mgL),
         VW_Temp_C = sqrt(VW_Temp_C))

arima_epi_scale[,2:11] <- scale(arima_epi_scale[,2:11])

# Hypo
for (i in 2:13){
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
# Nothing: DOC, DO, SW Radiation, Thermo, Anoxia_time, rain
# Log: Chla_hypo, Inflow, Chla_epi, TN, TP
# Cube root: DOC_processing
# sqrt: Temp

# Transform and scale data
arima_hypo_scale <- arima_hypo %>% 
  mutate(DOC_processing_mgL = Math.cbrt(DOC_processing_mgL),
         VW_Chla_ugL = log(VW_Chla_ugL),
         Inflow_m3s = log(Inflow_m3s),
         VW_Temp_C = sqrt(VW_Temp_C),
         Epi_Chla_ugL = log(Epi_Chla_ugL),
         hypo_TP = log(hypo_TP),
         hypo_TN = log(hypo_TN))

arima_hypo_scale[,2:13] <- scale(arima_hypo_scale[,2:13])

###############################################################################

## ARIMA modeling - DOC concentrations!
# Following MEL code : )
# THINGS TO CHANGE: 'cols' (change to the environmental variables); best fit!

###############################################################################
## Check Assumption about number of AR terms
png("./Figs/Fig_SX_PACF.png", width = 750, height = 650)
par(mfrow=c(2,2))

## [DOC] Epi - limit to AR(2)
pacf(arima_epi$DOC_mgL,xlim=c(1,20),na.action=na.pass,main="Epi. DOC (mg per L)")

# [DOC] Hypo - limit to AR(2)
pacf(arima_hypo$DOC_mgL,xlim=c(1,20),na.action=na.pass,main="Hypo. DOC (mg per L)")

# DOC processing epi - limit to AR(1)
pacf(arima_epi$DOC_processing_mgL,xlim=c(1,20),na.action=na.pass,main="Epi. DOC Processing (kg per d)")

# DOC processing hypo - limit to AR(1)
pacf(arima_hypo$DOC_processing_mgL,xlim=c(1,20),na.action=na.pass,main="Hypo. DOC Processing (kg per d)")

dev.off()

###############################################################################
# Epi
colnames(arima_epi_scale)

cols <- c(4:11) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_epi_scale[,2] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE}) # Constrain to 1 AR term
    
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

best

best.vars <- colnames(arima_epi_scale)[combn(cols,4)[,24]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,24] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(arima_epi_scale[,3]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_epi_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_epi_scale[,3])-unlist(fitted(best.fit))), na.rm = TRUE)

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
  print(forecast::accuracy(good.fit.1))
  
  
}

###############################################################################
# Hypo
colnames(arima_hypo_scale)

cols <- c(4:13) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,2] # UPDATE THIS TO DOC CONCENTRATION

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,my.combn[,j]]),max.p = 2, max.P = 2), error = function(e) { skip_to_next <<- TRUE})
    
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

fit <- auto.arima(y, max.p = 2, max.P = 2)
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

best

best.vars <- colnames(arima_hypo_scale)[combn(cols,4)[,37]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,37] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 2, max.P = 2)
best.fit
hist(resid(best.fit))
forecast::accuracy(best.fit)
hist(unlist(arima_hypo_scale[,3]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,3])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 2, max.P = 2)
  print(good.fit.1)
  print(forecast::accuracy(good.fit.1))
  
  
}

###############################################################################
## Epi DOC Processing
colnames(arima_epi_scale)

cols <- c(4:11) # UPDATE THIS TO THE ENV. VARIABLES - include TP only; remove WRT!
sub.final <- NULL
final <- NULL

y <- arima_epi_scale[,3] # UPDATE THIS TO DOC PROCESSING

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

best

best.vars <- colnames(arima_epi_scale)[combn(cols,4)[,45]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,4)[,45] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_epi_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
forecast::accuracy(best.fit)
hist(unlist(arima_epi_scale[,4]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_epi_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_epi_scale[,4])-unlist(fitted(best.fit))), na.rm = TRUE)

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
  print(forecast::accuracy(good.fit.1))
  
  
}

###############################################################################
## Hypo DOC processing
colnames(arima_hypo_scale)

cols <- c(4:13) # UPDATE THIS TO THE ENV. VARIABLES
sub.final <- NULL
final <- NULL

y <- arima_hypo_scale[,3] # UPDATE THIS TO DOC PROCESSING

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

best

best.vars <- colnames(arima_hypo_scale)[combn(cols,3)[,56]] # UPDATE THIS FOLLOWING 'BEST'
best.vars.cols <- combn(cols,3)[,56] # UPDATE THIS FOLLOWING 'BEST'

best.fit <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
forecast::accuracy(best.fit)
hist(unlist(arima_hypo_scale[,4]))
plot_fit <- as.numeric(fitted(best.fit))
plot_x <- as.numeric(unlist(arima_hypo_scale[,1]))
plot(plot_x,plot_fit)
abline(a = 0, b = 1)
median((unlist(arima_hypo_scale[,4])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc >= as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){0.16
  good.vars.1 <- colnames(arima_hypo_scale)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(arima_hypo_scale[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(forecast::accuracy(good.fit.1))
  
  
}

###############################################################################
