###############################################################################

### Script to explore use of optical data (EEMs) to support authochthonous
### production in FCR - for Summer 2019 ONLY
### 2 Jul 2025, A. Hounshell

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
## Download optical data from EDI
## For Summer 2019

# FCR Optical data - https://portal.edirepository.org/nis/mapbrowse?packageid=edi.841.1
# Last Downloaded: 02 Jul 2025
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/841/1/f63272976bcd151f8e879cbd14d9a9ce"
#infile1 <- paste0(getwd(),"/Data/OpticalData.csv")
#download.file(inUrl1,infile1,method="curl")

## Select for FCR only and weir (100), wetlands (200), and dam (50)
## For Sta 50 - only select surface (0.1 m) and bottom (9 m)
eems_data <- read.csv("./Data/OpticalData.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%m/%d/%Y", tz="EST"))) %>% 
  filter(Reservoir == "FCR" & Site %in% c(100,200,50) & Depth_m %in% c(0.1,9) & DateTime<as.POSIXct("2020-01-01"))

## Combine replicates, when available
avg_eems <- eems_data %>% 
  select(-Reservoir) %>% 
  group_by(Site,DateTime,Depth_m) %>% 
  summarise_all(list(mean, sd),na.rm=TRUE)

## Define groups
avg_eems <- avg_eems %>% 
  mutate(loc = ifelse(Site == 50 & Depth_m == 0.1, "Surface",
                      ifelse(Site == 50 & Depth_m == 9, "Bottom",
                             ifelse(Site == 100, "TB",
                                    ifelse(Site == 200, "FC", NA))))) %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320)

###############################################################################
## Load in Chla for visualization w/ EEMs
chla_ugL <- read.csv("./Data/Formated_Chla_ugL.csv",header=T) %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
  filter(DateTime>as.POSIXct("2019-01-01") & DateTime<as.POSIXct("2019-12-01")) %>% 
  mutate(doy = yday(DateTime)) %>% 
  filter(doy>=122 & doy<=320)

###############################################################################
## Plot - surface EEMs
peak_a <- avg_eems %>% 
  drop_na(A_fn1) %>% 
  ggplot(mapping=aes(x=DateTime,y=A_fn1,color=loc))+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = Inf,ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=A_fn1-A_fn2, ymax=A_fn1+A_fn2))+
  geom_line(size=1)+
  scale_color_manual(breaks=c('Surface','Bottom','TB','FC'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Epi","Hypo","TB","FC"))+
  ylab("Peak A (Allo., RFU)")+
  xlab("2019")+
  scale_x_continuous(breaks=c(as.POSIXct("2019-05-01"),as.POSIXct("2019-06-01"),as.POSIXct("2019-07-01"),as.POSIXct("2019-08-01"),as.POSIXct("2019-09-01"),as.POSIXct("2019-10-01"),as.POSIXct("2019-11-01")),
                     limits = c(as.POSIXct("2019-05-01"),as.POSIXct("2019-11-15")),
                     labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

peak_t <- avg_eems %>% 
  drop_na(T_fn1) %>% 
  ggplot(mapping=aes(x=DateTime,y=T_fn1,color=loc))+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = Inf,ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=T_fn1-T_fn2, ymax=T_fn1+T_fn2))+
  geom_line(size=1)+
  scale_color_manual(breaks=c('Surface','Bottom','TB','FC'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Epi","Hypo","TB","FC"))+
  ylab("Peak T (Auto., RFU)")+
  xlab("2019")+
  scale_x_continuous(breaks=c(as.POSIXct("2019-05-01"),as.POSIXct("2019-06-01"),as.POSIXct("2019-07-01"),as.POSIXct("2019-08-01"),as.POSIXct("2019-09-01"),as.POSIXct("2019-10-01"),as.POSIXct("2019-11-01")),
                     limits = c(as.POSIXct("2019-05-01"),as.POSIXct("2019-11-15")),
                     labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ratio <- avg_eems %>% 
  drop_na(A_T_fn1) %>% 
  ggplot(mapping=aes(x=DateTime,y=A_T_fn1,color=loc))+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = Inf,ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_hline(yintercept = 1,linetype="dotted",color="darkgrey")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=A_T_fn1-A_T_fn2, ymax=A_T_fn1+A_T_fn2))+
  geom_line(size=1)+
  scale_color_manual(breaks=c('Surface','Bottom','TB','FC'),values=c("#7EBDC2","#393E41","#F0B670","#E7804B"),labels=c("Epi","Hypo","TB","FC"))+
  ylab("Ratio A:T")+
  xlab("2019")+
  scale_x_continuous(breaks=c(as.POSIXct("2019-05-01"),as.POSIXct("2019-06-01"),as.POSIXct("2019-07-01"),as.POSIXct("2019-08-01"),as.POSIXct("2019-09-01"),as.POSIXct("2019-10-01"),as.POSIXct("2019-11-01")),
                     limits = c(as.POSIXct("2019-05-01"),as.POSIXct("2019-11-15")),
                     labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov"))+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

chla <- chla_ugL %>% 
  ggplot()+
  annotate(geom="rect",xmin = as.POSIXct("2019-06-03"), xmax = as.POSIXct("2019-06-17"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-07-08"),xmax = as.POSIXct("2019-07-19"), ymin=-Inf, ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-08-05"),xmax = as.POSIXct("2019-08-19"),ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on
  annotate(geom="rect",xmin = as.POSIXct("2019-09-02"), xmax = Inf,ymin=-Inf,ymax=Inf,alpha=0.2)+ # Oxygen on (technically turned-off on 2019-12-01)
  geom_vline(xintercept = as.POSIXct("2019-10-23"),linetype="dashed",color="darkgrey")+
  geom_line(mapping=aes(x=DateTime,y=epi_Chla,color="Epi_Chla"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=epi_Chla,color="Epi_Chla"),size=3)+
  geom_line(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo_Chla"),size=1)+
  geom_point(mapping=aes(x=DateTime,y=hypo_Chla,color="Hypo_Chla"),size=3)+
  scale_color_manual(breaks=c('Epi_Chla','Hypo_Chla'),values=c("#7EBDC2","#393E41"),labels=c("Epi","Hypo"))+
  scale_x_continuous(breaks=c(as.POSIXct("2019-05-01"),as.POSIXct("2019-06-01"),as.POSIXct("2019-07-01"),as.POSIXct("2019-08-01"),as.POSIXct("2019-09-01"),as.POSIXct("2019-10-01"),as.POSIXct("2019-11-01")),
                     limits = c(as.POSIXct("2019-05-01"),as.POSIXct("2019-11-15")),
                     labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov"))+
  ylab(expression(VW~Phyto~(~mu*g~L^-1)))+
  xlab("2019")+
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(peak_a,peak_t,ratio,chla,common.legend = TRUE,ncol=1,nrow=4,labels = c("A.", "B.","C.","D."),
          font.label=list(face="plain",size=15))

ggsave("./Figs/Fig_S7_2019_EEMs.jpg",width=10,height=13,units="in",dpi=320)


