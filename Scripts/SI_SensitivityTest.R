### Test to see differences in hypo bypass following Carey et al. 2018
### Including: 94, 74, and 54% bypass

### 'Rough cut' sensitivity test
### 12 Apr. 2024, A. Hounshell

### Updated: 8 July 2024, A. Hounshell
## To create separate code for just the senstivity tests

###############################################################################
# Set working directory
wd <- getwd()
setwd(wd)

# Load libraries
pacman::p_load(tidyverse,ggplot2,ggpubr,rMR,lme4,PerformanceAnalytics,astsa,cowplot,lubridate,dplR,zoo,naniar,
               DescTools,MuMIn,rsq,Metrics,truncnorm)
###############################################################################
## Load in necessary model data - from Eco_DOC_rlnorm.R - if needed
doc_dt_epi <- read.csv("./Data/doc_dt_epi.csv")

doc_inflow_mass <- read.csv("./Data/doc_inflow_mass.csv")

doc_hypo_mass_outflow <- read.csv("./Data/doc_hypo_mass_outflow.csv")

doc_epi_mass_outflow <- read.csv("./Data/doc_epi_mass_outflow.csv")

doc_entr<- read.csv("./Data/doc_entr.csv")

epi_vol <- read.csv("./Data/epi_vol.csv")

doc_dt_hypo <- read.csv("./Data/doc_dt_hypo.csv")

hypo_vol <- read.csv("./Data/hypo_vol.csv")

###############################################################################
## Load in data from original model (p=74)
p74_doc_inputs_g <- read.csv("./Data/26Apr24_final_doc_inputs.csv") %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

###############################################################################
## Calculate model results for p=94
# Calculate processing for epi
doc_epi_process_g <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # DOC epi processing for each time point
doc_epi_process_mgL <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))

p = 0.94 # For sensitivity test

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

p94_final_doc_inputs_g <- cbind(DateTime,doc_inputs_g)

p94_final_doc_inputs_g <- na.omit(p94_final_doc_inputs_g)

p94_final_doc_inputs_g <- p94_final_doc_inputs_g %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"))

## Save final model output
write.csv(p94_final_doc_inputs_g, "./Data/p94_final_doc_inputs.csv",row.names=FALSE)

###############################################################################
## Calculate model results for p=54
# Calculate processing for epi
doc_epi_process_g <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime))) # DOC epi processing for each time point
doc_epi_process_mgL <- as.data.frame(matrix(data=NA, ncol=1000, nrow=length(doc_box_full$DateTime)))

p = 0.54 # For sensitivity test

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

p54_final_doc_inputs_g <- cbind(DateTime,doc_inputs_g)

p54_final_doc_inputs_g <- na.omit(p54_final_doc_inputs_g)

p54_final_doc_inputs_g <- p54_final_doc_inputs_g %>% 
  filter(DateTime >= as.POSIXct("2017-01-01"))

## Save final model output
write.csv(p54_final_doc_inputs_g, "./Data/p54_final_doc_inputs.csv",row.names=FALSE)

###############################################################################
## Plot to see comparisons
epi <- ggplot()+
  geom_line(p74_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="p74"))+
  geom_ribbon(p74_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="p74"),alpha=0.5,fill="grey")+
  geom_line(p94_final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="p94"))+
  #geom_ribbon(p94_final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="p94"),alpha=0.5)+
  geom_line(p54_final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_epi_process_mgL,color="p54"))+
  #geom_ribbon(p54_final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_epi_process_mgL-sd_doc_epi_process_mgL,ymax=mean_doc_epi_process_mgL+sd_doc_epi_process_mgL,fill="p54"),alpha=0.5)+
  ylab("Epi Internal DOC (mg/L)")+
  xlab("")+
  theme_bw(base_size = 15)+
  theme(legend.title=element_blank())

hypo <- ggplot()+
  geom_line(p74_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="p74"))+
  geom_ribbon(p74_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="p74"),alpha=0.5,fill="grey")+
  geom_line(p94_final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="p94"))+
  #geom_ribbon(p94_final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="p94"),alpha=0.5)+
  geom_line(p54_final_doc_inputs_g,mapping=aes(x=DateTime,y=mean_doc_hypo_process_mgL,color="p54"))+
  #geom_ribbon(p54_final_doc_inputs_g,mapping=aes(x=DateTime,ymin=mean_doc_hypo_process_mgL-sd_doc_hypo_process_mgL,ymax=mean_doc_hypo_process_mgL+sd_doc_hypo_process_mgL,fill="p54"),alpha=0.5)+
  ylab("Hypo Internal DOC (mg/L)")+
  xlab("")+
  theme_bw(base_size = 15)+
  theme(legend.title=element_blank())

ggarrange(epi,hypo,common.legend = TRUE)

ggsave("./Figs/p_SensitivityTest.jpg",width=12,height=7,units="in",dpi=320)
