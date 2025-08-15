# FCR_DOC

Scripts to support Hounshell et al. in Review.


Run scripts in the following order: 

1_LakeAnalyzer_thermo.R - format sensor data (YSI, CTD) and calculate variable thermocline depth using R LakeAnalyzer;

1a_Q_Error.R - To estimate error from inflow measurements;

2_Eco_DOC_rlnorm_FC.R - Develop mass balance model to estimate internal DOC contributions to Falling Creek Reservoir;

3_Env_ARIMA_FC.R - Synthesize and visualize mass balance model results and conduct ARIMA time series modeling to identify key environmental parameters for [DOC] and internal DOC;

4_OpticalData.R - Script to visualize DOM quality data (fluorescent dissolved organic matter) for sumemr 2019 to support modeling and time series results.


Scripts used for the original manuscript submission can be found in /Scripts/Original. Key, generated data-sets, including from the developed DOC ecosystem model can be found in the /Data/ folder.
