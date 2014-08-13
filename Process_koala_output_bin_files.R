### Script for processing supercomputer output
# Edited by NB 14/8/14 to alter food intake calculations & calculate indexes of suitabiliy

args <- (commandArgs(TRUE))
simnum<-as.numeric(args[1])
bioreg<-as.numeric(args[2])
bioregion<-bioreg
sim<-as.numeric(args[3])

barcoo<-paste('/scratch/VR0212/bio',bioregion,'/',sep="")
load(paste(barcoo,'longlat.bin',sep=''))

numsites<-ceiling(nrow(data)/2/1000)
jstart<-numsites*(simnum-1)+1
jfinish<-numsites*(simnum-1)+numsites

for(jobnum in jstart:jfinish){
  
  
  #bioreg<-0
  #sim<-1
  i<-jobnum
  setwd(paste('/scratch/VR0212/koalas/bioreg',bioreg,'/sim',sim,sep=""))
  library(data.table)
  file<-paste(i,'sim_',sim,'.bin',sep="")
  load(file)  # Loads file with summed daily values
  
  
  ## fix up column names
  #colnames(Aust_sites)[20:23]<-c("WTR_G_H","RAINWET","RAINFALL","MAX_PCTWET")
  # use data table to speed things up
  library(data.table)
  AUST2<-data.table(Aust_sites) # make a new data table
  Aust_sites_week_by_year<-as.data.frame(AUST2[,lapply(.SD,mean),by=list(Sites,Morph,Behav,years,weeks)]) # get weekly averages for each year
  Aust_sites_week_by_year2<-as.data.frame(AUST2[,lapply(.SD,sum),by=list(Sites,Morph,Behav,years,weeks)]) # get weekly averages for each year
  A<-cbind(Aust_sites_week_by_year[,c(1:21,23)],Aust_sites_week_by_year2[,22],Aust_sites_week_by_year[,c(24:26)])
  colnames(A)[23] <- c('RAINFALL')
  AUST3<-data.table(A)
  
  # check for any weeks where either heat stressed (6 all week) or cold stressed (1.3 all week)
  heatstress<-with(A,ifelse(GMULT>=6.0,1,0))
  coldstress<-with(A,ifelse(GMULT<=1.3,1,0))
  A<-cbind(A,heatstress,coldstress)   
  
  # First calculate food intake to meet energy and water requirements
  Food_EN<-(A$Met_kJ_H/0.45/20.1) # Dry matter intake g/day to meet energy requirements
  Food_EN_NOMILK<-((A$Met_kJ_H-A$LACENRGY_kJ_H)/0.45/20.1) # Dry matter intake g/day to meet energy requirements without energy lost in milk
  Met_WT<-((((A$Met_kJ_H/0.45)*1000)/21.2)*0.00067) # Metabolic water production
  Req_WT<-(A$WTR_G_H-Met_WT)/0.66 # Amount of water required from foliage, 66% assimilation
  # Foliage water content 
  # low = 0.466 (mean moisture content selected species in winter, Ellis et al 1995)
  # med = 0.56 (mean moisture content of mature foliage of selected species - Nagy & Martin 1985)
  # high = 0.66 (mean moisture conted of young foliage of selected species - Nagy & martin 1985 
  # & similar to range (0.642-0.674) for species with new growth in summer Ellis et al 1995)   
  Food_WTL<-(Req_WT/0.466)*(1-0.466) # Dry matter intake to meet water req if % dry matter is 0.53
  Food_WTM<-(Req_WT/0.56)*(1-0.56) # Dry matter intake to meet water req if % dry matter is 0.41
  Food_WTH<-(Req_WT/0.66)*(1-0.66) # Dry matter intake to meet water req if % dry matter is 0.35
  ## Also calculate protein balance
  ## Assimilation is 0.032g per g dry matter ingested (Nagy & Martin 1985)
  ## From Kearney et al 2011 (greater glider) nitrogen= 1.4% dry matter *6.25 to convert to protein, digestive efficiency 48%
  ## which gives 0.042g per g dry matter ingested, similar to value of Cork et al 1983 that is 
  ## 1.1-1.5% dry matter & 45% digestive efficiency (gives 0.03-0.042g per g dry matter ingested) 
  Food_N<-(A$LACPRO/0.032) # Daily protein required for milk /digested per g dry matter (Nagy & Martin)
  Food_N_min<-(A$LACPRO/0.042) # Protein required for milk/digested per g dry matter (Cork) 
  A<-cbind(A,Food_EN,Food_EN_NOMILK,Food_WTL,Food_WTM,Food_WTH,Food_N,Food_N_min)
  
  # # Now check if required food intake exceeds average or maximum
  # # Average = 41.4 g/kg^0.75/day (Cork, Hume & Dawson)
  ## NB. edit 11/8/14 Average now calculated from Ellis et al 1995 so both mass^1
  # Average = 25.45g/kg/day 
  # # Max = 46.76 g/kg/day (Krockenberger 2003)
  # 
  LimE<-with(A,ifelse(Food_EN/MASS_KG>25.45,1,0))
  LimWL<-with(A,ifelse((Food_WTL/MASS_KG)>25.45,1,0))
  LimWM<-with(A,ifelse((Food_WTM/MASS_KG)>25.45,1,0))
  LimWH<-with(A,ifelse((Food_WTH/MASS_KG)>25.45,1,0))
  Lim2E<-with(A,ifelse(Food_EN/MASS_KG>46.76,1,0))
  Lim2WL<-with(A,ifelse((Food_WTL/MASS_KG)>46.76,1,0))
  Lim2WM<-with(A,ifelse((Food_WTM/MASS_KG)>46.76,1,0))
  Lim2WH<-with(A,ifelse((Food_WTH/MASS_KG)>46.76,1,0))
  LimE_nomilk<-with(A,ifelse(Food_EN_NOMILK/MASS_KG>25.45,1,0))
  Lim2E_nomilk<-with(A,ifelse(Food_EN_NOMILK/MASS_KG>46.76,1,0))
  LimN<-with(A,ifelse((Food_N/MASS_KG)>25.45,1,0))
  LimN_min<-with(A,ifelse((Food_N_min/MASS_KG)>25.45,1,0))
  Lim2N<-with(A,ifelse((Food_N/MASS_KG)>46.76,1,0))
  Lim2N_min<-with(A,ifelse((Food_N_min/MASS_KG)>46.76,1,0))
  
  ## if posture is 1.3 all week (curled in ball) can't feed so limit is 1
  LimE<-with(A,ifelse((coldstress)>0,1,LimE))
  Lim2E<-with(A,ifelse((coldstress)>0,1,Lim2E))
  LimE_nomilk<-with(A,ifelse((coldstress)>0,1,LimE_nomilk))
  Lim2E_nomilk<-with(A,ifelse((coldstress)>0,1,Lim2E_nomilk))
  LimN<-with(A,ifelse((coldstress)>0,1,LimN))
  LimN_min<-with(A,ifelse((coldstress)>0,1,LimN_min))
  Lim2N<-with(A,ifelse((coldstress)>0,1,Lim2N))
  Lim2N_min<-with(A,ifelse((coldstress)>0,1,Lim2N_min))
  
  
  # # create a new variable summing up rainfall in this week and week previous 
  rain<-as.data.frame(cbind(as.ts(A$RAINFALL),lag(A$RAINFALL,-1)))
  rain<-rain[1:1040,]
  colnames(rain)[1:2] <- c('RAIN','RAIN_PREV')
  rainsum<-rain$RAIN+rain$RAIN_PREV
  rainsum[1]<-rainsum[2]
  # rain[1040]<rain[1039]
  rain<-cbind(rain,rainsum)
  Water2<-with(rain,ifelse(rainsum>=1,1,0)) # if rainfall >=2mm water present in last two weeks
  Water<-with(rain,ifelse(rain$RAIN>=1,1,0)) # if rainfall >=1mm water present (last week)
  # 
  
  A<-cbind(A,LimE,LimWL,LimWM,LimWH,Lim2E,Lim2WL,Lim2WM,Lim2WH,LimE_nomilk,Lim2E_nomilk,LimN,LimN_min,Lim2N,Lim2N_min,Water,Water2)
  # # Now re-do water turnover variables taking into account if water is present & posture
  # # if water around limit is now set to 0, otherwise remains as calculated
  LimWLact<-with(A,ifelse((Water)>0,0,LimWL))
  LimWMact<-with(A,ifelse((Water)>0,0,LimWM))
  LimWHact<-with(A,ifelse((Water)>0,0,LimWH))
  Lim2WLact<-with(A,ifelse((Water)>0,0,Lim2WL))
  Lim2WMact<-with(A,ifelse((Water)>0,0,Lim2WM))
  Lim2WHact<-with(A,ifelse((Water)>0,0,Lim2WH))
  Lim2WHact_rain<-with(A,ifelse((Water2)>0,0,Lim2WH))
  A<-cbind(A,LimWLact,LimWMact,LimWHact,Lim2WLact,Lim2WMact,Lim2WHact,Lim2WHact_rain)
  # 
  # # Finally, re-calculate taking into account possible physiological changes that decrease water lost in faeces and urine
  Req_WT_min<-(A$WTR_G_H-Met_WT)/0.76 # Amount of water required from foliage
  # 
  Food_WTL_min<-(Req_WT_min/0.47)*(1-0.47) # Dry matter intake to meet water req if % dry matter is 0.53
  Food_WTM_min<-(Req_WT_min/0.59)*(1-0.59) # Dry matter intake to meet water req if % dry matter is 0.41
  Food_WTH_min<-(Req_WT_min/0.65)*(1-0.65) # Dry matter intake to meet water req if % dry matter is 0.35
  # 
  A<-cbind(A,Food_WTL_min,Food_WTM_min,Food_WTH_min)
  # # Now check if required food intake exceeds average or maximum
  # # Average = 42.2 g/kg^0.75/day (this study)
  # # Max = 46.76 g/kg/day (Krockenberger 2003) - make sure adjust for /kg instead of /kg^0.75
  # 
  LimWL_min<-with(A,ifelse((Food_WTL_min/MASS_KG)>25.45,1,0))
  LimWM_min<-with(A,ifelse((Food_WTM_min/MASS_KG)>25.45,1,0))
  LimWH_min<-with(A,ifelse((Food_WTH_min/MASS_KG)>25.45,1,0))
  Lim2WL_min<-with(A,ifelse((Food_WTL_min/MASS_KG)>46.76,1,0))
  Lim2WM_min<-with(A,ifelse((Food_WTM_min/MASS_KG)>46.76,1,0))
  Lim2WH_min<-with(A,ifelse((Food_WTH_min/MASS_KG)>46.76,1,0))
  # 
  
  A<-cbind(A,LimWL_min,LimWM_min,LimWH_min,Lim2WL_min,Lim2WM_min,Lim2WH_min)
  # # Now re-do water turnover variables taking into account if water is present
  # # if water around limit is now set to 0, otherwise remains as calculated
  LimWLact_min<-with(A,ifelse((Water)>0,0,LimWL_min))
  LimWMact_min<-with(A,ifelse((Water)>0,0,LimWM_min))
  LimWHact_min<-with(A,ifelse((Water)>0,0,LimWH_min))
  Lim2WLact_min<-with(A,ifelse((Water)>0,0,Lim2WL_min))
  Lim2WMact_min<-with(A,ifelse((Water)>0,0,Lim2WM_min))
  Lim2WHact_min<-with(A,ifelse((Water)>0,0,Lim2WH_min))
  Lim2WHact_rain_min<-with(A,ifelse((Water2)>0,0,Lim2WH_min))
  
  A<-cbind(A,LimWLact_min,LimWMact_min,LimWHact_min,Lim2WLact_min,Lim2WMact_min,Lim2WHact_min,Lim2WHact_rain_min)
  
  # calculate energy and water deficits in g dry food/mass (depends on limit)
  DefE<-25.45-(A$Food_EN/A$MASS_KG)
  DefWL<-25.45-(A$Food_WTL/A$MASS_KG)
  DefWM<-25.45-(A$Food_WTM/A$MASS_KG)
  DefWH<-25.45- (A$Food_WTH/A$MASS_KG)
  Def2E<-46.76-(A$Food_EN/A$MASS_KG)
  Def2WL<-46.76-(A$Food_WTL/A$MASS_KG)
  Def2WM<-46.76-(A$Food_WTM/A$MASS_KG)
  Def2WH<-46.76-(A$Food_WTH/A$MASS_KG)
  DefE_nomilk<-25.45-(A$Food_EN_NOMILK/A$MASS_KG)
  Def2E_nomilk<-46.76-(A$Food_EN_NOMILK/A$MASS_KG)
  
  A<-cbind(A,DefE,DefWL,DefWM,DefWH,Def2E,Def2WL,Def2WM,Def2WH,DefE_nomilk,Def2E_nomilk)
  
  # # Now re-do water turnover variables taking into account if water is present & posture
  # # if water around limit is now set to 0, otherwise remains as calculated
  DefWLact<-with(A,ifelse((Water)>0 & DefWL<0,0,DefWL))
  DefWMact<-with(A,ifelse((Water)>0 & DefWL<0,0,DefWM))
  DefWHact<-with(A,ifelse((Water)>0 & DefWL<0,0,DefWH))
  Def2WLact<-with(A,ifelse((Water)>0 & DefWL<0,0,Def2WL))
  Def2WMact<-with(A,ifelse((Water)>0 & DefWL<0,0,Def2WM))
  Def2WHact<-with(A,ifelse((Water)>0 & DefWL<0,0,Def2WH))
  Def2WHact_rain<-with(A,ifelse((Water2)>0 & DefWL<0,0,Def2WH))
  
  # Also calculate deficit using minimum required water values
  DefWL_min<-25.45-(A$Food_WTL_min/A$MASS_KG)
  DefWM_min<-25.45-(A$Food_WTM_min/A$MASS_KG)
  DefWH_min<-25.45- (A$Food_WTH_min/A$MASS_KG)
  Def2E_min<-46.76-(A$Food_EN_min/A$MASS_KG)
  Def2WL_min<-46.76-(A$Food_WTL_min/A$MASS_KG)
  Def2WM_min<-46.76-(A$Food_WTM_min/A$MASS_KG)
  Def2WH_min<-46.76-(A$Food_WTH_min/A$MASS_KG)
  
  A<-cbind(A,DefWL_min,DefWM_min,DefWH_min,Def2WL_min,Def2WM_min,Def2WH_min)
  
  # # Now re-do water deficits taking into account if water is present
  # # if water around limit is now set to 0, otherwise remains as calculated
  DefWLact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,DefWL_min))
  DefWMact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,DefWM_min))
  DefWHact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,DefWH_min))
  Def2WLact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,Def2WL_min))
  Def2WMact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,Def2WM_min))
  Def2WHact_min<-with(A,ifelse((Water)>0 & DefWL_min<0,0,Def2WH_min))
  Def2WHact_rain_min<-with(A,ifelse((Water2)>0 & DefWL_min<0,0,Def2WH_min))
  DefN<-25.45-(A$Food_N/A$MASS_KG)
  DefN_min<-25.45-(A$Food_N_min/A$MASS_KG)
  Def2N<-46.76-(A$Food_N/A$MASS_KG)
  Def2N_min<-46.76-(A$Food_N_min/A$MASS_KG)
  
  A<-cbind(A,DefWLact,DefWMact,DefWHact,Def2WLact,Def2WMact,Def2WHact,Def2WHact_rain)
  
  A<-cbind(A,DefWLact_min,DefWMact_min,DefWHact_min,Def2WLact_min,Def2WMact_min,Def2WHact_min,Def2WHact_rain_min,DefN,DefN_min,Def2N,Def2N_min) # Add min rain adjusted values
  
  
  # Summarize by year
  LIMITS<-data.table(A)
  Limits<-as.data.frame(LIMITS[,lapply(.SD,sum),by=list(Sites,Morph,Behav,years)]) # sum up how many weeks of the year are limited under each scenario
  Limits_mean<-as.data.frame(LIMITS[,lapply(.SD,mean),by=list(Sites,Morph,Behav,years)]) # get mean values for output (e.g. lat, long, rainfall)
  Limits_max<-as.data.frame(LIMITS[,lapply(.SD,min),by=list(Sites,Morph,Behav,years)]) 
  Koala_limits_yearly<-cbind(Limits_mean[,c(1:4,6,9:20,29:35,59:61)],Limits[,c(27:28,50:51,36:49,52:58,62:74)],Limits_mean[,c(75:108)],Limits_max[c(75:108)]) # paste output together & write to output file
  
  colnames(Koala_limits_yearly)[100:133]<-c("DefE_max","DefWL_max","DefWM_max","DefWH_max","Def2E_max","Def2WL_max",            
                                            "Def2WM_max","Def2WH_max","DefE_nomilk_max", "Def2E_nomilk_max","DefWLmin_max", "DefWMmin_max", "DefWHmin_max","Def2WLmin_max",       
                                            "Def2WMmin_max","Def2WHmin_max", "DefWLact_min_max", "DefWMact_max","DefWHact_max",
                                            "Def2WLact_max", "Def2WMact_max","Def2WHact_max","Def2WHact_rain_max", "DefWLact_min_max",      
                                            "DefWMact_min_max","DefWHact_min_max","Def2WLact_min_max", "Def2WMact_min_max",     
                                            "Def2WHact_min_max","Def2WHact_rain_min_max","DefN_max","DefN_min_max","Def2N_max","Def2N_min_max")
  
  
  ## Add columns calculating Water Index, Therm Index & limiting factors 
  #Firstcalculate habitat suitability index based on how easily can meet water requirements
  # 0 = under no circumstances
  # 1 = only with high food intake and high foliage water content
  # 2 = only with high food intake and medium foliage water content
  # 3 = only with medium food intake and medium foliage water content
  # 4 = with medium food intake and low foliage water content 
  
  WaterIndex<-with(Koala_limits_yearly,ifelse(Lim2WHact_min>1,0,1))
  WaterIndex<-with(Koala_limits_yearly,ifelse(Lim2WHact_min<=1,1,WaterIndex))
  WaterIndex<-with(Koala_limits_yearly,ifelse(Lim2WMact_min<=1,2,WaterIndex))
  WaterIndex<-with(Koala_limits_yearly,ifelse(LimWMact_min<=1,3,WaterIndex))
  WaterIndex<-with(Koala_limits_yearly,ifelse(LimWLact_min<=1,4,WaterIndex))
  
  Koala_limits_yearly<-cbind(Koala_limits_yearly,WaterIndex)
  
  ## Now check which factor is limiting across the year
  Limiting_low<-with(Koala_limits_yearly,ifelse(Food_EN>Food_WTL_min,0,1))
  Limiting_low<-with(Koala_limits_yearly,ifelse(Food_N>Food_WTL_min & Food_N>Food_EN,2,Limiting_low))
  Limiting_med<-with(Koala_limits_yearly,ifelse(Food_EN>Food_WTM_min,0,1))
  Limiting_med<-with(Koala_limits_yearly,ifelse(Food_N>Food_WTM_min & Food_N>Food_EN,2,Limiting_med))
  Limiting_high<-with(Koala_limits_yearly,ifelse(Food_EN>Food_WTH_min,0,1))
  Limiting_high<-with(Koala_limits_yearly,ifelse(Food_N>Food_WTH_min & Food_N>Food_EN,2,Limiting_high))
  
  Koala_limits_yearly<-cbind(Koala_limits_yearly,Limiting_low,Limiting_med,Limiting_high)
  
  ## Calculate habitat suitability index based on how easily can meet water AND energy requirements
  # including all lactation costs (no compensation)
  ThermIndex<-with(Koala_limits_yearly,ifelse(Lim2WHact_min>1,0,1))
  ThermIndex<-with(Koala_limits_yearly,ifelse(Lim2WHact_min<=1,1,ThermIndex))
  ThermIndex<-with(Koala_limits_yearly,ifelse(Lim2WMact_min<=1,2,ThermIndex))
  ThermIndex<-with(Koala_limits_yearly,ifelse(LimWMact_min<=1,3,ThermIndex))
  ThermIndex<-with(Koala_limits_yearly,ifelse(LimWLact_min<=1,4,ThermIndex))
  ThermIndex<-with(Koala_limits_yearly,ifelse(Lim2E>1,0,ThermIndex))
  ThermIndex<-with(Koala_limits_yearly,ifelse(LimE>1 & ThermIndex>2,2,ThermIndex))
  
  Koala_limits_yearly<-cbind(Koala_limits_yearly,ThermIndex)
  
  # excluding costs of energy exported in milk (compensation)
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(Lim2WHact_min>1,0,1))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(Lim2WHact_min<=1,1,ThermIndex_nomilk))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(Lim2WMact_min<=1,2,ThermIndex_nomilk))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(LimWMact_min<=1,3,ThermIndex_nomilk))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(LimWLact_min<=1,4,ThermIndex_nomilk))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(Lim2E_nomilk>1,0,ThermIndex_nomilk))
  ThermIndex_nomilk<-with(Koala_limits_yearly,ifelse(LimE_nomilk>1 & ThermIndex_nomilk>2,2,ThermIndex_nomilk))
  
  Koala_limits_yearly<-cbind(Koala_limits_yearly,ThermIndex_nomilk)
  
  ## Sum of weeks in negative energy or water balance/food intake above average
  UpIntake_wks<-sum(Koala_limits_yearly$LimE,Koala_limits_yearly$LimWMact_min)
  NegBalance_wks<-sum(Koala_limits_yearly$Lim2E,Koala_limits_yearly$Lim2WMact_min)
  UpIntake_wks_nomilk<-sum(Koala_limits_yearly$LimE_nomilk,Koala_limits_yearly$LimWMact_min)
  NegBalance_wks_nomilk<-sum(Koala_limits_yearly$Lim2E_nomilk,Koala_limits_yearly$Lim2WMact_min)
  
  Koala_limits_yearly<-cbind(Koala_limits_yearly,UpIntake_wks,UpIntake_wks_nomilk,NegBalance_wks,NegBalance_wks_nomilk)
    
  
  ## Output of max energy and water requirements (daily & weekly)
  Day_max<-as.data.frame(AUST2[,lapply(.SD,max),by=list(Sites,Morph,Behav)])
  Week_max<-as.data.frame(AUST3[,lapply(.SD,max),by=list(Sites,Morph,Behav)])
  Maximums<-cbind(Day_max[c(1:4,8:9,11:23)],Week_max[12:23])
  colnames(Maximums)[20:31] <- c("MET_W_WK","PCTWET_WK", "PCTSHD_WK","GMULT_WK","MET_BMR_WK", "Met_kJ_H_WK", "EVP_G_H_WK",
                                 "WCUT_G_H_WK","WTR_G_H_WK","RAINWET_WK","MAX_PCTWET_WK","RAINFALL_WK")
  
  Aust_sites_max_wet_hour<-as.data.frame(AUST2[,.SD[which.max(MAX_PCTWET)],by=list(Sites,Morph,Behav,Repro)])
  Aust_sites_max_wet_day<-as.data.frame(AUST2[,.SD[which.max(PCTWET)],by=list(Sites,Morph,Behav,Repro)])
  Aust_sites_max_water_day<-as.data.frame(AUST2[,.SD[which.max(WTR_G_H)],by=list(Sites,Morph,Behav,Repro)])
  Aust_sites_max_energy_day<-as.data.frame(AUST2[,.SD[which.max(Met_kJ_H)],by=list(Sites,Morph,Behav,Repro)])
  # # 
  Aust_sites_max_wet_week<-as.data.frame(AUST3[,.SD[which.max(PCTWET)],by=list(Sites,Morph,Behav,Repro)])
  Aust_sites_max_water_week<-as.data.frame(AUST3[,.SD[which.max(WTR_G_H)],by=list(Sites,Morph,Behav,Repro)])
  Aust_sites_max_energy_week<-as.data.frame(AUST3[,.SD[which.max(Met_kJ_H)],by=list(Sites,Morph,Behav,Repro)])
  # # 
  Aust_sites_max<-cbind(Aust_sites_max_wet_hour[c(1:5,7:11,23)],Aust_sites_max_wet_day[c(5,7,13)],Aust_sites_max_water_day[c(5,7,20)],Aust_sites_max_energy_day[c(5,7,16,17)],
                        Aust_sites_max_wet_week[c(5,6,13,14)],Aust_sites_max_water_week[c(5,6,20)],Aust_sites_max_energy_week[c(5,6,16:17)])
  
  names(Aust_sites_max)<-c("Sites","Morph","Behav","Repro","year_maxPCTWT","day_maxPCTWT", "latitude","longitude","weeks_maxPCTWT", "MASS_KG","MaxPCTWT",
                           "year_PCTWET_day","day_PCTWET_day","PCTWET_day","year_maxWTR_day","day_maxWTR_day","maxWTR_day","year_maxMET_day","day_maxMET_day",
                           "maxMET_BMR_day","maxMET_day","year_maxPCTWET","week_maxPCTWET","maxPCTWET","maxPCTSHD",
                           "year_maxWTR","week_maxWTR","maxWTR","year_maxMET", "week_maxMET","maxMET_BMR","maxMET")
  setwd(paste('/scratch/VR0212/koalas/bioreg',bioreg,'/',sep=""))
  write.table(Koala_limits_yearly, file = paste("Koala_limits_yearly_bio",bioreg,"_sim",sim,".csv",sep=""), sep = ",", col.names = F, qmethod = "double", append = T)
  write.table(Aust_sites_max, file = paste("Aust_sites_max_bio",bioreg,"_sim",sim,".csv",sep=""), sep = ",", col.names = F, qmethod = "double", append = T)
  
}
