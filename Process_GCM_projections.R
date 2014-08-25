### Script for processing individual GCM climate projection files ###########
library(adehabitat)
library(raster)
# Set Global circulation model 
GCM<-"HadGEM2-CC"
GCM_name<-"HadGEM2-CC"
setwd(paste("F:/",GCM,sep=""))

## TMIN##

# read in 2050 predictions & stack
minTst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 1:12){ 
  m<-months[i]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvTMin ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
#  plot(r)
  minTst_2050<-stack(minTst_2050,r)
}
# read in 2070 predictions and stack
minTst_2070<-stack()
for(i in 1:12){ 
  m<-months[i]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvTMin ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
#  plot(r)
  minTst_2070<-stack(minTst_2070,r)
}
# read in "baselne" data (2000)
minTst_base<-stack()
for(i in 1:12){ 
m<-months[i]  
 file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvTMin ",GCM_name," RCP8.5 ",m,".asc",sep="")
 b <- import.asc(file) 
 r <- raster(b)
# plot(r)
 minTst_base<-stack(minTst_base,r)
}
# subtract baseline to get absolute monthly changes
minT_diff_2050<-minTst_2050-minTst_base
plot(minT_diff_2050)

minT_diff_2070<-minTst_2070-minTst_base
plot(minT_diff_2070)

for(i in 1:12){
writeRaster(minT_diff_2050[[i]], paste("F:/",GCM,"/2050/2050_Tmin",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(minT_diff_2070[[i]], paste("F:/",GCM,"/2070/2070_Tmin",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))])

## TMAX##

# read in 2050 predictions & stack
maxTst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 13:24){ 
  m<-months[i-12]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvTMax ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  maxTst_2050<-stack(maxTst_2050,r)
}
# read in 2070 predictions and stack
maxTst_2070<-stack()
for(i in 13:24){ 
  m<-months[i-12]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvTMax ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  maxTst_2070<-stack(maxTst_2070,r)
}
# read in "baseline" data (2000)
maxTst_base<-stack()
for(i in 13:24){ 
  m<-months[i-12]  
  file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvTMax ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  # plot(r)
  maxTst_base<-stack(maxTst_base,r)
}
  # subtract baseline to get absolute monthly changes
  maxT_diff_2050<-maxTst_2050-maxTst_base
  plot(maxT_diff_2050)


  maxT_diff_2070<-maxTst_2070-maxTst_base
  plot(maxT_diff_2070)
  
for(i in 1:12){
  writeRaster(maxT_diff_2050[[i]], paste("F:/",GCM,"/2050/2050_Tmax",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(maxT_diff_2070[[i]], paste("F:/",GCM,"/2070/2070_Tmax",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))])

## RH ###########
RHst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvRelHum ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  RHst_2050<-stack(RHst_2050,r)
}
# read in 2070 predictions and stack
RHst_2070<-stack()
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvRelHum ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  RHst_2070<-stack(RHst_2070,r)
}
# read in "baseline" data (2000)
RHst_base<-stack()
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvRelHum ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  # plot(r)
  RHst_base<-stack(RHst_base,r)
}
# subtract baseline to get absolute monthly changes
RH_diff_2050<-RHst_2050-RHst_base
plot(RH_diff_2050)


RH_diff_2070<-RHst_2070-RHst_base
plot(RH_diff_2070)

for(i in 1:12){
  writeRaster(RH_diff_2050[[i]], paste("F:/",GCM,"/2050/2050_RelHum",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(RH_diff_2070[[i]], paste("F:/",GCM,"/2070/2070_RelHum",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))])


### Wind ###
VELst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 49:60){ 
  m<-months[i-48]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvWind ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  VELst_2050<-stack(VELst_2050,r)
}
# read in 2070 predictions and stack
VELst_2070<-stack()
for(i in 49:60){ 
  m<-months[i-48]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvWind ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  VELst_2070<-stack(VELst_2070,r)
}
# read in "baseline" data (2000)
VELst_base<-stack()
for(i in 49:60){ 
  m<-months[i-48]  
  file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvWind ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  # plot(r)
  VELst_base<-stack(VELst_base,r)
}

# divide by baseline to get proportional monthly changes? 
VEL_prop_2050<-VELst_2050/VELst_base
plot(VEL_prop_2050)


VEL_prop_2070<-VELst_2070/VELst_base
plot(VEL_prop_2070)

for(i in 1:12){
  writeRaster(VEL_prop_2050[[i]], paste("F:/",GCM,"/2050/2050_Wind",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(VEL_prop_2070[[i]], paste("F:/",GCM,"/2070/2070_Wind",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))])

### RAIN ###
RAINst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 25:36){ 
  m<-months[i-24]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvPrecip ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  RAINst_2050<-stack(RAINst_2050,r)
}
# read in 2070 predictions and stack
RAINst_2070<-stack()
for(i in 25:36){ 
  m<-months[i-24]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvPrecip ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  RAINst_2070<-stack(RAINst_2070,r)
}
# read in "baseline" data (2000)
RAINst_base<-stack()
for(i in 25:36){ 
  m<-months[i-24]  
  file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvPrecip ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  # plot(r)
  RAINst_base<-stack(RAINst_base,r)
}


# divide by baseline to get proportional monthly changes? Q. What happens if one is 0?
RAIN_prop_2050<-RAINst_2050/RAINst_base
plot(RAIN_prop_2050)


RAIN_prop_2070<-RAINst_2070/RAINst_base
plot(RAIN_prop_2070)

for(i in 1:12){
  writeRaster(RAIN_prop_2050[[i]], paste("F:/",GCM,"/2050/2050_Precip",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(RAIN_prop_2070[[i]], paste("F:/",GCM,"/2070/2070_Precip",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))])


### Solar ###
SOLst_2050<-stack()
months<-c("January","February", "March","April","May","June","July","August","September","October","November","December")
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2050/",i," Linked Australia 2050 cvSolar ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  SOLst_2050<-stack(SOLst_2050,r)
}
# read in 2070 predictions and stack
SOLst_2070<-stack()
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2070/",i," Linked Australia 2070 cvSolar ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  #  plot(r)
  SOLst_2070<-stack(SOLst_2070,r)
}
# read in "baseline" data (2000)
SOLst_base<-stack()
for(i in 37:48){ 
  m<-months[i-36]  
  file <-  paste("F:/",GCM,"/",GCM, " 2000/",i," Linked Australia 2000 cvSolar ",GCM_name," RCP8.5 ",m,".asc",sep="")
  b <- import.asc(file) 
  r <- raster(b)
  # plot(r)
  SOLst_base<-stack(SOLst_base,r)
}
# subtract baseline to get absolute monthly changes
SOL_diff_2050<-SOLst_2050-SOLst_base
plot(SOL_diff_2050)


SOL_diff_2070<-SOLst_2070-SOLst_base
plot(SOL_diff_2070)

for(i in 1:12){
  writeRaster(SOL_diff_2050[[i]], paste("F:/",GCM,"/2050/2050_Solar",i,".tif",sep=""))
}
for(i in 1:12){
  writeRaster(SOL_diff_2070[[i]], paste("F:/",GCM,"/2070/2070_Solar",i,".tif",sep=""))
}

rm(list= ls()[!(ls() %in% c('GCM','GCM_name'))]

