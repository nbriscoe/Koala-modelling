#############Comparison of NicheMapper & Maxent models for koalas #########################
# load useful libraries
library(raster)
library(rgdal)
library(adehabitat)
library(pROC)
library(dismo)
library(mgcv)

setwd("C:/Users/nbriscoe.UNIMELB/Documents/Koala modelling/github/Comparison NicheMapper & Maxent")
wd<-getwd()

# get bioregion outlines & mask
shfilename<-("bg_bioregs_single")
bioregions<-readOGR(dsn="C:/Users/nbriscoe.UNIMELB/Documents/Koala modelling/github/Comparison NicheMapper & Maxent",layer=shfilename)
myshp <- spTransform(bioregions, CRS=CRS("+proj=longlat +ellps=WGS84 +no_defs +towgs84=0,0,0,0,0,0,0"))
r <- raster(extent(myshp))
res(r)=0.05
r[]<-0
r <- rasterize(myshp, field=1, r)


shfilename<-("bg_bioregs_single")
Aust_outline<-readOGR(dsn="C:/Users/nbriscoe.UNIMELB/Documents/Koala modelling/Climate refugia/koala_maxent",layer=shfilename)
myshp_Aust <- spTransform(Aust_outline, CRS=CRS("+proj=longlat +ellps=WGS84 +no_defs +towgs84=0,0,0,0,0,0,0"))
Aus <- raster(extent(myshp_Aust))
res(Aus)=0.05
Aus[]<-0
Aus <- rasterize(myshp_Aust, field=1, Aus)

# Read in koala points for analysis
kpts_sub<-read.csv("kpts_sub.csv") # sub-set of points used to build Maxent model (temporally & spatially filtered - accuracy <=1km & after 1960)
kpts_all<-read.csv("Koala_points_1km_accuracy.csv") # all koala observation records accurate to <=1km
Pres_abs_fam<-read.csv("Ale_presabs_fam_10km_crop.csv") # Data file with presence points of other arboreal mammals 
# (includes species in the families Phalangeridae, Petauridae and Pseudocheiridae as well as the tree kangaroos - cropped to exclude records within 10km of a koala obs)

## First read in and plot predictions
# make a stack of NicheMapper yearly predictions
NM.st<-stack()
names <- list.files(path=file.path(wd, '/UpIntake_Index_nomilk_current'), full.names=TRUE)
#now make the stacK:
NM.st <- stack(names) 
#get mean across years
meanTI<-mean(NM.st)

# make a stack of NicheMapper yearly predictions - 2050
NM_2050.st<-stack()
names <- list.files(path=file.path(wd, '/UpIntake_Index_nomilk_2050'), full.names=TRUE)
#now make the stacK:
NM_2050.st <- stack(names[1:20]) 
#get mean across years
names(NM_2050.st)
meanTI_2050<-mean(NM_2050.st)


# make a stack of NicheMapper yearly predictions - 2050
NM_2070.st<-stack()
names <- list.files(path=file.path(wd, '/UpIntake_Index_nomilk_2070'), full.names=TRUE)
#now make the stacK:
NM_2070.st <- stack(names[1:20]) 
#get mean across years
names(NM_2070.st)
meanTI_2070<-mean(NM_2070.st)

# current
mod7_pred_curr<-raster(paste(wd,"/Maxent_predictions/mod7_current.tif",sep=""))
mod7E_pred_curr<-raster(paste(wd,"/Maxent_predictions/mod7E_current.tif",sep=""))
mod9E_pred_curr<-raster(paste(wd,"/Maxent_predictions/mod9E_current.tif",sep=""))
# 2050 - GCM = ACCESS 1.3
mod7_pred_2050<-raster(paste(wd,"/Maxent_predictions/mod7_ACCESS1.3_2050.tif",sep=""))
mod7E_pred_2050<-raster(paste(wd,"/Maxent_predictions/mod7E_ACCESS1.3_2050.tif",sep=""))
mod9E_pred_2050<-raster(paste(wd,"/Maxent_predictions/mod9E_ACCESS1.3_2050.tif",sep=""))
# 2070 - GCM = ACCESS 1.3
mod7_pred_2070<-raster(paste(wd,"/Maxent_predictions/mod7_ACCESS1.3_2070.tif",sep=""))
mod7E_pred_2070<-raster(paste(wd,"/Maxent_predictions/mod7E_ACCESS1.3_2070.tif",sep=""))
mod9E_pred_2070<-raster(paste(wd,"/Maxent_predictions/mod9E_ACCESS1.3_2070.tif",sep=""))

######################################################################################
### Also try calculating a sliding window mean - if unsuitable in more than 2 years in a row becomes unsuitable##
## better to get sum of poor years in each 5, then multiply these together. If there are 3 poor years in a row then cannot be there. 
gen<-5
seq_poor_years.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM.st[[st:b]]
  gen_suit<-sum(layers)/gen
  
  # plot(gen_suit)
  seq_poor_years.st<-stack(seq_poor_years.st,gen_suit)
}
plot(seq_poor_years.st)

min_poor_yrs<-min(seq_poor_years.st)
plot(min_poor_yrs)
points(kpts_sub,cex=0.1)

## Use average habitat suitability but set places with >2 years in a row unsuitable to 0
min_poor_yrs[min_poor_yrs>0]<-1
plot(min_poor_yrs)

meanTI_poor<-min_poor_yrs*meanTI

## Also calculate geometric mean
for (n in 1:(20-gen+1)){  
  if (n==1){
    gmean<-seq_poor_years.st[[n]]
  }else{
    gmean<-gmean*seq_poor_years.st[[n]]
  }
  }
  gmean<-gmean^(1/(20-gen+1))
plot(gmean)
plot(meanTI_poor)
points(kpts_sub,cex=0.1)
### Also calculate for 2050 ################
gen<-5
seq_poor_years.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM_2050.st[[st:b]]
  gen_suit<-sum(layers)/gen
  # plot(gen_suit)
  seq_poor_years.st<-stack(seq_poor_years.st,gen_suit)
}
plot(seq_poor_years.st)

min_poor_yrs_2050<-min(seq_poor_years.st)
plot(min_poor_yrs_2050)
points(kpts_sub,cex=0.1)

## Also calculate geometric mean
for (n in 1:(20-gen+1)){  
  if (n==1){
    gmean_2050<-seq_poor_years.st[[n]]
  }else{
    gmean_2050<-gmean_2050*seq_poor_years.st[[n]]
  }
}
gmean_2050<-gmean_2050^(1/(20-gen+1))
plot(gmean_2050)
plot(meanTI_2050_poor)
points(kpts_sub,cex=0.1)

## Use average habitat suitability but set places with >2 years in a row unsuitable to 0
min_poor_yrs_2050[min_poor_yrs_2050>0]<-1
plot(min_poor_yrs_2050)

meanTI_2050_poor<-min_poor_yrs_2050*meanTI_2050
plot(meanTI_2050_poor)
points(kpts_sub,cex=0.1)

#########Now do for 2070######################
gen<-5
seq_poor_years.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM_2070.st[[st:b]]
  gen_suit<-sum(layers)/gen
#  gen_suit<-(layers[[1]]*layers[[2]]*layers[[3]])^(1/3)
  # plot(gen_suit)
  seq_poor_years.st<-stack(seq_poor_years.st,gen_suit)
}
plot(seq_poor_years.st)

min_poor_yrs_2070<-min(seq_poor_years.st)
plot(min_poor_yrs_2070)
points(kpts_sub,cex=0.1)

## Use average habitat suitability but set places with >2 years in a row unsuitable to 0
min_poor_yrs_2070[min_poor_yrs_2070>0]<-1
plot(min_poor_yrs_2070)

meanTI_2070_poor<-min_poor_yrs_2070*meanTI_2070
plot(meanTI_2070_poor)
points(kpts_sub,cex=0.1)

## Also calculate geometric mean
for (n in 1:(20-gen+1)){  
  if (n==1){
    gmean_2070<-seq_poor_years.st[[n]]
  }else{
    gmean_2070<-gmean_2070*seq_poor_years.st[[n]]
  }
}
gmean_2070<-gmean_2070^(1/(20-gen+1))
plot(gmean_2070)
plot(meanTI_2070_poor)
points(kpts_sub,cex=0.1)

################################################################################################
# mean habitat suitability (average across years - but exclude where can't survive for 5 years)
par(mfrow=c(2,2))
plot(meanTI_poor)
plot(meanTI_2050_poor)
plot(meanTI_2070_poor)
plot(meanTI_poor-meanTI_2070_poor,col=bpy.colors(50))
#################################################################################################
par(mfrow=c(2,2))
plot(gmean)
plot(gmean_2050)
plot(gmean_2070)
plot(gmean-gmean_2070,col=bpy.colors(50))
################################### Plot predictions #######################################
meanTI_plot<-crop(meanTI,mod7_pred_curr)
meanTI_2050_plot<-crop(meanTI_2050,mod7_pred_curr)
meanTI_2070_plot<-crop(meanTI_2070,mod7_pred_curr)

gmean_plot<-crop(gmean,mod7_pred_curr)
gmean_2050_plot<-crop(gmean_2050,mod7_pred_curr)
gmean_2070_plot<-crop(gmean_2070,mod7_pred_curr)

meanTI_poor_plot<-crop(meanTI_poor,mod7_pred_curr)
meanTI_2050_poor_plot<-crop(meanTI_2050_poor,mod7_pred_curr)
meanTI_2070_poor_plot<-crop(meanTI_2070_poor,mod7_pred_curr)



windows(10,8)
par(mfrow=c(3,4),mar = c(0, 0, 0, 0))
breaks<-seq(0,1,length.out=100)
brks<-c(seq(0,0.7,length.out=100))
cols<-terrain.colors(100)
cols<-rev(cols)
plot(gmean_plot,breaks=breaks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
#plot(mod10_pred_curr,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7_pred_curr,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7E_pred_curr,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod9E_pred_curr,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)

plot(gmean_2050_plot,breaks=breaks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
#plot(mod10_pred_2050,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7_pred_2050,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7E_pred_2050,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod9E_pred_2050,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)

plot(gmean_2070_plot,breaks=breaks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
#plot(mod10_pred_2050,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7_pred_2070,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod7E_pred_2070,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)
plot(mod9E_pred_2070,xlim=c(130,160),breaks=brks,col=cols,legend=FALSE,box=FALSE,axes=FALSE)

dev.print(tiff, "ModelPredictions_hires_NM_milk_poor.tiff", compression = "lzw", res=600, height=11, width=17, units="cm")

## First calculate AUC for NicheMapper (with milk included)##########
## First use 10,000 randomly placed absence points & koala presences with average output across years
rpoints<- randomPoints(r,10000) 
names(kpts_sub)<-c("x","y")
test<-rbind(kpts_sub,rpoints)
S<-extract(meanTI,test)
S_poor<-extract(min_poor_yrs,test)
S_mean_poor<-extract(meanTI_poor,test)
S_gmean<-extract(gmean,test)
test_out<-cbind(cbind(c(rep(1,4387),rep(0,10000))),test,S,S_poor,S_mean_poor,S_gmean)
colnames(test_out)[1]<-"Type"

Unsuitable<-subset(test_out,S==0)
Unsuitable<-subset(Unsuitable,Type==1)

#NicheMapper - milk
ROC_S_milk_av<-roc(Type~S,data=test_out,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_S_milk_av

#NicheMapper - milk - min years
ROC_NM_poor_yrs<-roc(Type~S_poor,data=test_out,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_poor_yrs

#NicheMapper - milk - min years
ROC_NM_mean_poor<-roc(Type~S_mean_poor,data=test_out,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_mean_poor

#NicheMapper - milk - min years
ROC_NM_gmean<-roc(Type~S_gmean,data=test_out,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_gmean


## Also calculate AUC with year of or year prior to observation ###############
## First get a distribution of background samples that match koala observations in time#######
kpts_yrs<-subset(kpts_all,year>1989)
kpts_yrs<-subset(kpts_yrs,year<2010)
names(kpts_yrs)[1:2]<-c("x","y")
yrs<-kpts_yrs$year

bg_yr<-sample(yrs,100000,replace=TRUE)
rpoints_yr<- randomPoints(r,100000) 
bg_yrs<-as.data.frame(cbind(rpoints_yr,bg_yr))

yrs<-seq(1990,2009,1)
for(j in 2:20){
  print(j)
  yr<-yrs[[j]]
  WI_yr <- NM.st[[j]]
  WI_yr_p <- NM.st[[j-1]]
  points<-subset(bg_yrs,bg_yr==yr)
  locs<-cbind(points$x,points$y) 
  locs<- gridSample(locs, r, n=1) # sample one point in each 0.05 x 0.05 resolution cell (to match climate data) 
  years<-rep(yr,length(locs[,1]))
  WaterIndexP<-extract(WI_yr_p,locs)
  WaterIndex<-extract(WI_yr,locs)
  Predict_yr_prev<-cbind(locs,years,WaterIndex,WaterIndexP)
  #  paste("Predict_",yr,sep="")<-Predict_yr
  if(j==2){
    Predict_all_bg_yr<-Predict_yr_prev}else{
      Predict_all_bg_yr<-rbind(Predict_all_bg_yr,Predict_yr_prev)
    }
}

Predict_all_bg_yr<-as.data.frame(Predict_all_bg_yr)
Predict_all_bg_yr_crop<-subset(Predict_all_bg_yr,WaterIndex!="NA")
Predict_all_bg_yr_crop$WIav<-(Predict_all_bg_yr_crop$WaterIndexP+Predict_all_bg_yr_crop$WaterIndex)/2


## Now run for koala points
yrs<-seq(1990,2009,1)
for(j in 2:20){
  print(j)
  yr<-yrs[[j]]
  WI_yr <- NM.st[[j]]
  WI_yr_p <- NM.st[[j-1]]
  points<-subset(kpts_yrs,year==yr)
  locs<-cbind(points$x,points$y) 
  locs<- gridSample(locs, r, n=1) # sample one point in each 0.05 x 0.05 resolution cell (to match climate data) 
  years<-rep(yr,length(locs[,1]))
  WaterIndexP<-extract(WI_yr_p,locs)
  WaterIndex<-extract(WI_yr,locs)
  Predict_yr_prev<-cbind(locs,years,WaterIndex,WaterIndexP)
  #  paste("Predict_",yr,sep="")<-Predict_yr
  if(j==2){
    Predict_all_kpts_yr<-Predict_yr_prev}else{
      Predict_all_kpts_yr<-rbind(Predict_all_kpts_yr,Predict_yr_prev)
    }
}

Predict_all_kpts_yr<-as.data.frame(Predict_all_kpts_yr)
Predict_all_kpts_yr_crop<-subset(Predict_all_kpts_yr,WaterIndex!="NA")
Predict_all_kpts_yr_crop$WIav<-(Predict_all_kpts_yr_crop$WaterIndexP+Predict_all_kpts_yr_crop$WaterIndex)/2

#write.table(Predict_all_bg_yr_crop,"Yearly_NicheMapper_bg_WIscores_av_yrs_spatial.csv",sep=",",col.names=TRUE,row.names=FALSE)
NM_roc_year<-rbind(Predict_all_kpts_yr_crop,Predict_all_bg_yr_crop)
l1<-length(Predict_all_kpts_yr_crop)
l2<-length(Predict_all_bg_yr_crop)
NM_roc_year$Type<-c(rep(1,6471),rep(0,96428))
head(NM_roc_year)
## first calc AUC with average of prev & that year
ROC_NM_yr_av<-roc(Type~WIav,data=NM_roc_year,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_yr_av

## Also check just using prev year = slightly lower than average 
ROC_NM_yr_prev<-roc(Type~WaterIndexP,data=NM_roc_year,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_yr_prev

## Also check just using year = slightly lower than average 
ROC_NM_yr<-roc(Type~WaterIndex,data=NM_roc_year,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_yr

############Also use presence 'absence' points to test models ###############
abs_sub<- gridSample(Pres_abs_fam, r, n=1) # sample one point in each 0.05 x 0.05 resolution cell (to match climate data) 
# also get rid of points that fall outside the bioregion modelled (r)
extract(r,abs_sub)->x
abs_sub[which(!is.na(x)),]->abs_sub 
names(abs_sub)<-c("x","y")
# see where the points are
plot(myshp)
points(abs_sub,col="red",cex=0.1)
points(kpts_sub,col="black",cex=0.1)

test_abs<-rbind(abs_sub,kpts_sub)
NM<-extract(meanTI,test_abs)
NM_poor<-extract(meanTI_poor,test_abs)
NM_gmean<-extract(gmean,test_abs)
Max_7<-extract(mod7_pred_curr,test_abs)
Max_7E<-extract(mod7E_pred_curr,test_abs)
Max_9E<-extract(mod9E_pred_curr,test_abs)
test_out_abs<-cbind(cbind(c(rep(0,3285),rep(1,4387))),test_abs,NM,NM_poor,NM_gmean,Max_7,Max_7E,Max_9E)
colnames(test_out_abs)[1]<-"Type"

#NicheMapper - milk
ROC_NM_abs<-roc(Type~NM,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_abs

#NicheMapper - milk
ROC_NM_poor_abs<-roc(Type~NM_poor,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_poor_abs

#NicheMapper - gmean
ROC_NM_gmean_abs<-roc(Type~NM_gmean,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_NM_gmean_abs

#maxent averages (mod7)
ROC_Max7_abs<-roc(Type~Max_7,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_Max7_abs

#maxent extremes (mod7E)
ROC_Max7E_abs<-roc(Type~Max_7E,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_Max7E_abs

#maxent extremes (mod9E)
ROC_Max9E_abs<-roc(Type~Max_9E,data=test_out_abs,AUC=TRUE,ci=TRUE,plot=TRUE)
ROC_Max9E_abs

##################################################################################
# Get correlations between predictions
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (r/2))
}

NicheMapper<-extract(meanTI,rpoints,buffer=NULL)
NicheMapper_poor<-extract(meanTI_poor,rpoints,buffer=NULL)
Nichemapper_gmean<-extract(gmean,rpoints,buffer=NULL)
Maxent_av<-extract(mod7_pred_curr,rpoints,buffer=NULL)
Maxent_ext7<-extract(mod7E_pred_curr,rpoints,buffer=NULL)
Maxent_ext9<-extract(mod9E_pred_curr,rpoints,buffer=NULL)
Model_testR<-as.data.frame(cbind(NicheMapper,NicheMapper_poor,Nichemapper_gmean,Maxent_av,Maxent_ext7,Maxent_ext9))
Model_testR<-na.omit(Model_testR)
pairs(Model_testR, lower.panel=panel.smooth, upper.panel=panel.cor)  

NicheMapper<-extract(meanTI_2050,rpoints,buffer=NULL)
NicheMapper_poor<-extract(meanTI_2050_poor,rpoints,buffer=NULL)
NicheMapper_gmean<-extract(gmean_2050,rpoints,buffer=NULL)
Maxent_av<-extract(mod7_pred_2050,rpoints,buffer=NULL)
Maxent_ext7<-extract(mod7E_pred_2050,rpoints,buffer=NULL)
Maxent_ext9<-extract(mod9E_pred_2050,rpoints,buffer=NULL)
Model_testR_2050<-cbind(NicheMapper,NicheMapper_poor,NicheMapper_gmean,Maxent_av,Maxent_ext7,Maxent_ext9)
Model_testR_2050<-na.omit(Model_testR_2050)
pairs(Model_testR_2050, lower.panel=panel.smooth, upper.panel=panel.cor)  

NicheMapper<-extract(meanTI_2070,rpoints,buffer=NULL)
NicheMapper_poor<-extract(meanTI_2070_poor,rpoints,buffer=NULL)
NicheMapper_gmean<-extract(gmean_2070,rpoints,buffer=NULL)
Maxent_av<-extract(mod7_pred_2070,rpoints,buffer=NULL)
Maxent_ext7<-extract(mod7E_pred_2070,rpoints,buffer=NULL)
Maxent_ext9<-extract(mod9E_pred_2070,rpoints,buffer=NULL)
Model_testR_2070<-cbind(NicheMapper,NicheMapper_poor,NicheMapper_gmean,Maxent_av,Maxent_ext7,Maxent_ext9)
Model_testR_2070<-na.omit(Model_testR_2070)
pairs(Model_testR_2070, lower.panel=panel.smooth, upper.panel=panel.cor) 

## Pretty plots of relationship
mod<-gam(NicheMapper_gmean~s(Maxent_av),data=Model_testR)
plot(mod, se=TRUE)
dat<-as.data.frame(seq(0,0.7,length.out=1000))
colnames(dat)<-"Maxent_ext7"
p<-predict(mod, dat, type = "link", se.fit = TRUE)
upr <- mod$family$linkinv(upr)
lwr <- mod$family$linkinv(lwr)

##############Possibly the number of sequential weeks that koalas need to increase food intake to meet thermoreg is better indicator of suitability?#########
# make a stack of NicheMapper yearly predictions
NM_seq.st<-stack()
names <- list.files(path=file.path(wd, '/Seq_wks_negW_current'), full.names=TRUE)
#now make the stacK:
NM_seq.st <- stack(names) 
#get mean across years
meanTI_seq<-mean(NM_seq.st)

# make a stack of NicheMapper yearly predictions - 2050
NM_seq_2050.st<-stack()
names <- list.files(path=file.path(wd, '/Seq_wks_negW_2050'), full.names=TRUE)
#now make the stacK:
NM_seq_2050.st <- stack(names[1:20]) 
#get mean across years
names(NM_seq_2050.st)
meanTI_seq_2050<-mean(NM_seq_2050.st)


# make a stack of NicheMapper yearly predictions - 2050
NM_seq_2070.st<-stack()
names <- list.files(path=file.path(wd, '/Seq_wks_negW_2070'), full.names=TRUE)
#now make the stacK:
NM_seq_2070.st <- stack(names[1:20]) 
#get mean across years
names(NM_seq_2070.st)
meanTI_seq_2070<-mean(NM_seq_2070.st)

################Again plot sliding window years
gen<-5
seq_wks_negW.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM_seq.st[[st:b]]
  gen_suit<-sum(layers)/gen
   # plot(gen_suit)
  seq_wks_negW.st<-stack( seq_wks_negW.st,gen_suit)
}

gen<-5
seq_wks_negW_2050.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM_seq_2050.st[[st:b]]
  gen_suit<-sum(layers)/gen
  # plot(gen_suit)
  seq_wks_negW_2050.st<-stack( seq_wks_negW_2050.st,gen_suit)
}
plot(seq_wks_negW_2050.st)

gen<-5
seq_wks_negW_2070.st<-stack()
for(b in gen:20){
  st<-b-(gen-1)
  end<-b
  layers<-NM_seq_2070.st[[st:b]]
  gen_suit<-sum(layers)/gen
  # plot(gen_suit)
  seq_wks_negW_2070.st<-stack( seq_wks_negW_2070.st,gen_suit)
}
plot(seq_wks_negW_2070.st)


