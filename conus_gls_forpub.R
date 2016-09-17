##Conus data set
##predicting log max egg diameter and log max number of eggs per capsule; running separate models for different types of data. 
##Note: these analyses use log-transformed data
##1.) morphometry model (log.SL.max, rep.LW.ratio, rep.PMD)
##2) predation model (mol.pred, verm.pred, pisc.pred)
##3.) substratum model (sand, reef, rock, rubble, mud)
##4.) depth model (min.depth, max.depth, intertidal)
##5.) geography model (maxLat, minLong, maxLong)
##6.) raw climate model (sstmean.MEAN, chlomean.MEAN, salinity.MEAN)
##7.) climate variability model (sstmean.STD, chlomean.STD, salinity.STD, ph.STD)

setwd("~/Documents/NESCent Project/Main conus files/Latest versions") #set working directory

library(ape) #required for pruning tree
library (MuMIn) # required for multimodel inference
library(arm) # required for standardizing effect sizes
library (nlme) #required for gls
library (phylolm) #required for models including phylogeny
library(HH) #required for vif

conus_56spp=read.csv("conus_56spp_forpub.csv", row.names = 1) #attach full Conus dataset
attach(conus_56spp)


############Produce pruned Conus tree############
conetree <- read.nexus("conetree.nex")

remove <- conetree$tip.label[-which(conetree$tip.label %in% rownames(conus_56spp))] #find mismatched tips

pru.conetree <- drop.tip(conetree, tip=remove) #drop names not in conus_56spp

#scale all continuous numerical variables used:
log.max.ed <- scale(log.max.ed)
log.max.cap <- scale(log.max.cap)
resid <- scale(resid)
rep.LW.ratio <- scale(rep.LW.ratio)
log.SL.max <- scale(log.SL.max)
rep.PMD <- scale(rep.PMD)
min.depth <- scale(min.depth)
max.depth <- scale(max.depth)
minLong <- scale(minLong)
maxLong <- scale(maxLong)
maxLat <- scale(maxLat)
salinity.MEAN <- scale(salinity.MEAN)
sstmean.MEAN <- scale(sstmean.MEAN)
chlomean.MEAN <- scale(chlomean.MEAN)
chlomean.STD <- scale(chlomean.STD)
salinity.STD <- scale(salinity.STD)
sstmean.STD <- scale(sstmean.STD)
ph.STD <- scale(ph.STD)
silicate.MEAN <- scale(silicate.MEAN)
nitrate.MEAN <- scale(nitrate.MEAN)
phosphate.MEAN <- scale(phosphate.MEAN) 
calcite.MEAN <- scale(calcite.MEAN)

##########gls models for predicting log maximum egg diameter##########

#####morphometric model with no phylogeny#####
egg.diam.morph<- gls(log.max.ed ~      
                       log.SL.max + rep.LW.ratio + rep.PMD,
                     method = "ML")
AICc(egg.diam.morph) #get AIC of morph model
summary(egg.diam.morph) #summary of morph model 

dr_egg.diam.morph <- dredge(egg.diam.morph,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]])) #dredge

#plot(dr_egg.diam.morph)
dr_egg.diam.morph #shows top models
modelaverage_egg.diam.morph<-model.avg(get.models(dr_egg.diam.morph, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.morph #summary effects from top models
summary(modelaverage_egg.diam.morph) #show best models with delta <2


#####predation model with no phylogeny#####
egg.diam.pred<- gls(log.max.ed ~
                  mol.pred + pisc.pred + verm.pred, 
                  method = "ML")
AICc(egg.diam.pred) #get AIC of predation model
summary(egg.diam.pred) #summary of predation model

dr_egg.diam.pred <- dredge(egg.diam.pred,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.pred)
dr_egg.diam.pred #shows top models
modelaverage_egg.diam.pred<-model.avg(get.models(dr_egg.diam.pred, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.pred #summary effects from top models
summary(modelaverage_egg.diam.pred) #show best models with delta <2


#####substratum model with no phylogeny#####
egg.diam.sub<- gls(log.max.ed ~      
                     sand + reef + rock + rubble + mud,
                   method = "ML")
AICc(egg.diam.sub) #get AIC of substratum model
summary(egg.diam.sub) #summary of substratum model

dr_egg.diam.sub <- dredge(egg.diam.sub,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.sub)
dr_egg.diam.sub #shows top models
modelaverage_egg.diam.sub<-model.avg(get.models(dr_egg.diam.sub, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.sub #summary effects from top models
summary(modelaverage_egg.diam.sub) #show best models with delta <2


#####depth model with no phylogeny#####
egg.diam.depth<- gls(log.max.ed ~      
                       min.depth + max.depth + intertidal,
                     method = "ML")
AICc(egg.diam.depth) #get AIC of depth model
summary(egg.diam.depth) #summary of depth model

dr_egg.diam.depth <- dredge(egg.diam.depth,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.depth)
dr_egg.diam.depth #shows top models
modelaverage_egg.diam.depth<-model.avg(get.models(dr_egg.diam.depth, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.depth #summary effects from top models
summary(modelaverage_egg.diam.depth) #show best models with delta <2


#####geography model with no phylogeny#####
egg.diam.geog<- gls(log.max.ed ~      
                      maxLat + minLong + maxLong + Atlantic, #removed area.km2
                    method = "ML")
AICc(egg.diam.geog) #get AIC of geography model
summary(egg.diam.geog) #summary of geography model

dr_egg.diam.geog <- dredge(egg.diam.geog,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.geog)
dr_egg.diam.geog #shows top models
modelaverage_egg.diam.geog<-model.avg(get.models(dr_egg.diam.geog, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.geog #summary effects from top models
summary(modelaverage_egg.diam.geog) #show best models with delta <2


#####raw climate model with no phylogeny#####
egg.diam.rawclim<- gls(log.max.ed ~    
                         chlomean.MEAN + salinity.MEAN + sstmean.MEAN, #dissox.MEAN high neg corr w/ sstmean.MEAN; same for ph.MEAN w/chlomean.MEAN
                       method = "ML")
AICc(egg.diam.rawclim) #get AIC of raw clim model
summary(egg.diam.rawclim) #summary of raw clim model

dr_egg.diam.rawclim <- dredge(egg.diam.rawclim,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.rawclim)
dr_egg.diam.rawclim #shows top models
modelaverage_egg.diam.rawclim<-model.avg(get.models(dr_egg.diam.rawclim, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.rawclim #summary effects from top models
summary(modelaverage_egg.diam.rawclim) #show best models with delta <2


#####climate variability model with no phylogeny#####
egg.diam.climvar<- gls(log.max.ed ~      
                         chlomean.STD + salinity.STD + sstmean.STD + ph.STD, #dissox.STD high corr w/sstmean.STD
                       method = "ML")
AICc(egg.diam.climvar) #get AIC of climvar model
summary(egg.diam.climvar) #summary of climvar model

dr_egg.diam.climvar <- dredge(egg.diam.climvar,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.climvar)
dr_egg.diam.climvar #shows top models
modelaverage_egg.diam.climvar<-model.avg(get.models(dr_egg.diam.climvar, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.climvar #summary effects from top models
summary(modelaverage_egg.diam.climvar) #show best models with delta <2

#####nutrient model with no phylogeny#####
egg.diam.nuts<- gls(log.max.ed ~      
                         calcite.MEAN + nitrate.MEAN + phosphate.MEAN + silicate.MEAN,
                       method = "ML")
AICc(egg.diam.nuts) #get AIC of nutrient model
summary(egg.diam.nuts) #summary of nutrient model

dr_egg.diam.nuts <- dredge(egg.diam.nuts,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.nuts)
dr_egg.diam.nuts #shows top models
modelaverage_egg.diam.nuts<-model.avg(get.models(dr_egg.diam.nuts, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.nuts #summary effects from top models
summary(modelaverage_egg.diam.nuts) #show best models with delta <2

#####PCA climate model with no phylogeny#####
egg.diam.pcaclim<- gls(log.max.ed ~      
                         Factor1 + Factor2 + Factor3 + Factor4 + Factor5,
                       method = "ML")
AICc(egg.diam.pcaclim) #get AIC of PCA clim model
summary(egg.diam.pcaclim) #summary of PCA clim model

dr_egg.diam.pcaclim <- dredge(egg.diam.pcaclim,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_egg.diam.pcaclim)
dr_egg.diam.pcaclim #shows top models
modelaverage_egg.diam.pcaclim<-model.avg(get.models(dr_egg.diam.pcaclim, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.pcaclim #summary effects from top models
summary(modelaverage_egg.diam.pcaclim) #show best models with delta <2

vif(dr_egg.diam.pcaclim) #all equal to 1

#####full egg diam model with no phylogeny#####
egg.diam.full <- gls(log.max.ed ~      
                       rep.LW.ratio + intertidal + maxLong + salinity.MEAN + Atlantic +
                       chlomean.STD + salinity.STD + silicate.MEAN + nitrate.MEAN,
                     method = "ML")
AICc(egg.diam.full) #get AIC of full model
summary(egg.diam.full) #summary of full model

dr_egg.diam.full <- dredge(egg.diam.full,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
#plot(dr_egg.diam.full)
#dr_egg.diam.full #shows top models
head(dr_egg.diam.full, n = 40L)
modelaverage_egg.diam.full<-model.avg(get.models(dr_egg.diam.full, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_egg.diam.full #summary effects from top models
summary(modelaverage_egg.diam.full) #show full models with delta <2


#####full egg diam model WITH phylogeny#####
conus_56spp$log.max.ed <- log.max.ed
conus_56spp$rep.LW.ratio <- rep.LW.ratio
conus_56spp$minLong <- minLong
conus_56spp$salinity.MEAN <- salinity.MEAN
conus_56spp$chlomean.STD <- chlomean.STD
conus_56spp$salinity.STD <- salinity.STD
conus_56spp$silicate.MEAN <- silicate.MEAN
conus_56spp$nitrate.MEAN <- nitrate.MEAN
egg.diam.phyl <- phylolm(log.max.ed ~      
                           rep.LW.ratio + intertidal + maxLong + salinity.MEAN + Atlantic +
                           chlomean.STD + salinity.STD + silicate.MEAN + nitrate.MEAN, data=conus_56spp, 
                         phy=pru.conetree,
                         model="lambda")

egg.diam.phyl
summary(egg.diam.phyl)


##########gls models for predicting log max number of eggs per capsule##########

#####morphometric model with no phylogeny#####
logmax.cap.morph <- gls(log.max.cap ~      
                          log.SL.max + rep.LW.ratio + rep.PMD,
                        method = "ML")
AICc(logmax.cap.morph) #get AIC of morph model
summary(logmax.cap.morph) #summary of morph model 

dr_logmax.cap.morph <- dredge(logmax.cap.morph,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]])) #dredge

#plot(dr_logmax.cap.morph)
dr_logmax.cap.morph #shows top models
modelaverage_logmax.cap.morph<-model.avg(get.models(dr_logmax.cap.morph, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.morph #summary effects from top models
summary(modelaverage_logmax.cap.morph) #show best models with delta <2


#####predation model with no phylogeny#####
logmax.cap.pred<- gls (log.max.cap ~
                         mol.pred + pisc.pred + verm.pred,
                       method = "ML")
AICc(logmax.cap.pred) #get AIC of predation model
summary(logmax.cap.pred) #summary of predation model

dr_logmax.cap.pred <- dredge(logmax.cap.pred,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.pred)
dr_logmax.cap.pred #shows top models
modelaverage_logmax.cap.pred<-model.avg(get.models(dr_logmax.cap.pred, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.pred #summary effects from top models
summary(modelaverage_logmax.cap.pred) #show best models with delta <2


#####substratum model with no phylogeny#####
logmax.cap.sub<- gls(log.max.cap ~      
                       sand + reef + rock + rubble + mud,
                     method = "ML")
AICc(logmax.cap.sub) #get AIC of substratum model
summary(logmax.cap.sub) #summary of substratum model

dr_logmax.cap.sub <- dredge(logmax.cap.sub,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.sub)
dr_logmax.cap.sub #shows top models
modelaverage_logmax.cap.sub<-model.avg(get.models(dr_logmax.cap.sub, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.sub #summary effects from top models
summary(modelaverage_logmax.cap.sub) #show best models with delta <2


#####depth model with no phylogeny#####
logmax.cap.depth<- gls(log.max.cap ~      
                         min.depth + max.depth + intertidal,
                       method = "ML")
AICc(logmax.cap.depth) #get AIC of depth model
summary(logmax.cap.depth) #summary of depth model

dr_logmax.cap.depth <- dredge(logmax.cap.depth,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.depth)
dr_logmax.cap.depth #shows top models
modelaverage_logmax.cap.depth<-model.avg(get.models(dr_logmax.cap.depth, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.depth #summary effects from top models
summary(modelaverage_logmax.cap.depth) #show best models with delta <2


#####geography model with no phylogeny#####
logmax.cap.geog<- gls(log.max.cap ~      
                        maxLat + minLong + maxLong + Atlantic, #took out area.km2
                      method = "ML")
AICc(logmax.cap.geog) #get AIC of geography model
summary(logmax.cap.geog) #summary of geography model

dr_logmax.cap.geog <- dredge(logmax.cap.geog,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.geog)
dr_logmax.cap.geog #shows top models
modelaverage_logmax.cap.geog<-model.avg(get.models(dr_logmax.cap.geog, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.geog #summary effects from top models
summary(modelaverage_logmax.cap.geog) #show best models with delta <2

#####raw climate model with no phylogeny#####
logmax.cap.rawclim<- gls(log.max.cap ~      
                           chlomean.MEAN + salinity.MEAN + sstmean.MEAN,
                       method = "ML")
AICc(logmax.cap.rawclim) #get AIC of raw clim model
summary(logmax.cap.rawclim) #summary of raw clim model

dr_logmax.cap.rawclim <- dredge(logmax.cap.rawclim,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.rawclim)
dr_logmax.cap.rawclim #shows top models
modelaverage_logmax.cap.rawclim<-model.avg(get.models(dr_logmax.cap.rawclim, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.rawclim #summary effects from top models
summary(modelaverage_logmax.cap.rawclim) #show best models with delta <2


#####climate variability model with no phylogeny#####
logmax.cap.climvar<- gls(log.max.cap ~      
                         chlomean.STD + ph.STD + salinity.STD + sstmean.STD,
                       method = "ML")
AICc(logmax.cap.climvar) #get AIC of climvar model
summary(logmax.cap.climvar) #summary of climvar model

dr_logmax.cap.climvar <- dredge(logmax.cap.climvar,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.climvar)
dr_logmax.cap.climvar #shows top models
modelaverage_logmax.cap.climvar<-model.avg(get.models(dr_logmax.cap.climvar, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.climvar #summary effects from top models
summary(modelaverage_logmax.cap.climvar) #show best models with delta <2

#####nutrient model with no phylogeny#####
logmax.cap.nuts<- gls(log.max.cap ~      
                      calcite.MEAN + nitrate.MEAN + silicate.MEAN + phosphate.MEAN,
                    method = "ML")
AICc(logmax.cap.nuts) #get AIC of nutrient model
summary(logmax.cap.nuts) #summary of nutrient model

dr_logmax.cap.nuts <- dredge(logmax.cap.nuts,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_logmax.cap.nuts)
dr_logmax.cap.nuts #shows top models
modelaverage_logmax.cap.nuts<-model.avg(get.models(dr_logmax.cap.nuts, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.nuts #summary effects from top models
summary(modelaverage_logmax.cap.nuts) #show best models with delta <2


#####full max eggs per cap model with no phylogeny#####
logmax.cap.full <- gls(log.max.cap ~      
                         log.SL.max + maxLat + Atlantic + salinity.MEAN + chlomean.STD + 
                         salinity.STD + calcite.MEAN + phosphate.MEAN,
                       method = "ML")
AICc(logmax.cap.full) #get AIC of full model
summary(logmax.cap.full) #summary of full model

dr_logmax.cap.full <- dredge(logmax.cap.full,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
#plot(dr_logmax.cap.full)
dr_logmax.cap.full #shows top models
#head(dr_logmax.cap.full, n = 40L) #prints first 40 lines of output so that all delta < 2 are visible
modelaverage_logmax.cap.full<-model.avg(get.models(dr_logmax.cap.full, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_logmax.cap.full #summary effects from top models
summary(modelaverage_logmax.cap.full) #show best models with delta <2


#####full max eggs per cap model WITH phylogeny#####
conus_56spp$log.max.cap <- log.max.cap
conus_56spp$log.SL.max <- log.SL.max
conus_56spp$minLong <- minLong
conus_56spp$salinity.MEAN <- salinity.MEAN
conus_56spp$chlomean.STD <- chlomean.STD
conus_56spp$salinity.STD <- salinity.STD
conus_56spp$calcite.MEAN <- calcite.MEAN
conus_56spp$phosphate.MEAN <- phosphate.MEAN

logmax.cap.phyl <- phylolm(log.max.cap ~      
                             log.SL.max + maxLat + Atlantic + salinity.MEAN + chlomean.STD + 
                             salinity.STD + calcite.MEAN + phosphate.MEAN, data=conus_56spp, phy=pru.conetree,
                        model="lambda")
logmax.cap.phyl
summary(logmax.cap.phyl)


##########gls models for predicting residuals of egg diam vs. max eggs cap##########

#####morphometric model with no phylogeny#####
resid.morph<- gls(resid ~      
                    log.SL.max + rep.LW.ratio + rep.PMD,
                  method = "ML")
AICc(resid.morph) #get AIC of morph model
summary(resid.morph) #summary of morph model 

dr_resid.morph <- dredge(resid.morph,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]])) #dredge

#plot(dr_resid.morph)
dr_resid.morph #shows top models
modelaverage_resid.morph<-model.avg(get.models(dr_resid.morph, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.morph #summary effects from top models
summary(modelaverage_resid.morph) #show best models with delta <2


#####predation model with no phylogeny#####
resid.pred<- gls(resid ~
                   mol.pred + pisc.pred + verm.pred, 
                 method = "ML")
AICc(resid.pred) #get AIC of predation model
summary(resid.pred) #summary of predation model

dr_resid.pred <- dredge(resid.pred,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.pred)
dr_resid.pred #shows top models
modelaverage_resid.pred<-model.avg(get.models(dr_resid.pred, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.pred #summary effects from top models
summary(modelaverage_resid.pred) #show best models with delta <2


#####substratum model with no phylogeny#####
resid.sub<- gls(resid ~      
                  sand + reef + rock + rubble + mud,
                method = "ML")
AICc(resid.sub) #get AIC of substratum model
summary(resid.sub) #summary of substratum model

dr_resid.sub <- dredge(resid.sub,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.sub)
dr_resid.sub #shows top models
modelaverage_resid.sub<-model.avg(get.models(dr_resid.sub, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.sub #summary effects from top models
summary(modelaverage_resid.sub) #show best models with delta <2


#####depth model with no phylogeny#####
resid.depth<- gls(resid ~      
                    min.depth + max.depth + intertidal,
                  method = "ML")
AICc(resid.depth) #get AIC of depth model
summary(resid.depth) #summary of depth model

dr_resid.depth <- dredge(resid.depth,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.depth)
dr_resid.depth #shows top models
modelaverage_resid.depth<-model.avg(get.models(dr_resid.depth, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.depth #summary effects from top models
summary(modelaverage_resid.depth) #show best models with delta <2


#####geography model with no phylogeny#####
resid.geog<- gls(resid ~      
                   maxLat + minLong + maxLong + Atlantic, #removed area.km2
                 method = "ML")
AICc(resid.geog) #get AIC of geography model
summary(resid.geog) #summary of geography model

dr_resid.geog <- dredge(resid.geog,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.geog)
dr_resid.geog #shows top models
modelaverage_resid.geog<-model.avg(get.models(dr_resid.geog, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.geog #summary effects from top models
summary(modelaverage_resid.geog) #show best models with delta <2


#####raw climate model with no phylogeny#####
resid.rawclim<- gls(resid ~    
                      chlomean.MEAN + salinity.MEAN + sstmean.MEAN, #dissox.MEAN high neg corr w/ sstmean.MEAN; same for ph.MEAN w/chlomean.MEAN
                    method = "ML")
AICc(resid.rawclim) #get AIC of raw clim model
summary(resid.rawclim) #summary of raw clim model

dr_resid.rawclim <- dredge(resid.rawclim,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.rawclim)
dr_resid.rawclim #shows top models
modelaverage_resid.rawclim<-model.avg(get.models(dr_resid.rawclim, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.rawclim #summary effects from top models
summary(modelaverage_resid.rawclim) #show best models with delta <2


#####climate variability model with no phylogeny#####
resid.climvar<- gls(resid ~      
                      chlomean.STD + salinity.STD + sstmean.STD + ph.STD, #dissox.STD high corr w/sstmean.STD
                    method = "ML")
AICc(resid.climvar) #get AIC of climvar model
summary(resid.climvar) #summary of climvar model

dr_resid.climvar <- dredge(resid.climvar,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.climvar)
dr_resid.climvar #shows top models
modelaverage_resid.climvar<-model.avg(get.models(dr_resid.climvar, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.climvar #summary effects from top models
summary(modelaverage_resid.climvar) #show best models with delta <2

#####nutrient model with no phylogeny#####
resid.nuts<- gls(resid ~      
                   calcite.MEAN + nitrate.MEAN + phosphate.MEAN + silicate.MEAN,
                 method = "ML")
AICc(resid.nuts) #get AIC of nutrient model
summary(resid.nuts) #summary of nutrient model

dr_resid.nuts <- dredge(resid.nuts,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))

#plot(dr_resid.nuts)
dr_resid.nuts #shows top models
modelaverage_resid.nuts<-model.avg(get.models(dr_resid.nuts, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.nuts #summary effects from top models
summary(modelaverage_resid.nuts) #show best models with delta <2


#####full residual model with no phylogeny#####
resid.full <- gls(resid ~      
                    rep.LW.ratio + log.SL.max + mol.pred,
                  method = "ML")
AICc(resid.full) #get AIC of full model
summary(resid.full) #summary of full model

dr_resid.full <- dredge(resid.full,extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))
#plot(dr_resid.full)
dr_resid.full #shows top models
modelaverage_resid.full<-model.avg(get.models(dr_resid.full, subset= delta <2)) #retreives models with delta <2 of best model
modelaverage_resid.full #summary effects from top models
summary(modelaverage_resid.full) #show full models with delta <2


#####full residual model WITH phylogeny#####
conus_56spp$resid <- resid
conus_56spp$rep.LW.ratio <- rep.LW.ratio
conus_56spp$log.SL.max <- log.SL.max
resid.phyl <- phylolm(resid ~      
                        rep.LW.ratio + log.SL.max + mol.pred, data=conus_56spp,
                      phy=pru.conetree,
                      model="lambda")

resid.phyl
summary(resid.phyl)


#########Regression code##########
####To put Greek letters in: 
# for x-axis: plot(log.min.ed,log.max.cap, pch=19, col="blue", xlab = expression(paste("Log Min Egg Diameter (", mu, "m)")), ylab = "Log Max # Eggs Per Capsule")
# for y-axis: plot(rep.LW.ratio,log.max.ed, pch=19, col="blue", xlab = "Log Max # Eggs Per Capsule"), ylab = expression(paste("Log Max Egg Diameter (", mu, "m)"))

#plot for max egg diam regressed with max eggs per cap
plot(log.max.ed,log.max.cap, pch=21, bg="gray45", cex=1.5, xlab = expression(paste("Log"[10]*" Max Egg Diameter (", mu, "m)")), 
     ylab = expression(paste("Log"[10]*" Max # Eggs/Capsule")))
reg1<-lm(log.max.cap~log.max.ed)
summary (reg1)
abline(reg1, col = "black", lwd=1.5)
#legend("topright", c(expression(paste("Adj. R"^"2"," = xxx")), "p = xxx"), bty="o")


##### for log.max.ed ####
par(mfrow=c(2,2))

plot(chlomean.STD,log.max.ed, pch=21, bg="gray45", cex=1.5, xlab = expression(paste("StDev of Mean Chlorophyll (mg/","m"^"3",")")), 
     ylab = expression(paste("Log"[10]*" Max Egg Diameter (", mu, "m)")))#used darkorchid4 before
reg1<-lm(log.max.ed~chlomean.STD)
summary (reg1)
abline(reg1, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.028")), "p = 0.114"), bty="o")

plot(salinity.STD,log.max.ed, pch=21, bg="gray45", cex=1.5, xlab = "StDev of Mean Salinity (psu)", 
     ylab = expression(paste("Log"[10]*" Max Egg Diameter (", mu, "m)")))#used cornflowerblue before
reg2<-lm(log.max.ed~salinity.STD)
summary (reg2)
abline(reg2, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.072")), "p = 0.026"), bty="o")

plot(rep.LW.ratio,log.max.ed, pch=21, bg="gray45", cex=1.5, xlab = "Shell Length/Width Ratio", 
     ylab = expression(paste("Log"[10]*" Max Egg Diameter (", mu, "m)")))#used palegreen4 before
reg3<-lm(log.max.ed~rep.LW.ratio)
summary (reg3)
abline(reg3, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.144")), "p = 0.002"), bty="o")

##### for log.max.cap ####
par(mfrow=c(2,2))

plot(salinity.STD,log.max.cap, pch=21, bg="gray45", cex=1.5, xlab = "StDev of Mean Salinity (psu)", 
     ylab = expression("Log"[10]*" Max # Eggs/Capsule"))
reg5<-lm(log.max.cap~salinity.STD)
summary (reg5)
abline(reg5, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.137")), "p = 0.003"), bty="o")

plot(log.SL.max,log.max.cap, pch=21, bg="gray45", cex=1.5, xlab = expression("Log"[10]*" Max Shell Length (mm)"), 
     ylab = expression("Log"[10]*" Max # Eggs/Capsule"))
reg6<-lm(log.max.cap~log.SL.max)
summary (reg6)
abline(reg6, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.190")), "p = 0.0005"), bty="o")

fit <- glm(Atlantic~log.max.cap,data=conus_56spp,family=binomial)
summary(fit) # display results
plot(log.max.cap,Atlantic,xlab= expression("Log"[10]*" Max # Eggs/Capsule"),
     ylab="Atlantic Range Only",pch=21,cex=1.5,bg="gray45", yaxt="n")# used darkslateblue before; plot with body size on x-axis and survival (0 or 1) on y-axis
axis(2,at=Atlantic)
curve(predict(fit, data.frame(log.max.cap=x),type="response"),lwd=1.5,col="black", add=TRUE)

boxplot(log.max.cap~Atlantic, xlab="Atlantic Range Only", 
        ylab = expression("Log"[10]*" Max # Eggs/Capsule"))


##### for residuals ####
par(mfrow=c(1,2))

plot(rep.LW.ratio,resid, pch=21, bg="gray45", cex=1.5, xlab = "Shell Length/Width Ratio", 
     ylab = "Residuals: Egg Diam v. Egg Num")
reg7<-lm(resid~rep.LW.ratio)
summary (reg7)
abline(reg7, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.108")), "p = 0.008"), bty="o")

plot(log.SL.max,resid, pch=21, bg="gray45", cex=1.5, xlab = expression("Log"[10]*" Max Shell Length (mm)"), 
     ylab = "Residuals: Egg Diam v. Egg Num")
reg8<-lm(resid~log.SL.max)
summary (reg8)
abline(reg8, col = "black")
#legend("topright", c(expression(paste("Adj. R"^"2"," = 0.344")), "p = 1.202e-06"), bty="o")

##########Boxplot code##########
par(mfrow=c(1,2))
##For log max egg diam vs. crawling juv
boxplot(log.max.ed~crawl.juv, xlab="Nonplanktonic Hatchling", 
        ylab = expression(paste("Log"[10]*" Max Egg Diameter (", mu, "m)")))
#legend("topleft", c("0 = No", "1 = Yes"), cex=0.9, bty="n")

##For log max cap vs. crawling juv
boxplot(log.max.cap~crawl.juv, xlab="Nonplanktonic Hatchling", 
        ylab = expression(paste("Log"[10]*" Max # Eggs/Capsule")))

##For Factor 1 vs. crawling juv
boxplot(Factor1~crawl.juv, xlab="Hatching as Crawling Juvenile", ylab = "Factor1")
legend("topright", c("0 = No", "1 = Yes"), cex=0.9, bty="n")


