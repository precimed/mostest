#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)

Path <- ...
setwd(paste0(Path,"MOST/DM/"))

Batch="33k"
Date="230519"

######################
# Merge all measures #
######################

DM <- fread(paste0(Path,"MOST/DM/UKB",Batch,"_baseGWAS_",Date,".txt"),data.table=F)
Dict <- fread(paste0(Path,"MOST/DM/FSdict_ver53_",Date,".txt"),data.table=F)

### Add FS aseg to base
Aseg <- fread(paste0(Path,"UKB/phenotypes/FS/UKB",Batch,"_subcortical_",Date,".txt"),data.table=F)
DM <- merge(DM, Aseg,by="ID")

### Add Desikan-Killiany
measures <- c("area","thickness","volume")
for(i in 1:3){
  templ <- fread(paste0(Path,"UKB/phenotypes/FS/UKB",Batch,"_DK.lh.",measures[i],"_",Date,".txt"),data.table=F)
  tempr <- fread(paste0(Path,"UKB/phenotypes/FS/UKB",Batch,"_DK.rh.",measures[i],"_",Date,".txt"),data.table=F)
  temp <- merge(templ,tempr,by="ID")
  DM <- merge(DM,temp,by.x="ID",by.y="ID")
} 

DM$lh_TotalCortVol_volume <- rowSums(DM[,grep("_volume",colnames(templ),value = T)])
DM$rh_TotalCortVol_volume <- rowSums(DM[,grep("_volume",colnames(tempr),value = T)])

colnames(DM) <- gsub("-","_",colnames(DM))

##########################
# Regress out covariates #
##########################

#### Regress out covariates for ICV, area, thickness and volume globals
Globals <- Dict$ROI_name[grep("Summary",Dict$ROI_type)]
for(g in Globals){
  TempModel <- summary(lm(reformulate(termlabels = c('Age','Sex',"Scanner","Euler",paste0("C",1:20)), response = g),data = DM)) 
  DM[,g] <- TempModel$coefficients[1,1] + TempModel$residuals
} 

#### For all aseg, residualize
Subcort <- Dict$ROI_name[grep("SubcorticalVolume\\b",Dict$ROI_type)]
for(s in Subcort){
  TempModel <- summary(lm(reformulate(termlabels = c('Age','Sex',"Scanner","Euler","EstimatedTotalIntraCranialVol",paste0("C",1:20)), response = s),data=DM)) 
  DM[,s] <- TempModel$coefficients[1,1] + TempModel$residuals
}

#### For all area, residualize
Area <- Dict$ROI_name[grep("CorticalArea\\b",Dict$ROI_type)]
for(hemi in c("lh","rh")){
HemiArea <- Area[grep(hemi,Area)]
for(a in HemiArea){
  TempModel <- summary(lm(reformulate(termlabels = c('Age','Sex',"Scanner","Euler",paste0(hemi,"_WhiteSurfArea_area"),paste0("C",1:20)), response = a),data=DM)) 
  DM[,a] <- TempModel$coefficients[1,1] + TempModel$residuals 
  }}

#### For all thickness, residualize
Thick <- Dict$ROI_name[grep("CorticalThickness\\b",Dict$ROI_type)]
for(hemi in c("lh","rh")){
  HemiThick <- Thick[grep(hemi,Thick)]
for(t in HemiThick){
  TempModel <- summary(lm(reformulate(termlabels = c('Age','Sex',"Scanner","Euler",paste0(hemi,"_MeanThickness_thickness"),paste0("C",1:20)), response = t),data=DM)) 
  DM[,t] <- TempModel$coefficients[1,1] + TempModel$residuals
  }}

#### For all volume, residualize
Vol <- Dict$ROI_name[grep("CorticalVolume\\b",Dict$ROI_type)]
for(hemi in c("lh","rh")){
  HemiVol <- Thick[grep(hemi,Vol)]
for(v in HemiVol){
  TempModel <- summary(lm(reformulate(termlabels = c('Age','Sex',"Scanner","Euler",paste0(hemi,"_TotalCortVol_volume"),paste0("C",1:20)), response = v),data=DM)) 
  DM[,v] <- TempModel$coefficients[1,1] + TempModel$residuals 
  }}

##### Write final DM
Features <- c(Globals,Subcort,Area,Thick,Vol)
for(f in Features) {
  DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
}

DM$FID <- DM$ID -> DM$IID

write.table(DM[,c("FID","IID",Features)],file=paste0(Path,"MOST/DM/UKB",Batch,"_GWAS_qnorm_ecdf_",Date,".txt"),
            sep = "\t",quote=FALSE,row.names=FALSE)

#### Select subset PLINK files
GenPath <- ""
system(paste0("plink --bfile ",GenPath,"UKB",Batch,"_QCed_",Date," --keep --out UKB26502_QCed_230519")






write.table(DM[,Features],file=paste0(Path,"MOST/DM/UKB",nrow(DM),"_GWAS_FeaturesMOST_",Date,".txt"),
            sep = "\t",quote=FALSE,row.names=FALSE)
