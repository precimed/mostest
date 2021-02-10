#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)

Path <- ...
setwd(paste0(Path,"MOST/DM/"))

Batch="39k"
Date="051119"

#################
# Prep DM GWAS  #
#################

############### Add pop comps and genetic batch to demog file
Demog <- fread(paste0(Path,"UKB/covariates/UKB",Batch,"_demogShort_",Date,".txt"),data.table=F)
PopComp <- fread(paste0(Path,"UKB/covariates/UKB500k_popComps_",Date,".txt"),data.table=F)[,c("ID",paste0("C",1:20))]
DM <- merge(Demog,PopComp, by = "ID")

############### Clean up and remove missings, non-caucasians, and brain Dx
DM$Sex <- as.data.frame(model.matrix(~ Sex - 1, data=DM))[,2]
DM <- DM[!duplicated(DM$ID),]
DM <- droplevels(DM)

Cauc <- fread(paste0(Path,"UKB/subject_lists/UKB500k_caucasians_",Date,".txt"),data.table=F)
DM <- DM[DM$ID %in% Cauc$V2,]
BrainDx <- fread(paste0(Path,"UKB/subject_lists/UKB500k_brainDxSubs_",Date,".txt"),data.table=F)
DM <- DM[!(DM$ID %in% BrainDx$V2),]
AvGen <- fread(paste0(Path,"UKB/subject_lists/UKB",Batch,"_genSubList_",Date,".txt"),data.table=F)

overlap <- intersect(DM$ID,AvGen$V1)
p <- match(overlap,DM$ID)
DM <- DM[p,]

############## Remove bad scans based on Euler numbers
Euler <- fread(paste0(Path,"UKB/covariates/UKB",Batch,"_Euler_",Date,".csv"),data.table=F)
colnames(Euler)[1] <- "ID"
DM <- merge(DM,Euler, by = "ID")
DM <- DM[complete.cases(DM),]

#### Write list of IDs with complete data for relatedness calculation
write.table(DM[,c("ID")],file=paste0(Path,"UKB/subject_lists/UKB",Batch,"_complete4relatedness_",Date,".txt"), sep = "\t", quote=FALSE, row.names=FALSE) 

Outliers <- c()
for(i in 1:length(unique(DM$Scanner))){
  current <- DM[DM$Scanner==unique(DM$Scanner)[i],]
  current$euler_lh_resid <- scale(lm(euler_lh ~ Sex + Age, data = current)$residuals)
  current$euler_rh_resid <- scale(lm(euler_rh ~ Sex + Age, data = current)$residuals)
  Outliers <- c(Outliers,current$ID[which(current$euler_lh_resid < -3 | current$euler_rh_resid < -3)])
}

DM <- DM[!(DM$ID %in% Outliers),]
DM$Euler <- DM$euler_lh + DM$euler_rh

################### Select base covariates and write to file
write.table(DM[,c("ID","Age","Sex","Scanner","Euler",paste0("C",1:20))],file=paste0(Path,"MOST/DM/UKB",Batch,"_baseGWAS_",Date,".txt"), sep = "\t", quote=FALSE, row.names=FALSE) 

# #### Make MATLAB friendly version
ScannerDums <- as.data.frame(model.matrix(~ Scanner - 1, data=DM))[,-1]
DM2 <- cbind(DM,ScannerDums)
write.table(DM2[,c("ID","Age","Sex",colnames(ScannerDums),"Euler",paste0("C",1:10))],file=paste0(Path,"MOST/DM/UKB",Batch,"_DummyBaseGWAS_",Date,".txt"), sep = "\t", quote=FALSE, row.names=FALSE) 

