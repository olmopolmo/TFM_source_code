library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(OmnipathR)
library(scales)
library(stringr)
library(ggplot2)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
samples_id <- list.files(path = "RDS/frontier/Cancer/")
samples_id<- samples_id[c(-5,-19)]
df <- readRDS("RDS/patients_metatdata.RDS")
df <- df[c(-5,-19:-20),]
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Pr <- unique(df$sample[which(df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD] )
HD <- HOT[HOT %in% DUTRE]
TF<- readRDS("TF_acts.rds")
triunvirates_DUTRENEO <- readRDS("Results_Correlation/TF_correlation_scores/TF_Pair_ListDUTRENEO.RDS")
triunvirates_Standar <- readRDS("Results_Correlation/TF_correlation_scores/TF_Pair_ListStandar.RDS")

DUTRENEO_df<- readRDS(paste0("Results_Correlation/TF_correlation_scores/BigDF_DUTRENEO.rds"))
Standar_df <- readRDS(paste0("Results_Correlation/TF_correlation_scores/BigDF_Standar.rds"))

for(p in triunvirates_DUTRENEO){
  t <- strsplit(x = p,split = "[-_]")[[1]][1]
  p1 <- strsplit(x = p,split = "[-_]")[[1]][2]
  p2 <- strsplit(x = p,split = "[-_]")[[1]][3]
  for (samp in HOT[HOT %in% DUTRE]) {
    SUB_TF <- TF[which(TF$sample==samp),]
    rownames(SUB_TF) <- SUB_TF$spotid
    SUB_TF <- SUB_TF[,-c(1,length(colnames(SUB_TF)))]
    pvp <- readRDS(paste0("RDS/seurat_v5_SCT/",samp,".rds"))
    pvp <- AddMetaData(pvp,metadata = SUB_TF)
    
    border <- readRDS(paste0("RDS/frontier/Cancer/",samp,".rds"))
    border_TF <-  SUB_TF[which(rownames(SUB_TF) %in% colnames(border)),]
    border <- AddMetaData(border,border_TF)
    
    whole_view <- SpatialFeaturePlot(pvp,c(t,p1,p2),ncol = 3,pt.size.factor = 1.5)
    border_view <- SpatialFeaturePlot(border,c(t,p1,p2),ncol = 3,pt.size.factor = 1.5)
    if(samp %in% Nr){
      ggsave(paste0("Results_Correlation/TF_correlation_scores/Plots/Standar/",p,"/No_",p,"_",samp,"_whole.png"),whole_view,height = 8,width = 24)
      ggsave(paste0("Results_Correlation/TF_correlation_scores/Plots/Standar/",p,"/No_",p,"_",samp,"_Border.png"),border_view,height = 8,width = 24)
    }else{
      ggsave(paste0("Results_Correlation/TF_correlation_scores/Plots/Standar/",p,"/Complete_",p,"_",samp,"_whole.png"),whole_view,height = 8,width = 24)
      ggsave(paste0("Results_Correlation/TF_correlation_scores/Plots/Standar/",p,"/Complete_",p,"_",samp,"_Border.png"),border_view,height = 8,width = 24)
      
    }
  }
}


SpatialFeaturePlot(border,c("ABL1","MDK","SDC3"),ncol = 3,pt.size.factor = 1.5)


