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
samples_id <- list.files(path = "RDS/seurat_v5_SCT/")
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
df <- readRDS(paste0("Results_Correlation/RDS/Complete_vs_No_df_HotDUTRENEO.RDS"))

parejas <- as.vector(unique(df$parejas))
for(p in parejas[29:75]){
  p1 <- strsplit(x = p,split = "-")[[1]][1]
  p2 <- strsplit(x = p,split = "-")[[1]][2]
  for (samp in HOT[HOT %in% DUTRE]) {
    pvp <- readRDS(paste0("RDS/seurat_v5_SCT/",samp,".rds"))
    border <- readRDS(paste0("RDS/frontier/Cancer/",samp,".rds"))
    whole_view <- SpatialFeaturePlot(pvp,c(p1,p2),ncol = 1,pt.size.factor = 1.5)
    border_view <- SpatialFeaturePlot(border,c(p1,p2),ncol = 1,pt.size.factor = 1.5)
    if(samp %in% Nr){
      ggsave(paste0("Results_Correlation/SpatialPlots/HotDUTRENEO/",p,"/No_",p,"_",samp,"_whole.png"),whole_view,height = 15,width = 8)
      ggsave(paste0("Results_Correlation/SpatialPlots/HotDUTRENEO/",p,"/No_",p,"_",samp,"_Border.png"),border_view,height = 15,width = 8)
    }else{
      ggsave(paste0("Results_Correlation/SpatialPlots/HotDUTRENEO/",p,"/Complete_",p,"_",samp,"_whole.png"),whole_view,height = 15,width = 8)
      ggsave(paste0("Results_Correlation/SpatialPlots/HotDUTRENEO/",p,"/Complete_",p,"_",samp,"_Border.png"),border_view,height = 15,width = 8)
      
    }
  }
}
