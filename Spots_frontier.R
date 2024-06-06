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
ids <- c()
samples_id <- samples_id[-c(5,19,20)]
for(i in samples_id){1
  pvp<- readRDS(file = paste0("RDS/seurat_v5_SCT/",i))
  
  ids <- c(ids, paste0(i,"_",colnames(pvp[,which(pvp$compartment_NewBuff=="Cancer")])))
}
gsub(pattern = ".rds",replacement = "",x = ids)
saveRDS(object = ids,"RDS/frontier_SPOT_IDS_WholeCancer.RDS")
