library(STutility)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(OmnipathR)
library(scales)
# MISTy
library(mistyR)
library(future)
# data manipulation
library(Matrix)
library(tibble)
library(stringr)
library(dplyr)
library(purrr)
# normalization
library(sctransform)
# resource
library(decoupleR)
# plotting
library(ggplot2)
# setup parallel execution
plan(multisession)
# Seurat
library(Seurat)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
df <- readRDS("RDS/patients_metatdata.RDS")

HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])

STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])


Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Pr <- unique(df$sample[which(df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
n <- list.files(path = "RDS/seurat_v5_SCT/")
n <- n[c(5,19,20)]
for_misty <- readRDS("RDS/genes_LR.RDS")
for_misty <- c(for_misty$source_genesymbol,for_misty$target_genesymbol)
# WARINING: MISTY DOES NOT LIKE GENES THAT HAVE A "-" for the moment i will remove them
for_misty <- for_misty[-which(grepl(pattern = "-",x = for_misty))]
clusters <- c("Cancer")
for(clus in clusters){
  for(i in n){ # Resume from 14:30 LASTONE WAS 33
    du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
    DefaultAssay(du) <- "SCT"
    Idents(du) <- du$compartment_NewBuff
    rownames(du@tools$Staffli@meta.data) <- paste0(str_split(rownames(du@tools$Staffli@meta.data), "-", simplify = T)[,1], "-", str_split(colnames(du), "-", simplify = T)[1,2])
    du <- RegionNeighbours(du, id = clus,verbose = TRUE)
    du <- du[,which(du$nbs_Cancer==paste0("nbs_",clus)| du$nbs_Cancer==clus)]
    saveRDS(file = paste0("RDS/frontier/",clus,"/",i),object = du)
  }
}

