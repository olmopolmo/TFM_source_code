library(RColorBrewer)
library(magrittr)
library(spatstat)

library(scales)
# data manipulation
library(Matrix)
library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
df <- readRDS("RDS/patients_metatdata.RDS")
samples <- list.files(path = "RDS/seurat_v5_SCT/")
before_df <- data.frame()
after_df <- data.frame()
for(i in samples){
  du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
  du$compartment_NewBuff <- du$compartment
  # Find DE features 
  Idents(du) <- du$compartment
  de.markers <- FindMarkers(du,assay = "SCT", ident.1 = "Cancer", ident.2 = "TME")
  de.marker <- de.markers[de.markers$p_val_adj<5e-2,]
  TME.markers <- de.marker[which(de.marker$avg_log2FC<=-1),]
  Cancer.markers <- de.markers[which(de.markers$avg_log2FC>=1),]
  before <- as.data.frame(table(du$compartment))
  before_df <- rbind(before_df,before)
  # Build Buffer matrix
  buff_spots <- rownames(du@meta.data[which(du$compartment=="Buffer"),]) 
  buff <- as.data.frame(du@assays$SCT$data)
  buff <- buff[,which(colnames(buff) %in% buff_spots)]
  buff_tme <- as.data.frame(t(buff[which(rownames(buff) %in% rownames(TME.markers)),]))
  buff_cancer <- as.data.frame(t(buff[which(rownames(buff) %in% rownames(Cancer.markers)),]))
  buff_cancer$Exp_Can <- rowMeans(buff_cancer)
  buff_tme$Exp_TME <- rowMeans(buff_tme)
  buff <- cbind(buff_cancer,buff_tme)
  buff$Diff <- log2(buff$Exp_TME/buff$Exp_Can)# Positive values will be categoried as TME and negative as Cancer
  buff$id <- rownames(buff)
  buff <- buff[,c("id","Diff")]
  TME_id <- buff$id[buff$Dif>=1]
  du@meta.data$compartment_NewBuff[which(rownames(du@meta.data) %in% TME_id)] <- "TME"
  Cancer_id <- buff$id[buff$Diff<=-1]
  du@meta.data$compartment_NewBuff[which(rownames(du@meta.data) %in% Cancer_id)] <- "Cancer"
  saveRDS(du,file = paste0("RDS/seurat_v5_SCT/",i))
  after <- as.data.frame(table(du$compartment_NewBuff))
  after_df <- rbind(after_df,after)
}

# cancer_spots <- rownames(du@meta.data[which(du$compartment=="Cancer"),]) 
# TME_spots <- rownames(du@meta.data[which(du$compartment=="TME"),]) 
# SCT <- as.data.frame(du@assays$SCT$data)
# TME <- SCT[,which(colnames(SCT) %in% TME_spots)]
# TME<-as.data.frame(t(TME))
# cancer <- SCT[,which(colnames(SCT) %in% cancer_spots)]
# cancer<-as.data.frame(t(cancer))
# 
# log2(mean(cancer$DGCR6)/mean(TME$DGCR6))
# 
# buff <- buff[,which(colnames(buff) %in% buff_spots)]
# buff_tme <- as.data.frame(t(buff[which(rownames(buff) %in% rownames(TME.markers)),]))
# buff_cancer <- as.data.frame(t(buff[which(rownames(buff) %in% rownames(Cancer.markers)),]))
# buff_cancer$Exp_Can <- rowMeans(buff_cancer)
# buff_tme$Exp_TME <- rowMeans(buff_tme)
# buff <- cbind(buff_cancer,buff_tme)
# buff$Diff <- buff$Exp_TME/buff$Exp_Can # Positive values will be categoried as TME and negative as Cancer
# buff$id <- rownames(buff)
# buff <- buff[,c("id","Diff")]
