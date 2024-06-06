library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(liana)
library(tidyr)
setwd("/mnt/beegfs/mescobosa/DUTRENEO/LIANA")
samples <- list.files(path = "RDS/seurat_v5_SCT/")
names <- gsub(pattern = ".rds",replacement = "",x = samples)


# Load the needed libraries 
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(liana)
library(tidyr)
setwd("/mnt/beegfs/mescobosa/DUTRENEO/LIANA")
samples <- list.files(path = "../pipelines/seurat_v5_SCT")
names <- gsub(pattern = ".rds",replacement = "",x = list.files(path = "../pipelines/seurat_v5_SCT/"))

integrate <- merge(x = readRDS("../pipelines/seurat_v5_SCT/DU10.rds"),
                  y =c(readRDS("../pipelines/seurat_v5_SCT/DU13.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU14.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU15.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU17.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU18.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU19.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU2.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU21.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU22.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU23.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU24.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU25.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU26.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU27.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU28.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU29.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU3.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU30.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU31.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU33.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU34.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU35.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU37.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU38.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU39.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU40.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU41.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU42.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU5.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU50.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU8.rds"),
                       readRDS("../pipelines/seurat_v5_SCT/DU9.rds")),
                  project = "DUTRENEO_Integrated",add.cell.ids = names)
# Calculate the SCT assay
# integrate <- SCTransform(integrate, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)
# Normalize and scale the data (Harmony bases on the principal components to do the integration)
integrate <- integrate %>%
 Seurat::NormalizeData(verbose = FALSE) %>%
 FindVariableFeatures(selection.method = "SCT", nfeatures = 2000) %>%
 ScaleData(verbose = FALSE) %>%
 RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
# Add the column for the easier identification and later clustering
# NOTE: substitute the str_detect("X_") for the corresponding pattern in your data
for(sample in names){
 integrate$sample[colnames(integrate) %>% stringr::str_detect(sample)] <- sample
}
# Run harmony and then Clusterization on it.
integrate <- RunHarmony(object = integrate , plot_convergence = TRUE,group.by.vars = "sample")
integrate <- integrate %>%
 RunUMAP(reduction = "harmony", dims = 1:20) %>%
 FindNeighbors(reduction = "harmony", dims = 1:20) %>%
 FindClusters(resolution = 0.2)
# Save the object !!!
# saveRDS(file = paste0("Results/integrate.RDS"),object = integrate)
# integrate <- readRDS("Results_LIANA/RDS/Arms/integrate.RDS")
ids <- readRDS("Results_LIANA/RDS/Arms/frontier_SPOT_IDS_Cancer.RDS")
ids <- gsub(pattern = ".rds",replacement = "",x = ids)
integrate <- subset(integrate,cells = ids)
df <- readRDS("RDS/patients_metatdata.RDS")

HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])

STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])


Cr <- unique(df$sample[which(df$response=="COMPLETE" | df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD] )
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")
names <- c("Comp-No",
           "Comp-No",
           "Comp-No")
for(p in 1:3){
  p1 <- arm_samples[[p]][arm_samples[[p]] %in% Cr]
  p2 <- arm_samples[[p]][arm_samples[[p]] %in% Nr]

  n1 <- "Comp"
  n2 <- "No"
  # Anotar los nombres de las muestras junto al grupo que pertenecen
  integrate$sample[colnames(integrate) %>% stringr::str_detect(paste0(p1,collapse = "|"))] <- paste0(integrate$sample[colnames(integrate) %>% stringr::str_detect(paste0(p1,collapse = "|"))],"_",n1)
  integrate$sample[colnames(integrate) %>% stringr::str_detect(paste0(p2,collapse = "|"))] <- paste0(integrate$sample[colnames(integrate) %>% stringr::str_detect(paste0(p2,collapse = "|"))],"_",n2)
  # Anoto las muestras que son de cada grupo
  integrate$compartment_NewBuff[integrate$sample %>% stringr::str_detect(paste0(p1,collapse = "|"))] <- paste0(integrate$compartment_NewBuff[integrate$sample %>% stringr::str_detect(paste0(p1,collapse = "|"))],"_",n1)
  integrate$compartment_NewBuff[integrate$sample %>% stringr::str_detect(paste0(p2,collapse = "|"))] <- paste0(integrate$compartment_NewBuff[integrate$sample %>% stringr::str_detect(paste0(p2,collapse = "|"))],"_",n2)
  Idents(integrate) <-integrate$compartment_NewBuff
  ###################################################################
  # USE ALL SPOTS SURROUNDING CANCER BY ASIGNING TME TO ALL BUFFERS #
  ###################################################################
  integrate$compartment_NewBuff[which(integrate$compartment_NewBuff=="Buffer")]="TME"
  # Run liana
  liana_test <- liana_wrap(sce = integrate,resource = ,assay = "SCT")
  saveRDS(file = paste0("Results/LIANA_",arm[p],"_integrated.RDS"),object = liana_test)
  
  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate()
  #saveRDS(file = paste0("Results/LIANA_",arm[p],"_integrated.RDS"),object = liana_test) 
  gc()
  
  genes_du <- c()
  genes <- as.data.frame(AverageExpression(integrate,assays = "SCT"))
  genes$genes <- rownames(genes)
  genes_du <- genes$genes[which(genes[,1]>5 | genes[,2]>5 | genes[,3]>5 | genes[,4]>5 | genes[,5]>5 | genes[,6]>5)]
  
  genes_du <- unique(genes_du)
  big_df <- data.frame()
  plot1 <- liana_test %>%
    liana_dotplot(source_groups = c(paste0("Cancer_",n1),paste0("Cancer_",n2)), # CHANGE WHEN THE NEW DATA ARRIVES FOR THE NEW CLUSTERS
                  target_groups = c(paste0("TME_",n1),paste0("TME_",n2))) 
  # In this state the plot is not readable we have to filter
  df <- plot1$data
  # Change the symbol for clarity
  df$interaction <- sub(x = df$interaction,pattern = "-",replacement = "_")
  df$interaction <- sub(x = df$interaction,pattern = " _> ",replacement = "-")
  
  #df<- df[which(df$interaction %in% unique(CH$parejas)),]
  plot1$data <- df
  big_df <- rbind(big_df,df)
  big_df$interaction <- sub(x = big_df$interaction,pattern = " -> ",replacement = "-")
  
  # For loop in order to calculate the mean per target
  gp <- unique(big_df$interaction)
  targets <- unique(big_df$target)
  for (j in gp) {
    small_df <- big_df[which(big_df$interaction==j),]
    for(k in targets){
      big_df$magnitude[which(big_df$interaction==j & big_df$target == k)] <- mean(small_df$magnitude[which(small_df$target == k)])
    }
  }
  big_df$interaction<- as.factor(big_df$interaction)
  big_df$Group <- n1
  big_df$Group[big_df$source %>% stringr::str_detect(n2)] <- n2
  big_df <- big_df[!big_df$source_target==paste0("Cancer_",n1,"_TME_",n2),]
  big_df <- big_df[!big_df$source_target==paste0("Cancer_",n2,"_TME_",n1),]
  
  # REMOVE SINGLE INTERACTIONS
  interacciones <- unique(names(table(big_df$interaction)[which(table(big_df$interaction)==2)]))
  different_interactions <- c()
  filtered_df <- data.frame()
  for(i in interacciones){
    df <- big_df[which(big_df$interaction==i),]
    diferencia <- abs(df$magnitude[df$Group==n1] - df$magnitude[df$Group==n2])
    print(diferencia)
    if(diferencia>0.15){
      different_interactions <- c(different_interactions,i)
      filtered_df<-rbind(filtered_df,as.data.frame(df))
    }
  }
  saveRDS(object = liana_test, file = paste0("Results/liana_test_",arm[p],"_Cancer_",names[p],".RDS"))
  plot1 <- liana_test %>%
    liana_dotplot(source_groups = c(paste0("Cancer_",n1),paste0("Cancer_",n2)), 
                  target_groups = c(paste0("TME_",n1),paste0("TME_",n2))) 
  plot1$data <- filtered_df
  saveRDS(file = paste0("Results/filtered_df_",arm[p],"_",names[p],".RDS"),object = big_df)
}




