library(STutility)
library(tidyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(OmnipathR)
library(scales)
library(dplyr)
library(pheatmap)
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
n <- list.files(path = "RDS/seurat_v5_SCT/")
n <- n[-c(5,19,20)]
CS <- readRDS("Results_Correlation/RDS/Complete_vs_No_df_ColdStandar.RDS")
# Genes by whole tissue
for_misty <- c(CS$ligands,CS$receptors)
for_misty <- as.vector(for_misty[!duplicated(for_misty)])
compartments <- c("Cancer") 

for(c in compartments){ 
  for(i in n){
    du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
    DefaultAssay(du) <- "SCT"
    Idents(du) <- du$compartment_NewBuff
    du<- subset(du,ident=c(c))
    rownames(du@tools$Staffli@meta.data) <- paste0(str_split(rownames(du@tools$Staffli@meta.data), "-", simplify = T)[,1], "-", str_split(colnames(du), "-", simplify = T)[1,2])
    # Expression data
    expression <- GetAssayData(
      object = du,
      slot = "counts",
      assay = "SCT"
    )
    
    # Location data
    geometry <- GetTissueCoordinates(du, cols = c("row", "col"), scale = NULL
    )
    
    norm.data <- vst(as(expression, "dgCMatrix"), verbosity = 0)$y
    coverage <- rowSums(norm.data > 0) / ncol(norm.data)
    
    slide.markers <- for_misty
    slide.markers <- slide.markers[which(slide.markers %in% rownames(norm.data))]
    misty.views <- create_initial_view(t(norm.data)[, slide.markers] %>% as_tibble())
    misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 2)
    misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 3)
    misty.views <- misty.views %>% add_paraview(geometry,l=10)
    run_misty(misty.views, paste0("Results_MistyR/Correlation_Genes_Footprints/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")),num.threads = 1)
    # misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/",c,"/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")))
  }
}

du <- readRDS(paste0("RDS/frontier/Cancer/DU42.rds"))
du2 <- readRDS(paste0("RDS/frontier/Cancer/DU3.rds"))

SpatialFeaturePlot(du,features = c("PTN","PTPRS"))
SpatialFeaturePlot(du2,features = c("PTN","PTPRS"))

##############################################################
# NOW I BUILD A DATAFRAME WITH ALL THE DATA FROM EACH SAMPLE #
##############################################################
Cor_Foot <- list.files(path = "Results_MistyR/Correlation_Genes_Footprints/")

df <- readRDS("RDS/patients_metatdata.RDS")
n <- list.files(path = "RDS/seurat_v5_SCT/")
n <- n[-c(5,19,20)]
df <- df[-c(5,19,20),]
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Nr <- unique(df$sample[which(df$response=="NO")])
du <- c(HOT[1],COLD[1],STD[1],DUTRE[1],Cr[1],Nr[1])
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD] )
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")
MISTY <- data.frame()


for(w in Cor_Foot){
  sample <- gsub(x = w,pattern = "footprint_Cancer_",replacement = "")
  misty.results <- collect_results(paste0("Results_MistyR/Correlation_Genes_Footprints/",w))$importances.aggregated
  misty.results <- misty.results[-which(is.na(misty.results$Importance)),]  
  misty.results$sample <- sample
  misty.results$Type <- ""
  misty.results$Response <- ""
  misty.results$Treatment <- ""
  misty.results$Cluster <- "Whole"
  if(sample %in% HOT){
    misty.results$Type <- "HOT"
  } else{ misty.results$Type <- "COLD" }
  if(sample %in% STD){
    misty.results$Treatment <- "Standar"
  }else{ misty.results$Treatment <- "DUTRENEO" }
  if(sample %in% Nr){
    misty.results$Response <- "No"
  }else{misty.results$Response <- "Complete"}
  
  MISTY <- rbind(MISTY,misty.results)
}

# NOW ITERATE THROUGH ALL THE PAIRS THAT WERE OF INTEREST, MAKE A PLOT FOR EACH PAIR
# I SHOULD END UP WITH 10 FOR EACH GROUP (HOT,COLD,STD,DURENEO,NO, AND COMPLETE)
MISTY$Pairs <- paste0(MISTY$Predictor,"_",MISTY$Target)
saveRDS(object = MISTY,"Results_MistyR/MISTY_CORRELATION.RDS")
MISTY<- readRDS("Results_MistyR/MISTY_CORRELATION.RDS")
top <- gsub(pattern = "-",replacement = "_",x = CS$parejas)
top <- unique(top)
top <- top[which(top %in% MISTY$Pairs )]
ColdSTD <- MISTY[which(MISTY$Type=="COLD" & MISTY$Treatment=="Standar"),]
HotDUTRe <- MISTY[which(MISTY$Type=="HOT" & MISTY$Treatment=="DUTRENEO"),]
HotSTD <- MISTY[which(MISTY$Type=="HOT" & MISTY$Treatment=="Standar"),]
# library(ggpubr)
# for(p in top){
#   df_plot <- ColdSTD[which(ColdSTD$Pairs %in% p),]
#   # pair_HC<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Type)) +
#   #   geom_boxplot()+
#   #   facet_grid(rows = vars(Cluster)) +
#   #   labs(title = p)+
#   #   theme_classic()+
#   #   stat_compare_means(aes(group = Type))
#   # pair_SD<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Treatment)) +
#   #   geom_boxplot()+
#   #   facet_grid(rows = vars(Cluster)) +
#   #   labs(title = p)+
#   #   theme_classic()+
#   #   stat_compare_means(aes(group = Treatment))
#   pair_CN<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Response)) +
#     geom_boxplot()+
#     facet_grid(rows = vars()) +
#     labs(title = p)+
#     theme_classic()+
#     stat_compare_means(aes(group = Response))
#   # ggsave(filename = paste0("Results_MistyR/Plots/COLD_Standar/",p,"/",p,"_Cold_vs_Hot.png"),plot = pair_HC,height = 15,width = 10)
#   # ggsave(filename = paste0("Results_MistyR/Plots/COLD_Standar/",p,"_Standar_vs_Dutreneo.png"),plot = pair_SD,height = 15,width = 10)
#   ggsave(filename = paste0("Results_MistyR/Plots/COLD_Standar/",p,"_Complete_vs_No.png"),plot = pair_CN,height = 15,width = 10)
# }
vistas <- unique(MISTY$view)
sub_Misty<-MISTY[MISTY$Pairs %in% top,]
sub_Misty$Combination <- paste0(sub_Misty$sample,"_",sub_Misty$Type,"_",sub_Misty$Treatment,"_",sub_Misty$Response)

# SEPARATED BY ARMS
arms <- paste0(sub_Misty$Type,"_",sub_Misty$Treatment)
sub_Misty$Arms <- arms
arms <- unique(arms)
for(arm in arms){
  sub2_misty <- sub_Misty[which(sub_Misty$Arms==arm),]
  for(v in vistas){
    sub <- sub2_misty[sub2_misty$view == v,]
    # GOOD PLOT:
    heatmap_data <- sub %>%
      select(Pairs, Combination, Importance) %>%
      spread(key = Combination, value = Importance)
    rows <- as.factor(heatmap_data$Pairs)
    heatmap_data[,1] <- NULL
    heatmap_data <- as.data.frame(heatmap_data)
    rownames(heatmap_data)<- rows
    heatmap_data <- as.matrix(heatmap_data)
    labels <- strsplit(x = colnames(heatmap_data),split = "_")
    # Extract fourth elements
    Response <- sapply(labels, "[[", 4)
    
    # create a `df` with the samples grouped in the same way you want to show
    anno <- data.frame(Response=Response)
    
    # Sort columns based on Type
    sorted_columns <- order(anno$Response)
    heatmap_data <- heatmap_data[, sorted_columns]
    anno <- as.data.frame(anno[sorted_columns,] )  
    # set rownames so that anno and your data can be matched
    rownames(anno) <- colnames(heatmap_data)
    
    Heat<-pheatmap(heatmap_data,scale = "row", main = "Column Cluster Heatmap",annotation_col = anno,
                   cluster_cols = F)
    ggsave(filename = paste0("Results_MistyR/Plots/",arm,"/",v,"_Heatmap.png"),plot = Heat,height = 10,width = 15)
  }
}

