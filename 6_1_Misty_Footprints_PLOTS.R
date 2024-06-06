library(STutility)
library(tidyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
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


# THE POINT OF THIS SCRIPT IS TO TAKE THE TOP MARKERS FOR EACH CLSUTER AND CATEGORY
# AND FEED THEM AGAIN IN MISTY SO I GET THE VALUES FOR THE TOP MARKERS FOR ALL OF 
# THE CONDITIONS AND THEREFORE I CAN PLOT THEM. BECAUSE RIGHT NOW THE TOP MARKERS
# DO NOT HAVE ANY VALUE IN THE OTHER GROUPS.

# ¡¡¡ I WILL TAKE THE TOP 10 PAIRS OF GENES THAT HAVE THE HIGHEST MEAN VALUE IGNORING COLAGEN PAIRS !!!


#########
# WHOLE #
#########
misty_top <- readRDS("Results_MistyR/Misty_Top.rds")
misty_top$Views <- paste0(misty_top$view)
misty_top$Group <- as.factor(misty_top$Group)
misty_top <- misty_top %>% 
  filter(!grepl("COL", Pairs)) %>%  
  arrange(Group, desc(mean)) %>%  
  group_by(Group) %>% 
  slice(1:10)
##########
# BUFFER #
##########
misty_top_Buffer <- readRDS("Results_MistyR/Misty_Top_Buffer.rds")
misty_top_Buffer$Views <- paste0("Buffer_",misty_top_Buffer$view)
misty_top_Buffer$Group <- as.factor(misty_top_Buffer$Group)
misty_top_Buffer <- misty_top_Buffer %>% 
  filter(!grepl("COL", Pairs)) %>%  
  arrange(Group, desc(mean)) %>%  
  group_by(Group) %>% 
  slice(1:10)
#######
# TME #
#######
misty_top_TME <- readRDS("Results_MistyR/Misty_Top_TME.rds")
misty_top_TME$Views <- paste0("TME_",misty_top_TME$view)
misty_top_TME$Group <- as.factor(misty_top_TME$Group)
misty_top_TME <- misty_top_TME %>% 
  filter(!grepl("COL", Pairs)) %>%  
  arrange(Group, desc(mean)) %>%  
  group_by(Group) %>% 
  slice(1:10)
##########
# Cancer #
##########
misty_top_Cancer <- readRDS("Results_MistyR/Misty_Top_Cancer.rds")
misty_top_Cancer$Views <- paste0("Cancer_",misty_top_Cancer$view)
misty_top_Cancer$Group <- as.factor(misty_top_Cancer$Group)
misty_top_Cancer <- misty_top_Cancer %>% 
  filter(!grepl("COL", Pairs)) %>%  
  arrange(Group, desc(mean)) %>%  
  group_by(Group) %>% 
  slice(1:10)
# Variable with all the top gene pairs
top <- unique(c(misty_top_Buffer$Pairs,misty_top_Cancer$Pairs,misty_top_TME$Pairs,misty_top$Pairs))
top <- unique(misty_top$Pairs)
# NOW LETS RUN MISTY WITH THE NEW GENES WITHOUT FILTERING SO WE KEEP EVEN THE SMALL VALUES
# Genes by cluster
for_misty <- c(misty_top_Cancer$Predictor,misty_top_Cancer$Target,
                    misty_top_Buffer$Predictor,misty_top_Buffer$Target,
                    misty_top_TME$Predictor,misty_top_TME$Target)
for_misty <- for_misty[!duplicated(for_misty)]
# Genes by whole tissue
for_misty <- c(misty_top$Predictor,misty_top$Target)
for_misty <- for_misty[!duplicated(for_misty)]

compartments <- c("Cancer","Buffer","TME") 
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
    run_misty(misty.views, paste0("Results_MistyR/Latest_11_April_footprints/Whole/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")),num.threads = 1)
    # misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/",c,"/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")))
  }
}
#######################
# WHOLE TISSUE MISTYR #
#######################

for(i in n){
  du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
  DefaultAssay(du) <- "SCT"
  Idents(du) <- du$compartment_NewBuff
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
  run_misty(misty.views, paste0("Results_MistyR/Latest_11_April_footprints/Whole/footprint_",gsub(x = i,pattern = ".rds",replacement = "")),num.threads = 1)
  # misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/",c,"/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")))
}

##############################################################
# NOW I BUILD A DATAFRAME WITH ALL THE DATA FROM EACH SAMPLE #
##############################################################
all <- list.files(path = "Results_MistyR/Latest_11_April_footprints/All/")
buf <- list.files(path = "Results_MistyR/Latest_11_April_footprints/Buffer/")
can <- list.files(path = "Results_MistyR/Latest_11_April_footprints/Cancer/")
tme <- list.files(path = "Results_MistyR/Latest_11_April_footprints/TME/")
whole <- list.files(path = "Results_MistyR/Latest_11_April_footprints/Whole/")
df <- readRDS("RDS/patients_metatdata.RDS")
n <- list.files(path = "RDS/seurat_v5_SCT/")
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Nr <- unique(df$sample[which(df$response=="NO")])
du <- c(HOT[1],COLD[1],STD[1],DUTRE[1],Cr[1],Nr[1])
All <- list(HOT,COLD,STD,DUTRE,Cr,Nr)
All_names <- c("Hot","Cold","Standar","Dutreneo","Complete","No")
MISTY <- data.frame()

# for(b in buf){
#   sample <- gsub(x = b,pattern = "footprint_Buffer_",replacement = "")
#   misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/All/",b))$importances.aggregated
#   misty.results <- misty.results[-which(is.na(misty.results$Importance)),]  
#   
#   misty.results$Type <- ""
#   misty.results$Response <- ""
#   misty.results$Treatment <- ""
#   misty.results$Cluster <- "Buffer"
#   if(sample %in% HOT){
#     misty.results$Type <- "HOT"
#   } else{ misty.results$Type <- "COLD" }
#   if(sample %in% STD){
#     misty.results$Treatment <- "Standar"
#   }else{ misty.results$Treatment <- "DUTRENEO" }
#   if(sample %in% Nr){
#     misty.results$Response <- "No"
#   }else{misty.results$Response <- "Complete"}
#   
#   MISTY <- rbind(MISTY,misty.results)
# }
# for(c in can){
#   sample <- gsub(x = c,pattern = "footprint_Cancer_",replacement = "")
#   misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/All/",c))$importances.aggregated
#   misty.results <- misty.results[-which(is.na(misty.results$Importance)),]  
#   
#   misty.results$Type <- ""
#   misty.results$Response <- ""
#   misty.results$Treatment <- ""
#   misty.results$Cluster <- "Cancer"
#   if(sample %in% HOT){
#     misty.results$Type <- "HOT"
#   } else{ misty.results$Type <- "COLD" }
#   if(sample %in% STD){
#     misty.results$Treatment <- "Standar"
#   }else{ misty.results$Treatment <- "DUTRENEO" }
#   if(sample %in% Nr){
#     misty.results$Response <- "No"
#   }else{misty.results$Response <- "Complete"}
#   
#   MISTY <- rbind(MISTY,misty.results)
# }
# for(t in tme){
#   sample <- gsub(x = t,pattern = "footprint_TME_",replacement = "")
#   misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/All/",t))$importances.aggregated
#   misty.results <- misty.results[-which(is.na(misty.results$Importance)),]  
#   
#   misty.results$Type <- ""
#   misty.results$Response <- ""
#   misty.results$Treatment <- ""
#   misty.results$Cluster <- "TME"
#   if(sample %in% HOT){
#     misty.results$Type <- "HOT"
#   } else{ misty.results$Type <- "COLD" }
#   if(sample %in% STD){
#     misty.results$Treatment <- "Standar"
#   }else{ misty.results$Treatment <- "DUTRENEO" }
#   if(sample %in% Nr){
#     misty.results$Response <- "No"
#   }else{misty.results$Response <- "Complete"}
#   
#   MISTY <- rbind(MISTY,misty.results)
# }

for(w in whole){
  sample <- gsub(x = w,pattern = "footprint_",replacement = "")
  misty.results <- collect_results(paste0("Results_MistyR/Latest_11_April_footprints/Whole/",w))$importances.aggregated
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
#saveRDS(object = MISTY,"Results_MistyR/MISTY.RDS")
MISTY<- readRDS("Results_MistyR/MISTY.RDS")
top <- top[which(top %in% MISTY$Pairs )]

# du_hot <- readRDS("RDS/seurat_v5_SCT/DU24.rds")
# du_cold <- readRDS("RDS/seurat_v5_SCT/DU29.rds")
# du_hot<-AddModuleScore(du_hot,features = list(c("FN1","CCN2")),assay = "SCT",name = "COLD")
# du_cold<-AddModuleScore(du_cold,features = list(c("FN1","CCN2")),assay = "SCT",name = "COLD")
# SpatialFeaturePlot(du_hot,features = "COLD1") + labs(title = "Sample 24 - HOT")
# SpatialFeaturePlot(du_cold,features = "COLD1") + labs(title = "Sample 29 - COLD")

for(p in top){
  df_plot <- MISTY[which(MISTY$Pairs %in% p),]
  pair_HC<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Type)) +
    geom_boxplot()+
    facet_grid(rows = vars(Cluster)) +
    labs(title = p)+
    theme_classic()+
    stat_compare_means(aes(group = Type))
  pair_SD<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Treatment)) +
    geom_boxplot()+
    facet_grid(rows = vars(Cluster)) +
    labs(title = p)+
    theme_classic()+
    stat_compare_means(aes(group = Treatment))
  pair_CN<- ggplot(data = df_plot,aes(x=view,y=Importance, col=Response)) +
    geom_boxplot()+
    facet_grid(rows = vars(Cluster)) +
    labs(title = p)+
    theme_classic()+
    stat_compare_means(aes(group = Response))
  ggsave(filename = paste0("Results_MistyR/Latest_11_April_footprints/Plots/",p,"/",p,"_Cold_vs_Hot.png"),plot = pair_HC,height = 15,width = 10)
  ggsave(filename = paste0("Results_MistyR/Latest_11_April_footprints/Plots/",p,"/",p,"_Standar_vs_Dutreneo.png"),plot = pair_SD,height = 15,width = 10)
  ggsave(filename = paste0("Results_MistyR/Latest_11_April_footprints/Plots/",p,"/",p,"_Complete_vs_No.png"),plot = pair_CN,height = 15,width = 10)
}
vistas <- unique(MISTY$view)
sub_Misty<-MISTY[MISTY$Pairs %in% top,]
sub_Misty$Combination <- paste0(sub_Misty$sample,"_",sub_Misty$Type,"_",sub_Misty$Treatment,"_",sub_Misty$Response)
vis <- c("Intraview","Juxtaview","Juxtaview","Paraview")
# GENERAL PLOTS
i<- 1
for(v in vistas){
  sub <- sub_Misty[sub_Misty$view == v,]
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
  # Extract second elements
  Type <- sapply(labels, "[[", 2)
  # Extract third elements
  Treatment <- sapply(labels, "[[", 3)
  # Extract fourth elements
  Response <- sapply(labels, "[[", 4)
  
  # create a `df` with the samples grouped in the same way you want to show
  anno <- data.frame(Treatment=Treatment,
                     Type =Type,
                     Response=Response)
  
  # Sort columns based on Type
  sorted_columns <- order(anno$Response,anno$Type,anno$Treatment)
  heatmap_data <- heatmap_data[, sorted_columns]
  anno <- anno[sorted_columns, ]
  
  # set rownames so that anno and your data can be matched
  rownames(anno) <- colnames(heatmap_data)
  
  Heat<-pheatmap(heatmap_data, main = vis[i],annotation_col = anno,
           cluster_cols = F)
  ggsave(filename = paste0("Results_MistyR/Latest_11_April_footprints/Plots/",v,"_DefHeatmap.png"),plot = Heat,height = 10,width = 15)
  i <- i+1
}

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
    
    Heat<-pheatmap(heatmap_data,scale = "row", main = v,annotation_col = anno,
                   cluster_cols = F)
    ggsave(filename = paste0("Results_MistyR/Latest_11_April_footprints/Plots/",arm,"/",v,"_Heatmap.png"),plot = Heat,height = 10,width = 15)
  }
}


# misty_top <- readRDS("Results_MistyR/Misty_Top.rds")
# misty_top$Pairs <- as.character(misty_top$Pairs)
# misty_top$Pairs <- as.factor(misty_top$Pairs)
# # Finding unique pairs for each category of "Views"
# pairs_by_views <- misty_top %>%
#   group_by(Views) %>%
#   distinct(Pairs) %>%
#   ungroup()
# 
# pair_counts <- pairs_by_views %>%
#   count(Pairs) %>% filter(n == n_distinct(pairs_by_views$Views))
# # Extracting pairs present in every category of "Views"
# pairs_present_in_all <- pair_counts$Pairs
# 
# # Viewing the pairs present in every category of "Views"
# print(pairs_present_in_all)
# 
# misty_top <- misty_top[which(misty_top$Pairs %in% good_pairs),]
# # HOT VS COLD
# Hot_Cold <- misty_top[which(misty_top$Group=="Hot" | misty_top$Group=="Cold"),]
# # Filter those who are present in all views
# good_pairs <- c()
# for(pa in Hot_Cold$Pairs){
#   if(length(unique(Hot_Cold$Views[Hot_Cold$Pairs==pa]))==6){
#     print(pa)
#     good_pairs <- c(good_pairs, pa)
#   }
# }
# Hot_Cold <- Hot_Cold[which(Hot_Cold$Pairs %in% good_pairs),]
# Top_HC <- Hot_Cold %>%
#   group_by(Group) %>%
#   top_n(10, wt = mean) %>%
#   arrange(desc(mean))
# 
# ggplot(data = Hot_Cold[which(Hot_Cold$Pairs %in% Top_HC$Pairs[1:5]),],aes(x=view,y=mean, col=Group)) +
#   geom_boxplot()+
#   facet_grid(rows = vars(Pairs)) +
#   labs(title = "Hot vs Cold")+
#   theme_classic()
# 
# Hot <- Hot_Cold[Hot_Cold$Group=="Hot",]
# Cold <- Hot_Cold[Hot_Cold$Group=="Cold",]
# # Select Pairs
# Hot <- Hot[!Hot$Pairs %in% Cold$Pairs,]
# Hot <- Hot %>%
#   group_by(Views) %>%
#   top_n(10, wt = mean) %>%
#   ungroup() %>%
#   arrange(Views, desc(mean))
# Hot_unique_pairs <- unique(Hot$Pairs)
# Hot <- misty_top[misty_top$Pairs %in% Hot_unique_pairs ,]
# Hot_unique <- unique(unlist(strsplit(as.character(Hot_unique_pairs),split = "_")))
# Cold <- Cold[!Cold$Pairs %in% Hot$Pairs,]
# Cold <- Cold %>%
#   group_by(Views) %>%
#   top_n(10, wt = mean) %>%
#   ungroup() %>%
#   arrange(Views, desc(mean))
# Cold_unique_pairs <- unique(Cold$Pairs)
# Cold_unique <- unique(unlist(strsplit(as.character(Cold_unique_pairs),split = "_")))
# Cold<- Cold[Cold$Pairs %in% Cold_unique_pairs,]
# ggplot(data = Hot[Hot$Pairs %in% Hot_unique_pairs,],aes(x=Views,y=mean, col=Pairs, group=Pairs)) +
#   geom_point() +
#   geom_line() +
#   labs(title = "Hot")+
#   theme_classic()
# ggplot(data = Cold[Cold$Pairs %in% Cold_unique_pairs[1:10],],aes(x=view,y=mean, col=Pairs, group=Pairs)) +
#   geom_point() +
#   geom_line() +
#   labs(title = "Cold")+
#   theme_classic()

