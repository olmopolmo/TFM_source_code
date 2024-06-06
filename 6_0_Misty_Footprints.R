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
saveRDS()
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
for_misty <- readRDS("RDS/genes_LR.RDS")
for_misty <- c(for_misty$source_genesymbol,for_misty$target_genesymbol)
# WARINING: MISTY DOES NOT LIKE GENES THAT HAVE A "-" for the moment i will remove them
for_misty <- for_misty[-which(grepl(pattern = "-",x = for_misty))]
compartments <- c("Cancer","Buffer","TME") # "TME" fnished
for(c in compartments){ # Du 50
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
    
    slide.markers <- names(which(coverage[which(names(coverage)%in%for_misty)]>0.05))
    slide.markers <- slide.markers[which(slide.markers %in% rownames(norm.data))]
    misty.views <- create_initial_view(t(norm.data)[, slide.markers] %>% as_tibble())
    misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 2)
    misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 3)
    misty.views <- misty.views %>% add_paraview(geometry,l=10)
    run_misty(misty.views, paste0("Results_MistyR/New_footprints/",c,"/footprint_",c,"_",gsub(x = i,pattern = ".rds",replacement = "")),num.threads = 1)
    #misty.results <- collect_results( paste0("Results_MistyR/footprint_Cancer_",gsub(x = i,pattern = ".rds",replacement = "")))
  }
}
##########################################################################################
# CALCULATE FOOTPRINTS FOR NEIGHBOURS: CANCER-BUFFER, TME-BUFFER
##########################################################################################
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
All <- list(HOT,COLD,STD,DUTRE,Cr,Nr)
All_names <- c("Hot","Cold","Standar","Dutreneo","Complete","No")
comp <- c("CH","CH","DS","DS","CN","CN")
views <- c("intra","juxta.2","juxta.3","Para")
HOT_and_COLD <- c(HOT,COLD)
misty_top <- data.frame() 
Misty_big <- data.frame()
for(c in compartments){
  for(j in 1:6){
    misty <- data.frame()
    parejas <- c()
    for (i in All[[j]]) {
      mist <- collect_results( paste0("Results_MistyR/New_footprints/",c,"/footprint_",c,"_",i,"/"))
      mist <- mist$importances.aggregated
      mist$sample <- i
      misty <- rbind(misty,mist)
    }
    parejas<- paste0(misty$Predictor,"_",misty$Target)
    misty$Pairs <- parejas
    #misty<- misty[-grep("^COL", misty$Pairs),]
    parejas <- unique(names(table(parejas)[which(table(parejas)==max(table(parejas)))]))
    misty <- misty[which(misty$Pairs %in% unique(parejas)),]
    
    misty$Pairs <- as.factor(misty$Pairs)
    parejas<- unique(misty$Pairs[which(misty$Importance>10)])
    # CALCULATE THE MEANS THAT THE PAIRS SHOULD HAVE
    misty$mean <- 0
    for(pair in parejas){
      #Intra
      intra <- which(misty$Pairs==pair & misty$view=="intra")
      mean_to_change <- mean(misty$Importance[intra])
      misty$mean[intra] <- mean_to_change
      #Juxta 2
      juxta <- which(misty$Pairs==pair & misty$view=="juxta.2")
      mean_to_change <- mean(misty$Importance[juxta])
      misty$mean[juxta] <- mean_to_change
      #Juxta 3
      juxta <- which(misty$Pairs==pair & misty$view=="juxta.3")
      mean_to_change <- mean(misty$Importance[juxta])
      misty$mean[juxta] <- mean_to_change
      #Paraview
      juxta <- which(misty$Pairs==pair & misty$view=="para.10")
      mean_to_change <- mean(misty$Importance[juxta])
      misty$mean[juxta] <- mean_to_change
  
    }
    
    misty$Importance[is.na(misty$Importance)] = 0
    misty <- misty[which(misty$mean>0),]
    misty <- misty[,-c(4,5,6)]
    misty <- distinct(misty)
    misty$Group <- All_names[j]
    mis<- misty %>% group_by(view) %>% arrange(desc(mean))
    mis <- misty[misty$Pairs %in% unique(mis$Pairs),]
    # genes <- as.vector(unique(mis$Pairs))
    # genes <- as.data.frame(genes)
    # genes$Sample <- All_names[j]
    # gens <- rbind(gens,genes)
    misty_top <- rbind(misty_top,mis)
    #Misty_big <- rbind(Misty_big,misty)
  }
  saveRDS(object = misty_top,file = paste0("Results_MistyR/Misty_Top_",c,".rds"))
  misty_top <- data.frame() 
}  
# COMPLETE VS NO
Complete <- misty_top[misty_top$Group=="Complete",]
No <- misty_top[misty_top$Group=="No",]

Complete <- Complete[!Complete$Pairs %in% No$Pairs,]
Complete <- Complete %>% group_by(view) %>% arrange(desc(mean))
Complete_unique_pairs <- unique(Complete$Pairs)[1:10]
Complete_unique <- unique(unlist(strsplit(as.character(Complete_unique_pairs),split = "_")))
No <- No[!No$Pairs %in% Complete$Pairs,]
No <- No %>% group_by(view) %>% arrange(desc(mean))
No_unique_pairs <- unique(No$Pairs)[1:10]
No_unique <- unique(unlist(strsplit(as.character(No_unique_pairs),split = "_")))

ggplot(data = Complete[Complete$Pairs %in% Complete_unique_pairs,],aes(x=view,y=mean, col=Pairs, group=Pairs)) +
  geom_point() +
  geom_line() +
  labs(title = "Complete")+
  theme_classic()
ggplot(data = No[No$Pairs %in% No_unique_pairs,],aes(x=view,y=mean, col=Pairs, group=Pairs)) +
  geom_point() +
  geom_line() +
  labs(title = "No")+
  theme_classic()
du_Complete <- readRDS("RDS/seurat_v5_SCT/DU24.rds")
du_No <- readRDS("RDS/seurat_v5_SCT/DU29.rds")
du_Complete<-AddModuleScore(du_Complete,features = list(No_unique),assay = "SCT",name = "No")
du_No<-AddModuleScore(du_No,features = list(No_unique),assay = "SCT",name = "No")
SpatialFeaturePlot(du_Complete,features = "No1") + labs(title = "Sample 24 - Complete")
SpatialFeaturePlot(du_No,features = "No1") + labs(title = "Sample 29 - No")



# DUTRENEO VS STANDAR
Dutreneo <- misty_top[misty_top$Group=="Dutreneo",]
Standar <- misty_top[misty_top$Group=="Standar",]

Dutreneo <- Dutreneo[!Dutreneo$Pairs %in% Standar$Pairs,]
Dutreneo <- Dutreneo %>% group_by(view) %>% arrange(desc(mean))
Dutreneo_unique_pairs <- unique(Dutreneo$Pairs)[1:10]
Dutreneo_unique <- unique(unlist(strsplit(as.character(Dutreneo_unique_pairs),split = "_")))
Standar <- Standar[!Standar$Pairs %in% Dutreneo$Pairs,]
Standar <- Standar %>% group_by(view) %>% arrange(desc(mean))
Standar_unique_pairs <- unique(Standar$Pairs)[1:10]
Standar_unique <- unique(unlist(strsplit(as.character(Standar_unique_pairs),split = "_")))

ggplot(data = Dutreneo[Dutreneo$Pairs %in% Dutreneo_unique_pairs,],aes(x=view,y=mean, col=Pairs, group=Pairs)) +
  geom_point() +
  geom_line() +
  labs(title = "Dutreneo")+
  theme_classic()
ggplot(data = Standar[Standar$Pairs %in% Standar_unique_pairs,],aes(x=view,y=mean, col=Pairs, group=Pairs)) +
  geom_point() +
  geom_line() +
  labs(title = "Standar")+
  theme_classic()
############
# HEATMAPS #
############

Cold <- misty_top[misty_top$Group=="Cold",]
Cold <- Cold[!Cold$Pairs %in% Hot$Pairs,]
Cold <- Cold %>% group_by(view) %>% arrange(desc(mean))
Cold_unique_pairs <- unique(Cold$Pairs)[1:5]
Cold_unique <- unique(unlist(strsplit(as.character(Cold_unique_pairs),split = "_")))
Cold<- Cold[which(Cold$Predictor %in% Cold_unique & Cold$Target %in% Cold_unique),]

misty.results <- collect_results( paste0("Results_MistyR/footprints/vignette_model_footprints_DU15"))
misty.results <-  misty_results
misty.results$importances.aggregated <- Cold

colnames(misty.results$importances.aggregated)<- c("view","Predictor","Target","Pairs","Importance","Group")
png(paste0("MistyR/Plots/",All_names[j],"_intra.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.1.5", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.2.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.2", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.2.5.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.2.5", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.3.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.3", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.3.5.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.3.5", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.4.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.4", cutoff = 0.1)
dev.off()
png(paste0("MistyR/Plots/",All_names[j],"_juxta.4.5.png"),width = 1000,height = 1000)
misty.results %>% plot_interaction_heatmap(view = "juxta.4.5", cutoff = 0.1)
dev.off()




