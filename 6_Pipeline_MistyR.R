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
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD])
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")
tis <- c("TME","Cancer")

ligands <- readRDS("RDS/DEFINITIVE_LIGANDS.rds")
receptor <- readRDS("RDS/DEFINITIVE_RECEPTORS.rds")



# RUNNING MISTY WITH THE WHOLE SET OF LIANA + OMNIPATH GENE PAIRS
for_misty <- c(ligands,receptor)
for(i in n[-1]){
  du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
  for_misty <- c(ligands,receptor)
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

  for_misty <- names(which(coverage[which(names(coverage)%in%for_misty)]>0.1))
  for_misty <- for_misty[-which(grepl(pattern = "-",x = for_misty))]
  # WARINING: MISTY DOES NOT LIKE GENES THAT HAVE A "-" for the moment i will remove them
  slide.markers <- for_misty
  slide.markers <- slide.markers[which(slide.markers %in% rownames(norm.data))]
  misty.views <- create_initial_view(t(norm.data)[, slide.markers] %>% as_tibble()) %>%
    add_paraview(geometry, l=5)
  misty.views <- misty.views %>% add_views(create_view("paraview.misty.5", misty.views[["paraview.5"]]$data, "para.misty.5"))
  misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 1)
  misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 2)
  misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 3)
  misty.views <- misty.views %>% add_juxtaview(geometry,neighbor.thr = 4)
  run_misty(misty.views, paste0("Results_MistyR/Final_footprints/1%_footprints_",gsub(pattern = ".rds",replacement = "",x = i)),num.threads = 2)

  misty.results <- collect_results( paste0("Results_MistyR/Final_footprints/1%_footprints_",gsub(pattern = ".rds",replacement = "",x = i)))

  saveRDS(object = misty.results$importances.aggregated,file = paste0("Results_MistyR/RDS/Misty_DU",gsub(pattern = ".rds",replacement = "",x = i),".RDS"))
  
}


# Generalized
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Nr <- unique(df$sample[which(df$response=="NO")])
du <- c(HOT[1],COLD[1],STD[1],DUTRE[1],Cr[1],Nr[1])

arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD])
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")

for(j in 1:3){
  misty <- data.frame()
  parejas <- c()
  # Crea una df con todos los datos de las muestras COLD
  for (i in arm_samples[[j]]) {
    mist <- readRDS(file = paste0("Results_MistyR/RDS/Misty_DU",gsub(pattern = ".rds",replacement = "",x = i),".RDS"))
    mist$sample <- i
    misty <- rbind(misty,mist)
  }
  
  parejas<- paste0(misty$Predictor,"_",misty$Target)
  misty$Pairs <- parejas
  parejas <- unique(names(table(parejas)[which(table(parejas)==42)]))
  misty <- misty[which(misty$Pairs %in% unique(parejas)),]
  
  misty$Pairs <- as.factor(misty$Pairs)

  #misty <- misty[which(misty$Pairs==parejas),]
  # CALCULATE THE MEANS THAT THE PAIRS SHOULD HAVE
  misty$mean <- 0
  for(i in parejas){
    #Intra
    intra <- which(misty$Pairs==i & misty$view=="intra")
    mean_to_change <- mean(misty$Importance[intra])
    misty$mean[intra] <- mean_to_change

        #Juxta 2
    juxta <- which(misty$Pairs==i & misty$view=="juxta.2")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change

    #Juxta 3
    juxta <- which(misty$Pairs==i & misty$view=="juxta.3")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change

    #Juxta 4
    juxta <- which(misty$Pairs==i & misty$view=="juxta.4")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change

  }

  # Divide la df en base al analysis de mistyR
  intra <- misty[which(misty$view=="intra"),]
  juxta <- misty[which(misty$view=="juxta.1.5"),]
  intra$Importance[is.na(intra$Importance)] = 0
  juxta$Importance[is.na(juxta$Importance)] = 0
  misty_data <- readRDS(paste0("Results_MistyR/RDS/Misty_frontier",du[j],".RDS"))
  # sustituye en el objeto con todos los genes el valor medio de la importancia
  for(i in 1:length(misty_data$view)){
    if(misty_data$view[i]=="intra"){
      misty_data$Importance[i] == intra$Importance[which(misty_data$Predictor[i]==intra$Predictor & misty_data$Target[i]==intra$Target)[1]]
    }
    else{
      misty_data$Importance[i] == juxta$Importance[which(misty_data$Predictor[i]==juxta$Predictor & misty_data$Target[i]==juxta$Target)[1]]
    }
  }
  misty.results <- collect_results( paste0("Results_MistyR/footprint_frontiers_",du[j]))
  misty.results$importances.aggregated$Importance <- misty_data$Importance
  Imp_Agg <- misty.results$importances.aggregated
  max_target <- Imp_Agg %>%
    group_by(Target) %>%
    summarise(max_importance = max(Importance, na.rm = TRUE))
  
  max_pred <- Imp_Agg %>%
    group_by(Predictor) %>%
    summarise(max_importance = max(Importance, na.rm = TRUE))
  # Filter gene pairs
  thresh<- sd(max_target$max_importance,na.rm = T)*2+mean(max_target$max_importance,na.rm = T)
  Targets <- max_target$Target[max_target$max_importance>thresh]
  thresh<- sd(max_pred$max_importance,na.rm = T)*2+mean(max_pred$max_importance,na.rm = T)
  Predictors <- max_pred$Predictor[max_pred$max_importance>thresh]
  # Intra
  Non_zero <- misty.results$importances.aggregated[which(misty.results$importances.aggregated$Importance>4.67),]
  intra_non <- Non_zero[which(Non_zero$view=="intra"),]
  intra_thresh <- summary(intra_non$Importance)[5]

  #Juxta 2
  Non_zero <- misty.results$importances.aggregated[which(misty.results$importances.aggregated$Importance>3.35),]
  jux2_non <- Non_zero[which(Non_zero$view=="juxta.2"),]
  jux2_thresh <- summary(jux2_non$Importance)[5]

  #Juxta 3
  Non_zero <- misty.results$importances.aggregated[which(misty.results$importances.aggregated$Importance>3),]
  jux3_non <- Non_zero[which(Non_zero$view=="juxta.3"),]
  jux3_thresh <- summary(jux3_non$Importance)[5]

  #Juxta 4
  Non_zero <- misty.results$importances.aggregated[which(misty.results$importances.aggregated$Importance>3.66),]
  jux4_non <- Non_zero[which(Non_zero$view=="juxta.4"),]
  jux4_thresh <- summary(jux4_non$Importance)[5]

  genes <- unique(Predictors,Targets)
  misty.results$importances.aggregated <- Imp_Agg[which(Imp_Agg$Predictor %in% genes & Imp_Agg$Target %in% genes),]
  png(paste0("Results_MistyR/Plots/frontier_",arm[j],"_intra.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = intra_thresh)
  dev.off()
  png(paste0("Results_MistyR/Plots/frontier_",arm[j],"_juxta.2.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.2", cutoff = jux2_thresh)
  dev.off()
  png(paste0("Results_MistyR/Plots/",arm[j],"_juxta.3.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.3", cutoff = jux3_thresh)
  dev.off()
  png(paste0("Results_MistyR/Plots/",arm[j],"_juxta.4.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.4", cutoff = jux4_thresh)
  dev.off()
}  

















