# CHANGE WHEN THE NEW DATA ARRIVES
# HOT -> DU9
# COLD -> DU1
# DUTRENEO -> DU18
# STANDAR -> DU1
# COMPLETE -> DU11
# NO -> DU14
# Generalized
df <- readRDS("RDS/patients_metatdata.RDS")

HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Nr <- unique(df$sample[which(df$response=="NO")])

All <- list(HOT,COLD,STD,DUTRE,Cr,Nr)
All_names <- c("Hot","Cold","Standar","Dutreneo","Complete","No")
comp <- c("CH","CH","DS","DS","CN","CN")
du <- c("DU9","DU1","DU1","DU18","DU11","DU14")
HOT_and_COLD <- c(HOT,COLD)

misty.results <- collect_results( paste0("Results_MistyR/vignette_model_footprints_DU1/"))
misty.results %>%
  plot_improvement_stats("gain.R2") 
misty.results %>%
  plot_improvement_stats("gain.RMSE")
mist <- misty.results$importances.aggregated

Imp <- mist[which(mist$Importance>10),]
Imp <- mist[which(mist$Predictor %in% unique(Imp$Predictor)),]
misty.results$importances.aggregated <- Imp
misty.results %>% plot_interaction_heatmap(view = "intra")
misty.results %>% plot_interaction_heatmap(view = "juxta.1.5")
misty.results %>% plot_interaction_heatmap(view = "juxta.2")
misty.results %>% plot_interaction_heatmap(view = "juxta.2.5")
misty.results %>% plot_interaction_heatmap(view = "juxta.3")
misty.results %>% plot_interaction_heatmap(view = "juxta.3.5")
misty.results %>% plot_interaction_heatmap(view = "juxta.4")
misty.results %>% plot_interaction_heatmap(view = "juxta.4.5")
#mist <- mist[which(mist$Target %in% gp),]

max(mist$Importance)
mist$Importance[which(is.na(mist$Importance))]=0
misty.results %>% plot_view_contributions()
parejas <- c()

# Guarda las parejas en un vector
for(i in 1:length(mist$Predictor)){
  p <- paste0(mist$Predictor[i],"_",mist$Target[i])
  parejas <- c(parejas, p)
}
mist$Pairs <- parejas
mist$Pairs <- as.factor(mist$Pairs)
parejas<-unique(mist$Pairs)
# mist$Importance[which(is.na(mist$Importance))]=0
# means_p <- mist %>% group_by(Predictor) %>%
#   summarize(mean_importance = mean(Importance))
# 
# gp<-means_p$Predictor[which(means_p$mean_importance>0.1)]
# mist <- mist[which(mist$Predictor %in% gp),]
# 
# means_t <- mist %>% group_by(Target) %>%
#   summarize(mean_importance = mean(Importance))
# tp <- means_t$Target[which(means_t$mean_importance>0.1)]
# mist <- mist[which(mist$Target %in% tp),]
# misty.results$importances.aggregated <- mist
# misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.1)


for(j in 1:6){
  misty <- data.frame()
  # Crea una df con todos los datos de las muestras COLD
  for (i in All[[j]]) {
    mist <- readRDS(file = paste0("Results_MistyR/RDS/CAFs/Misty_",i,"_",comp[j],".RDS"))
    mist$sample <- i
    misty <- rbind(misty,mist)
  }
  parejas <- c()
  # Guarda las parejas en un vector
  for(i in 1:length(misty$Predictor)){
    p <- paste0(misty$Predictor[i],"_",misty$Target[i])
    parejas <- c(parejas, p)
  }
  misty$Pairs <- parejas
  misty$Pairs <- as.factor(misty$Pairs)
  parejas<-unique(misty$Pairs)
  
  # CALCULATE THE MEANS THAT THE PAIRS SHOULD HAVE
  misty$mean <- 0
  for(i in parejas){
    #Intra
    intra <- which(misty$Pairs==i & misty$view=="intra")
    mean_to_change <- mean(misty$Importance[intra])
    misty$mean[intra] <- mean_to_change
    #Juxta 1.5
    juxta <- which(misty$Pairs==i & misty$view=="juxta.1.5")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 2
    juxta <- which(misty$Pairs==i & misty$view=="juxta.1.5")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 2.5
    juxta <- which(misty$Pairs==i & misty$view=="juxta.2")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 3 
    juxta <- which(misty$Pairs==i & misty$view=="juxta.2.5")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 3.5
    juxta <- which(misty$Pairs==i & misty$view=="juxta.3.5")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 4
    juxta <- which(misty$Pairs==i & misty$view=="juxta.4")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
    #Juxta 4.5
    juxta <- which(misty$Pairs==i & misty$view=="juxta.4.5")
    mean_to_change <- mean(misty$Importance[juxta])
    misty$mean[juxta] <- mean_to_change
  }
  # Divide la df en base al analysis de mistyR
  intra <- misty[which(misty$view=="intra"),]
  juxta <- misty[which(misty$view=="juxta.1.5"),]
  intra$Importance[is.na(intra$Importance)] = 0
  juxta$Importance[is.na(juxta$Importance)] = 0
  misty_data <- readRDS(paste0("Results_MistyR/RDS/",k,"/Misty_",du[j],"_",comp[j],".RDS"))
  # sustituye en el objeto con todos los genes el valor medio de la importancia
  for(i in 1:length(misty_data$view)){
    if(misty_data$view[i]=="intra"){
      misty_data$Importance[i] == intra$Importance[which(misty_data$Predictor[i]==intra$Predictor & misty_data$Target[i]==intra$Target)[1]]
    }
    else{
      misty_data$Importance[i] == juxta$Importance[which(misty_data$Predictor[i]==juxta$Predictor & misty_data$Target[i]==juxta$Target)[1]]
    }
  }
  misty.results <- collect_results( paste0("MistyR/vignette_model_footprints_",du[j],"_",comp[j]))
  misty.results$importances.aggregated$Importance <- misty_data$Importance
  png(paste0("Results_MistyR/RDS/",All_names[j],"_intra.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.1.5", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.2.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.2", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.2.5.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.2.5", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.3.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.3", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.3.5.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.3.5", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.4.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.4", cutoff = 0.1)
  dev.off()
  png(paste0("Results_MistyR/RDS/",All_names[j],"_juxta.4.5.png"),width = 1000,height = 1000)
  misty.results %>% plot_interaction_heatmap(view = "juxta.4.5", cutoff = 0.1)
  dev.off()
  
}  