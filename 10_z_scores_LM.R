library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(scales)
library(stringr)
library(ggplot2)
library(SpatialExperiment)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
samples_id <- list.files(path = "RDS/frontier/Cancer/")
common_Ligands<- readRDS("Results_Correlation/Ligan_list.RDS")
common_receptors<- readRDS("Results_Correlation/Receptor_list.RDS")
lr <- readRDS("RDS/genes_LR.RDS")
lr <- lr[,c(3,4)]
more <- read.csv("RDS/LIANA_genes.csv")
more <- more[,2:3]
colnames(more) <- colnames(lr)
good_lr <- rbind(lr,more)
df <- readRDS("RDS/patients_metatdata.RDS")
#df <- df[c(-5,-19:-20),]
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Pr <- unique(df$sample[which(df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD])
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")
samples <- unique(df$sample)
for(b in 1:3){
  big_df<- readRDS(paste0("Results_Correlation/DF_1%_",arm[b],".rds"))
  
  big_df$correlacion[is.na(big_df$correlacion)] <- 0
  big_df$parejas <- as.factor(big_df$parejas)
  big_df$receptors <- as.factor(big_df$receptors)
  big_df$ligands <- as.factor(big_df$ligands)
  big_df$Sample <- as.factor(big_df$Sample)
  big_df$Reaction <- as.factor(big_df$Reaction)
  
  state <- 1:length(big_df$correlacion)
  big_df$state <- state
  
  big_df$state[which(big_df$correlacion>=0.4)] <- "+ Correlation" 
  big_df$state[which(big_df$correlacion<=(-0.4))] <- "- Correlation" 
  big_df$state[which(big_df$correlacion>(-0.4) & big_df$correlacion<0.4)] <- "No Correlation" 
  
  #normalized_values <- ((big_df$nSpots - min(big_df$nSpots)) / (max(big_df$nSpots) - min(big_df$nSpots)) + 0.5)
  zero_df <-  big_df[which(big_df$correlacion==0),]
  big_df <- big_df[which(big_df$correlacion!=0),]
  #parejas <- c("APP-CD74","HLA-C-LILRB1","HLA-C-LILRB2","CLEC3A-CLEC10A","LYZ-ITGAL","HLA-C-CD8A","HLA-C-CD3D","EFNB1-ERBB2","HLA-C-CD3G","S100A8-ITGB2","COPA-CD74")
  # Mean calculation
  big_df$meandiff_CN <- 0
  lr<-unique(big_df$parejas)
  Exp_Matrices <- list()
  Frontier_Spots <- list()
  for(s in 1:length(arm_samples[[b]])){
    du <- readRDS(paste0("RDS/seurat_v5_SCT/",arm_samples[[b]][s],".rds"))
    Exp_Matrices[[s]] <-du@assays$RNA$data
    frontier <- readRDS(paste0("RDS/frontier/Cancer/",arm_samples[[b]][s],".rds"))
    Frontier_Spots[[s]] <- as.vector(colnames(frontier@assays$RNA$data))
  }
  Frontier_Spots<- as.list(Frontier_Spots)
  Exp_Matrices <- as.list(Exp_Matrices)
  saveRDS(object = Exp_Matrices,file = paste0("RDS/Expression_Matrices_",arm[b],".RDS"))
  saveRDS(object = Frontier_Spots,file = paste0("RDS/Border_ID_Spots_",arm[b],".RDS"))
  Exp_Matrices <- readRDS(paste0("RDS/Expression_Matrices_",arm[b],".RDS"))
  Frontier_Spots <- readRDS(paste0("RDS/Border_ID_Spots_",arm[b],".RDS"))
  LR <- data.frame()
  lr <- gsub(pattern = "HLA-",replacement = "HLA_",lr)
  for(i in 1:length(lr)){
    print(i)
    ligand <- strsplit(x = as.character(lr),split = "-")[[i]][1]
    ligand <-gsub(pattern = "_",replacement = "-",x = ligand)
    receptor <- strsplit(x = as.character(lr),split = "-")[[i]][2]
    receptor <-gsub(pattern = "_",replacement = "-",x = receptor)
    Z_df <- data.frame()
    for(s in 1:length(arm_samples[[b]])){
      z_scores_df <- data.frame()
      du <- Exp_Matrices[[s]]
      # First store the expression matrix for the specific pair
      expr_mat <- du[c(ligand, receptor),colnames(du)]
      # Calculate mean and standard deviation
      mean_ligand <- mean(expr_mat[1,])
      sd_ligand <- sd(expr_mat[1,])
      mean_receptor <- mean(expr_mat[2,])
      sd_receptor <- sd(expr_mat[2,])
      # Calculate z-scores for ligand
      z_scores_ligand <- (expr_mat[1,] - mean_ligand) / sd_ligand
      # Calculate z-scores for receptor
      z_scores_receptor <- (expr_mat[2,] - mean_receptor) / sd_receptor
      # Combine z-scores into a data frame
      z_scores_df <- data.frame(Spot = colnames(expr_mat),
                                Z_Score_ligand = z_scores_ligand,
                                Z_Score_receptor = z_scores_receptor)

      z_scores_df[which(z_scores_df$Spot %in% Frontier_Spots[[s]]),]
      z_scores_df$Sample <- arm_samples[[b]][s]

      if(arm_samples[[1]][s] %in% Nr){
        z_scores_df$Response <- "No"
      }else{
        z_scores_df$Response <- "Complete"
      }
      Z_df <- rbind(Z_df,z_scores_df)
    }
    model <- lm(Z_Score_receptor ~ Z_Score_ligand*Response, data = Z_df)
    sm <- summary(model)
    c_sm<- as.data.frame(coef(sm))

    c_sm$ligand <- ligand
    c_sm$receptor <- receptor
    c_sm$Info <- rownames(c_sm)
    #c_sm$pair <- paste0(c_sm$ligand,"-",c_sm$receptor)
    LR <- rbind(LR,c_sm)
    if(c_sm$`Pr(>|t|)`[4]<0.001 & c_sm$`Pr(>|t|)`[2]<0.001){
      # p <- ggplot(Z_df, aes(x = Z_Score_ligand, y = Z_Score_receptor, color = Response)) +
      #   geom_point() +
      #   geom_smooth(method = "lm", se = FALSE) +
      #   labs(x = paste0("Z-Score ",ligand), y = paste0("Z-Score ",receptor), color = "Treatment Response")
      # ggsave(paste0("Results_Correlation/Models_Zscores/Plots/",ligand,"_",receptor,".png"),plot = p,height = 10,width = 10)
      saveRDS(object = model,paste0("Results_Correlation/Models_Zscores/",arm[b],"/",ligand,"_",receptor,".RDS"))
      saveRDS(object = Z_df,paste0("Results_Correlation/RDS/z_scores/",arm[b],"/",ligand,"_",receptor,".RDS"))
    }
  }
  saveRDS(object = LR,file = paste0("Results_Correlation/RDS/",arm[b],"_LigandReceptor_Models_Pvalue.RDS"))
  LR<- readRDS(paste0("Results_Correlation/RDS/",arm[b],"_LigandReceptor_Models_Pvalue.RDS"))
  # VOLCANO 
  big_df$Reaction[which(big_df$Reaction=="Partial")] <- "Complete"
  # GET ONLY THE PAIRS PRESENT IN BOTH
  LR$pair <- paste0(LR$ligand,"-",LR$receptor)
  pvals_LR <- LR[which(LR$Info=="Z_Score_ligand:ResponseNo"),]
  big_df$pvals_Zscores <- 0
  for(p in lr){
    big_df$pvals_Zscores[which(big_df$parejas==p)] <- pvals_LR$`Pr(>|t|)`[which(pvals_LR$pair==p)]
  }
  big_df$Adj.pval<- p.adjust(p = big_df$pvals_Zscores,method = "BH")
  pairs<- gsub(pattern = ".RDS",replacement = "",x = list.files(path = paste0("Results_Correlation/Models_Zscores/",arm[b],"/")))
  pairs <- gsub(pattern = "_",replacement = "-",x = pairs)
  
  for (i in pairs){
    df <- big_df[which(big_df$parejas==i),]
    # Complete NO DIFF
    Comp_df <- df[which(df$Reaction=="Complete"),]
    No_df <- df[which(df$Reaction=="No"),]
    CN_diff<- mean(Comp_df$correlacion) - mean(No_df$correlacion)
    big_df$meandiff_CN[which(big_df$parejas==i)] <- CN_diff
    
  }
  
  big_df$diffexpressed <- "NO"
  big_df$diffexpressed[-log10(big_df$Adj.pval) > -log10(0.001) & big_df$pvals < 0.001 & abs(big_df$meandiff_CN)>0.2] <- "YES"
  big_df$delabel <- big_df$parejas
  big_df$delabel[which(big_df$diffexpressed == "NO")] <- NA
  
  v_plot<- ggplot(data = big_df, aes(x=meandiff_CN,y=-log10(Adj.pval),col=diffexpressed, label=delabel))+
    geom_point() + 
    theme_bw() +
    geom_text(na.rm = T)+
    geom_vline(xintercept = 0.2, col="green") + 
    #xlim(c(-0.5,0.5)) +
    geom_vline(xintercept = -0.2, col="green") +
    ggtitle(paste0(arm[b]," Arm - Cancer Border"))
  
  v_plot
  #ggsave(filename = paste0("Results_Correlation/Plots/Z_score_1%_Results/1%_Volcano_CompleteNo_",arm[b],".png"),plot = v_plot,width = 15,height = 15)

  high_corr_pairs <- c()
  for(pa in pairs){
    df <- big_df[which(big_df$parejas==pa),]
    if(dim(big_df[which(big_df$parejas==pa),])[1]>3){
      df_No <- df[which(df$Reaction=="No"),]
      df_comp <- df[which(df$Reaction=="Complete"),]
      m_pa_No <- mean(df_No$correlacion)
      m_pa_Comp <- mean(df_comp$correlacion)
      if(m_pa_No>0.4 | m_pa_Comp>0.4){
        high_corr_pairs <- c(high_corr_pairs,pa)
      }
    }
  }
  big_df<-big_df[which(abs(big_df$meandiff_CN)>0.2),]
  CN <- big_df[which(big_df$parejas %in% high_corr_pairs),]
  saveRDS(object = CN,paste0("Results_Correlation/RDS/Complete_vs_No_ZSCORE_",arm[b],".RDS"))
  ggplot(data = CN,aes(x=parejas,y=correlacion, fill=Reaction)) + 
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," - Cancer Border"))
  
  g2 <- ggplot(data = CN,aes(x=parejas,y=correlacion, fill=Reaction)) + 
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," Arm - Cancer Border"))
  IMPORTANT_PAIRS<- unique(g2$data$parejas)
  saveRDS(object = as.vector(IMPORTANT_PAIRS),file = paste0("Results_Correlation/DEFENITIVEPAIRS_",arm[b],".rds"))
  ggsave(filename = paste0("Results_Correlation/Plots/1%_Results/1%_CompleteNo_ZSCORE_",arm[b],".png"),plot = g2,width = 18,height = 10)
}

