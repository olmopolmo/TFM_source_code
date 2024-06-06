library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
#library(OmnipathR)
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
#df <- df[c(5,19:20),]
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Pr <- unique(df$sample[which(df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
arm_samples <-list(HOT[HOT %in% DUTRE],HOT[HOT %in% STD],COLD[COLD %in% STD])
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")
big_df <- data.frame()

for(b in 1:2){
  gc()
  samples_id <- arm_samples[[b]]
  brazo <- readRDS(paste0("Results_Correlation/BigDF_",arm[b],".rds"))
  for(i in samples_id){
    du <-  readRDS(paste0("RDS/frontier/Cancer/",i,".rds"))
    DefaultAssay(du) <- "SCT"
    Idents(du) <- du$compartment_NewBuff
    #du <- du[,which(du$compartment==paste0("Cancer") | du$compartment==paste0("TME"))]
    rownames(du@tools$Staffli@meta.data) <- paste0(str_split(rownames(du@tools$Staffli@meta.data), "-", simplify = T)[,1], "-",  str_split(colnames(du), "-", simplify = T)[1,2])
    sub_du <- du
    # Run the RegionNeighbours() function from STUtility to select the adjacent spots
    sub_du <- SetIdent(sub_du, value = "nbs_Cancer")
    lr <- good_lr[which(good_lr$source_genesymbol %in% rownames(sub_du)),]
    lr <- lr[which(lr$target_genesymbol %in% rownames(sub_du)),]
    w_matrix <- as.data.frame(sub_du@assays$RNA@data)
    parejas <- c()
    correlacion <- c()
    receptors <- c()
    ligands <- c()
    pvals <- c()
    porcentage_expression_L <- c()
    porcentage_expression_R <- c()
    lr<- unique(lr)
    for (j in 1:length(lr$source_genesymbol)){
      ligand_gene <- lr$source_genesymbol[j]
      receptor_gene <- lr$target_genesymbol[j]
      ligand <- w_matrix[which(rownames(w_matrix)==ligand_gene),]
      receptor <- w_matrix[which(rownames(w_matrix)==receptor_gene),]
      percentage_ligand <- length(which(ligand>0))/length(ligand)
      percentage_receptor <- length(which(receptor>0))/length(receptor)

      if(percentage_receptor<0.01 | percentage_ligand<0.01){ #  & (mean_receptor<1 | mean_ligand<1)
        correlacion <- c(correlacion,0)
        pvals <- c(pvals,1)
      }else{
        correlation <- cor.test(unlist(ligand),unlist(receptor))
        correlacion <- c(correlacion,correlation$estimate)
        pvals <- c(pvals,correlation$p.value)
      }
      ligands <- c(ligands,ligand_gene)
      receptors <- c(receptors, receptor_gene)
      parejas <- c(parejas,paste0(ligand_gene,"-",receptor_gene))
      porcentage_expression_L <- c(porcentage_expression_L,percentage_ligand)
      porcentage_expression_R <- c(porcentage_expression_R,percentage_receptor)
    }
    
    cor_pairs<- data.frame(correlacion,parejas,ligands,receptors, pvals,porcentage_expression_R,porcentage_expression_L)
    #write.table(cor_pairs,file = paste0("Results_Correlation/",i,"_Correlation_5percent_scores_",arm[b],".txt"))
    df <- cor_pairs
    df$Sample <- paste0(i)
    if(gsub(x = i,pattern = ".rds",replacement = "") %in% Nr){
      df$Reaction <- "No"
    }else{
      df$Reaction <- "Complete"
    }
    brazo<- rbind(brazo,df)
  }
  saveRDS(object = brazo,file = paste0("Results_Correlation/BigDF_",arm[b],".rds"))
}

for(b in 1:3){
  big_df<- readRDS(paste0("Results_Correlation/DF_1%_",arm[b],".rds"))
  
  big_df$nSpots <- 0
  for(i in unique(big_df$Sample)){
    du <- readRDS(paste0("RDS/frontier/Cancer/",i,".rds"))
    big_df$nSpots[which(big_df$Sample==i)]<- ncol(du)
  }
  
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
  
  big_df$CorrSpots <- big_df$nSpots*big_df$correlacion
  #normalized_values <- ((big_df$nSpots - min(big_df$nSpots)) / (max(big_df$nSpots) - min(big_df$nSpots)) + 0.5)
  zero_df <-  big_df[which(big_df$correlacion==0),]
  big_df <- big_df[which(big_df$correlacion!=0),]
  #parejas <- c("APP-CD74","HLA-C-LILRB1","HLA-C-LILRB2","CLEC3A-CLEC10A","LYZ-ITGAL","HLA-C-CD8A","HLA-C-CD3D","EFNB1-ERBB2","HLA-C-CD3G","S100A8-ITGB2","COPA-CD74")
  # Mean calculation
  big_df$meandiff_CN <- 0
  pairs<-unique(big_df$parejas)
  for (i in pairs){
    df <- big_df[which(big_df$parejas==i),]
    # Complete NO DIFF
    Comp_df <- df[which(df$Reaction=="Complete" | df$Reaction=="Partial"),]
    No_df <- df[which(df$Reaction=="No"),]
    CN_diff<- mean(Comp_df$CorrSpots) - mean(No_df$CorrSpots)
    big_df$meandiff_CN[which(big_df$parejas==i)] <- CN_diff
    
  }
  big_df <- big_df[-which(is.na(big_df$meandiff_CN)),]
  ####################
  #### GOOD PLOTS ####
  ####################
  big_df$Reaction[which(big_df$Reaction=="Partial")] <- "Complete"
  pairs <- unique(big_df$parejas)
  
  significant_CN <- c()
  big_df$pvals <- 1
  big_df$log10.pval <- 0
  for(x in pairs){
    p <- big_df[which(big_df$parejas==x),]
    if(sum(p$correlacion) > 0 & length(p$Reaction)!=length(p$Reaction[which(p$Reaction=="No")]) & 
       length(p$Reaction)!=length(p$Reaction[which(p$Reaction=="Complete")]) & 
       length(p$correlacion[which(p$Reaction=="No")])>0 &
       length(p$correlacion[which(p$Reaction=="Complete")])>0 ){
      tt_CN <- wilcox.test(p$correlacion[which(p$Reaction=="Complete")],p$correlacion[which(p$Reaction=="No")])
      big_df$pvals[which(big_df$parejas==x)] <- tt_CN$p.value
      big_df$log10.pval[which(big_df$parejas==x)] <- -log10(tt_CN$p.value)
      if (tt_CN$p.value<0.05){
        significant_CN <- c(significant_CN,x)
      }
    }
  }
  SIG_PAIRS<- significant_CN
  
  quantile(abs(big_df$meandiff_CN),0.99)

  big_df$diffexpressed <- "NO"
  big_df$diffexpressed[big_df$log10.pval > log10(0.01) & big_df$pvals < 0.01 & abs(big_df$meandiff_CN)>quantile(abs(big_df$meandiff_CN),0.9)] <- "YES"
  big_df$delabel <- big_df$parejas
  big_df$delabel[which(big_df$diffexpressed == "NO")] <- NA
  
  v_plot<- ggplot(data = big_df, aes(x=meandiff_CN,y=log10.pval,col=diffexpressed, label=delabel))+
    geom_point() + 
    theme_bw() +
    geom_text(na.rm = T)+
    geom_vline(xintercept = -quantile(abs(big_df$meandiff_CN),0.9), col="green") + 
    #xlim(c(-0.5,0.5)) +
    geom_vline(xintercept = quantile(abs(big_df$meandiff_CN),0.9), col="green") +
    geom_hline(yintercept=-log10(0.01), col="red") +
    ggtitle(paste0(arm[b]," Arm - Cancer Border"))
  v_plot
  ggsave(filename = paste0("Results_Correlation/Plots/1%_Results/1%_Volcano_CompleteNo_",arm[b],".png"),plot = v_plot,width = 15,height = 15)
  
  #Plots 
  high_corr_pairs <- c()
  for(pa in significant_CN){
    df <- big_df[which(big_df$parejas==pa),]
    df_No <- df[which(df$Reaction=="No"),]
    df_comp <- df[which(df$Reaction=="Complete"),]
    m_pa_No <- mean(df_No$correlacion)
    m_pa_Comp <- mean(df_comp$correlacion)
    if(m_pa_No>0.4 | m_pa_Comp>0.4){
      high_corr_pairs <- c(high_corr_pairs,pa)
    }
  }
  CN <- big_df[which(big_df$parejas %in% high_corr_pairs),]
  saveRDS(object = CN,paste0("Results_Correlation/RDS/Complete_vs_No_df_",arm[b],".RDS"))
  ggplot(data = CN,aes(x=parejas,y=correlacion, fill=Reaction)) + 
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," - Cancer Border"))
  
  # Complete vs No
  g2 <- ggplot(data = CN,aes(x=parejas,y=correlacion, fill=Reaction)) + 
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," Arm - Cancer Border"))
  
  ggsave(filename = paste0("Results_Correlation/Plots/1%_Results/1%_CompleteNo_",arm[b],".png"),plot = g2,width = 18,height = 10)
  
  
  big_df$FDR <-p.adjust(big_df$pvals, method = "fdr")
  CN$FDR <-  p.adjust(CN$pvals, method = "fdr")
  
} 

