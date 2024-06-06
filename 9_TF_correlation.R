library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
library(scales)
library(stringr)
library(ggplot2)
library(reshape2)
library(pheatmap)

setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
TF<- readRDS("TF_acts.rds")
samples_id <- list.files(path = "RDS/frontier/Cancer/")

big_df<- readRDS(paste0("Results_Correlation/BigDF_ColdStandar.rds"))

df <- readRDS("RDS/patients_metatdata.RDS")
df <- df[-c(5,19, 20),]
HOT <- unique(df$sample[which(df$type=="HOT")])
COLD <- unique(df$sample[which(df$type=="COLD")])
STD<- unique(df$sample[which(df$treatment=="STANDARD")])
DUTRE<- unique(df$sample[which(df$treatment!="STANDARD")])
Cr <- unique(df$sample[which(df$response=="COMPLETE")])
Pr <- unique(df$sample[which(df$response=="PARTIAL")])
Nr <- unique(df$sample[which(df$response=="NO")])
arm_samples <-list(HOT[HOT %in% DUTRE],COLD[COLD %in% STD],COLD[COLD %in% STD])
arm <- c("HotDUTRENEO","HotStandar","ColdStandar")

gene_pairs <- list(readRDS("Results_Correlation/DEFENITIVEPAIRS_HotDUTRENEO.rds"),
                   readRDS("Results_Correlation/DEFENITIVEPAIRS_HotStandar.rds"),
                   readRDS("Results_Correlation/DEFENITIVEPAIRS_ColdStandar.rds"))
HD_HC <- as.data.frame(table(unlist(list(readRDS("Results_Correlation/DEFENITIVEPAIRS_HotDUTRENEO.rds"),
                                         readRDS("Results_Correlation/DEFENITIVEPAIRS_HotStandar.rds")))))
HD_CC <- as.data.frame(table(unlist(list(readRDS("Results_Correlation/DEFENITIVEPAIRS_HotDUTRENEO.rds"),
                                         readRDS("Results_Correlation/DEFENITIVEPAIRS_ColdStandar.rds")))))
HC_CC <- as.data.frame(table(unlist(list(readRDS("Results_Correlation/DEFENITIVEPAIRS_HotStandar.rds"),
                                         readRDS("Results_Correlation/DEFENITIVEPAIRS_ColdStandar.rds")))))
for(b in 1:3){
  gc()
  samples_id <- arm_samples[[b]]
  brazo <- data.frame()
  pairs <- gene_pairs[[b]]# readRDS(paste0("Results_Correlation/StatisticallySignificantPairs_",arm[b],".RDS"))
  for(i in samples_id[8:10]){
    du <-  readRDS(paste0("RDS/frontier/Cancer/",i,".rds"))
    DefaultAssay(du) <- "SCT"
    Idents(du) <- du$nbs_Cancer
    rownames(du@tools$Staffli@meta.data) <- paste0(str_split(rownames(du@tools$Staffli@meta.data), "-", simplify = T)[,1], "-",  str_split(colnames(du), "-", simplify = T)[1,2])
    ids <- colnames(du)
    # CHANGE: FROM THE TF DF TAKE THE SPOTS FOR THE CORREPONDING SAMPLE AND RUN THE COR_TEST FOR BOTH, THEN STORE THAT IN A MATRIX
    # WHERE THE ROWS CAN BE THE GENE PAIRS THE COLS THE TF AND THE MEAN VALUE IS STORED IN EACH CELL.
    w_matrix <- as.data.frame(du@assays$SCT$scale.data)
    spots_id <- colnames(du)
    sub_TF <- TF[which(TF$sample==i),]
    sub_TF <- sub_TF[which(sub_TF$spotid %in% spots_id),]
    transcription_factors <- colnames(sub_TF[-c(1,length(sub_TF))])
    cor.ligand <- c()
    cor.receptor <- c()
    mean.corr <- c()
    diff.corr <- c()
    receptors <- c()
    ligands <- c()
    response <- c()
    sample.id <- c()
    parejas <- c()
    transcription.factor.name.for.df <- c()
    for (j in transcription_factors){
      for(k in pairs){
        ligand_gene <- strsplit(x = k,split = "-")[[1]][1]
        receptor_gene <- strsplit(x = k,split = "-")[[1]][2]
        ligand <- w_matrix[which(rownames(w_matrix)==ligand_gene),]
        receptor <- w_matrix[which(rownames(w_matrix)==receptor_gene),]
        spots <- sub_TF$spotid
        tf_mat <- sub_TF[,`j`]
        if(dim(receptor)[1]==0 | dim(ligand)[1]==0){
          cor.val.receptor <- 0
          cor.val.ligand <- 0
        }else{
          cor.val.ligand <- cor.test(unlist(ligand),tf_mat)
          cor.val.receptor <- cor.test(unlist(receptor),tf_mat)
          if(cor.val.receptor$p.value<0.05 & cor.val.ligand$p.value<0.05){
            cor.val.ligand <- cor.val.ligand$estimate
            cor.val.receptor <- cor.val.receptor$estimate
          }else{
            cor.val.receptor <- 0
            cor.val.ligand <- 0
          }
        }
        ligands <- c(ligands,ligand_gene)
        receptors <- c(receptors,receptor_gene)
        cor.ligand <- c(cor.ligand,cor.val.ligand)
        cor.receptor <- c(cor.receptor,cor.val.receptor)
        mean.corr <- c(mean.corr, (cor.val.ligand+cor.val.receptor)/2)
        diff.corr <- c(diff.corr, abs(cor.val.ligand-cor.val.receptor))
        parejas <- c(parejas,k)
        sample.id <- c(sample.id,i)
        transcription.factor.name.for.df <- c(transcription.factor.name.for.df,j)
        if(i %in% Nr){
          response <- c(response,"No")
        }else{
          response <- c(response,"Complete")
        }
      }
    }
    
    TF_cor_dataframe<- data.frame(cor.ligand,cor.receptor,mean.corr,diff.corr,receptors,ligands,parejas,response,sample.id,transcription.factor.name.for.df)
    saveRDS(TF_cor_dataframe,file = paste0("Results_Correlation/TF_correlation_scores/1%_",i,"_TF_scores_",arm[b],".rds"))
    brazo<- rbind(brazo,TF_cor_dataframe)
    
  }
  saveRDS(object = brazo,file = paste0("Results_Correlation/TF_correlation_scores/1%_BigDF_",arm[b],".rds"))
}
  
for(b in 1:3){
  big_df<- readRDS(paste0("Results_Correlation/TF_correlation_scores/1%_BigDF_",arm[b],".rds"))
  #big_df<- big_df[big_df$mean.corr!=0,]
  big_df$parejas <- as.factor(big_df$parejas)
  big_df$receptors <- as.factor(big_df$receptors)
  big_df$ligands <- as.factor(big_df$ligands)
  big_df$sample.id <- as.factor(big_df$sample.id)
  big_df$response <- as.factor(big_df$response)
  
  big_df$tf.pair <- paste0(big_df$transcription.factor.name.for.df,"_",big_df$parejas)
  # Mean calculation
  big_df$CorrMean <- 0
  triunvirates<-unique(big_df$tf.pair)

  #######################
  # Complete RESPONDERS #
  #######################
  Comp_df <- big_df[big_df$response=="Complete",]
  Comp_df$CorrMean <- 0
  Comp_df$LigandMean <- 0
  Comp_df$ReceptorMean <- 0
  for (i in triunvirates){
    df <- Comp_df[which(Comp_df$tf.pair==i),]
    Corr_mean <- mean(df$mean.corr)
    Comp_df$CorrMean[which(Comp_df$tf.pair==i)] <- Corr_mean
    Ligand_Mean <- mean(df$cor.ligand)
    Comp_df$LigandMean[which(Comp_df$tf.pair==i)] <- Ligand_Mean
    Receptor_Mean <- mean(df$cor.receptor)
    Comp_df$ReceptorMean[which(Comp_df$tf.pair==i)] <- Receptor_Mean
    
  }
  state <- 1:length(Comp_df$cor.ligand)
  Comp_df$state <- state
  ggplot(Comp_df,aes(x=CorrMean)) + geom_histogram()
  Comp_df$state[which(Comp_df$CorrMean>(-0.3) & Comp_df$CorrMean<0.3)] <- "No Correlation" 
  Comp_df$state[Comp_df$CorrMean>=0.3] <- "+ Correlation" 
  Comp_df$state[Comp_df$CorrMean<=(-0.3)] <- "- Correlation" 
  
  ####################
  #### GOOD PLOTS ####
  ####################
  tf_Complete<-unique(Comp_df$transcription.factor.name.for.df[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
  parejas_complete<-unique(Comp_df$parejas[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
  triunvirates_Complete<-unique(Comp_df$tf.pair[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
  Comp_df <- Comp_df[which(Comp_df$parejas %in% parejas_complete),]
  Comp_df <- Comp_df[which(Comp_df$transcription.factor.name.for.df %in% tf_Complete),]
  Comp_df <- Comp_df[which(Comp_df$sample.id!="DU33"),]
  Comp_df <- Comp_df[which(Comp_df$sample.id!="DU3"),]
  Comp_df <- Comp_df[which(Comp_df$sample.id!="DU9"),]
  Comp_df <- Comp_df[which(Comp_df$sample.id!="DU19"),]
  
  #saveRDS(object = Comp_df,paste0("Results_Correlation/TF_correlation_scores/TF_Complete_vs_No_df_",arm[b],".RDS"))
  ggplot(data = Comp_df,aes(x=transcription.factor.name.for.df,y=mean.corr)) +
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," Arm - Cancer Border")) + facet_grid(cols = vars(parejas))
  # 
  # Complete vs No
  if(dim(Comp_df)[1]>0){
    g2 <- ggplot(data = Comp_df,aes(x=transcription.factor.name.for.df,y=mean.corr)) + 
      geom_boxplot(position = "dodge") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle(paste0(arm[b]," Arm - Cancer Border - Complete Response")) + facet_grid(cols = vars(parejas))
    
    ggsave(filename = paste0("Results_Correlation/TF_correlation_scores/Plots/Complete_",arm[b],".png"),plot = g2,width = 18,height = 10)
  }
  #################
  # NO RESPONDERS #
  #################
  No_df <- big_df[big_df$response=="No",]
  No_df$CorrMean <- 0
  No_df$LigandMean <- 0
  No_df$ReceptorMean <- 0
  for (i in triunvirates){
    df <- No_df[which(No_df$tf.pair==i),]
    Corr_mean <- mean(df$mean.corr)
    No_df$CorrMean[which(No_df$tf.pair==i)] <- Corr_mean
    Ligand_Mean <- mean(df$cor.ligand)
    No_df$LigandMean[which(No_df$tf.pair==i)] <- Ligand_Mean
    Receptor_Mean <- mean(df$cor.receptor)
    No_df$ReceptorMean[which(No_df$tf.pair==i)] <- Receptor_Mean
    
  }
  state <- 1:length(No_df$cor.ligand)
  No_df$state <- state
  ggplot(No_df,aes(x=CorrMean)) + geom_histogram()
  No_df$state[which(No_df$CorrMean>(-0.3) & No_df$CorrMean<0.3)] <- "No Correlation" 
  No_df$state[No_df$CorrMean>=0.3] <- "+ Correlation" 
  No_df$state[No_df$CorrMean<=(-0.3)] <- "- Correlation" 
  
  ####################
  #### GOOD PLOTS ####
  ####################
  tf_No<-unique(No_df$transcription.factor.name.for.df[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
  triunvirates_No<-unique(No_df$tf.pair[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
  parejas_No<-unique(No_df$parejas[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
  No_df <- No_df[which(No_df$parejas %in% parejas_No),]
  No_df <- No_df[which(No_df$transcription.factor.name.for.df %in% tf_No),]
  #saveRDS(object = No_df,paste0("Results_Correlation/TF_correlation_scores/TF_Complete_vs_No_df_",arm[b],".RDS"))
  # ggplot(data = No_df,aes(x=tf.pair,y=mean.corr)) +
  #   geom_boxplot(position = "dodge") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   scale_fill_manual(values = c("green","blue")) +
  #   ggtitle(paste0(arm[b]," Arm - Cancer Border"))
  if(dim(No_df)[1]>0){
    # Complete vs No
    g2 <- ggplot(data = No_df,aes(x=transcription.factor.name.for.df,y=mean.corr)) + 
      geom_boxplot(position = "dodge") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle(paste0(arm[b]," Arm - Cancer Border - No Response")) + facet_grid(cols = vars(parejas))
    
    ggsave(filename = paste0("Results_Correlation/TF_correlation_scores/Plots/No_",arm[b],".png"),plot = g2,width = 18,height = 10)
  }
  
  
  Both_df <- big_df
  Both_df$CorrMean <- 0
  for (i in triunvirates){
    df <- Both_df[which(Both_df$tf.pair==i),]
    Corr_mean <- mean(df$mean.corr)
    Both_df$CorrMean[which(Both_df$tf.pair==i)] <- Corr_mean
    
  }
  ####################
  #### GOOD PLOTS ####
  ####################
  triunvirates<-c(triunvirates_No,triunvirates_Complete)
  Both_df <- Both_df[which(Both_df$tf.pair %in% triunvirates),]
  #saveRDS(object = Comp_df,paste0("Results_Correlation/TF_correlation_scores/TF_Complete_vs_No_df_",arm[b],".RDS"))
  g2  <- ggplot(data = Both_df,aes(x=transcription.factor.name.for.df,y=mean.corr, fill = response)) +
    geom_boxplot(position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("green","blue")) +
    ggtitle(paste0(arm[b]," Arm - Cancer Border")) + facet_grid(cols = vars(parejas))
  
  ggsave(filename = paste0("Results_Correlation/TF_correlation_scores/Plots/Complet_No_",arm[b],".png"),plot = g2,width = 18,height = 10)
  saveRDS(triunvirates,file = paste0("Results_Correlation/TF_correlation_scores/TF_Pair_List",arm[b],".RDS"))
}  

parejas <- unique(c(parejas_complete,parejas_No))
triunvirates <- unique(c(triunvirates_No,triunvirates_Complete))
Comp_df$response <- "Cr-Pr"
No_df$response <- "Nr"
Both_df <- rbind(Comp_df,No_df)
# Variable with all the transcription factors of interest
transcription_factors <- unique(tf_Complete,tf_No)
Both_df <- Both_df[which(Both_df$transcription.factor.name.for.df %in% transcription_factors),]
Both_df <- Both_df[which(Both_df$tf.pair %in% triunvirates),]

Both_df<- Both_df %>% arrange(parejas) #Ordena por parejas para que consiga un vector a llenar con NA por pareja
ligand_cor_values <- Both_df$LigandMean
receptor_cor_values <- Both_df$ReceptorMean
tf.pair <- Both_df$tf.pair 
pair <- Both_df$parejas
response <- Both_df$response
df <- data.frame(tf.pair, receptor_cor_values, ligand_cor_values,pair, response)
df<- df %>% arrange(pair) 
df <- unique(df)

library(tidyverse)
# Extract TF names
df$tf <- gsub("_.*", "", df$tf.pair)

# Pivot the dataframe
receptor_df <- df %>%
  select(tf, pair, receptor_cor_values,response) %>%
  spread(pair, receptor_cor_values)
receptor_df$fun <- "receptor"
ligand_df <- df %>%
  select(tf, pair, ligand_cor_values,response) %>%
  spread(pair, ligand_cor_values)
ligand_df$fun <- "ligand"
pivot_df <- rbind(ligand_df,receptor_df)
pivot_df <- pivot_df %>% arrange(response) %>% arrange(tf)

# paste0 con la tf fun y response para tener ids unicos
pivot_df$tf.interaction <- paste0(pivot_df$tf,".",pivot_df$fun,".",pivot_df$response)
pivot_df$fun <- factor(pivot_df$fun,levels = c("ligand","receptor"))
pivot_df$response <- factor(pivot_df$response,levels = c("Cr-Pr","Nr"))
molded_df <- data.frame()
for(i in unique(pivot_df$tf)){
  sub_df <- pivot_df[which(pivot_df$tf==i),c(3:(length(colnames(pivot_df))-2))] 
  new_row <- unlist(sub_df)
  sufixes <- rep(c("_Cr-Pr_ligand",
                   "_Cr-Pr_receptor",
                   "_Nr_ligand",
                   "_Nr_receptor"),4)
  names(new_row)<- paste0(gsub("\\d$", "", names(new_row)),sufixes)
  new_row<- t(as.data.frame(new_row))
  rownames(new_row) <- 1
  molded_df<- rbind(molded_df,new_row)
}
rownames(molded_df) <- unique(pivot_df$tf)
molded_matrix <- as.matrix(molded_df)

pairs_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[1])
response_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[2])
fun_ids <-sapply(strsplit(colnames(molded_df), "_"), function(x) x[3])
data.frame(pairs_ids,response_ids,fun_ids)
pheatmap(mat = molded_matrix)



annotation <- data.frame(fun_ids,response_ids,pairs_ids)
rownames(annotation) <- colnames(molded_matrix) # check out the row names of annotation
pheatmap(molded_matrix, annotation_col  = annotation,cluster_rows = F,cluster_cols = F,legend_breaks = T)

