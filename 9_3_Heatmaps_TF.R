library(Seurat)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(spatstat)
#library(OmnipathR)
library(scales)
library(stringr)
library(ggplot2)
library(reshape2)
library(devtools)
library(ComplexHeatmap)


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

b=1
big_df<- readRDS(paste0("Results_Correlation/TF_correlation_scores/1%_BigDF_",arm[b],".rds"))


#big_df<- big_df[big_df$mean.corr!=0,]
big_df$parejas <- as.factor(big_df$parejas)
big_df$receptors <- as.factor(big_df$receptors)
big_df$ligands <- as.factor(big_df$ligands)
big_df$Sample <- as.factor(big_df$Sample)
big_df$response <- as.factor(big_df$Reaction)

big_df$tf.pair <- paste0(big_df$transcription.factor.name.for.df,"_",big_df$parejas)
# big_df <- big_df[which(big_df$sample.id!="DU33"),]
# big_df <- big_df[which(big_df$sample.id!="DU3"),]
# big_df <- big_df[which(big_df$sample.id!="DU9"),]
# big_df <- big_df[which(big_df$sample.id!="DU19"),]

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
Comp_df$state[which(Comp_df$CorrMean>(-0.3) & Comp_df$CorrMean<0.3)] <- "No Correlation" 
Comp_df$state[Comp_df$CorrMean>=0.3] <- "+ Correlation" 
Comp_df$state[Comp_df$CorrMean<=(-0.3)] <- "- Correlation" 

tf_Complete<-unique(Comp_df$transcription.factor.name.for.df[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
parejas_complete<-unique(Comp_df$parejas[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
triunvirates_Complete<-unique(Comp_df$tf.pair[which(Comp_df$state=="+ Correlation" | Comp_df$state=="- Correlation")])
# Comp_df <- Comp_df[which(Comp_df$parejas %in% parejas_complete),]
# Comp_df <- Comp_df[which(Comp_df$transcription.factor.name.for.df %in% tf_Complete),]

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
No_df$state[which(No_df$CorrMean>(-0.3) & No_df$CorrMean<0.3)] <- "No Correlation" 
No_df$state[No_df$CorrMean>=0.3] <- "+ Correlation" 
No_df$state[No_df$CorrMean<=(-0.3)] <- "- Correlation" 

tf_No<-unique(No_df$transcription.factor.name.for.df[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
triunvirates_No<-unique(No_df$tf.pair[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
parejas_No<-unique(No_df$parejas[which(No_df$state=="+ Correlation" | No_df$state=="- Correlation")])
# No_df <- No_df[which(No_df$parejas %in% parejas_No),]
# No_df <- No_df[which(No_df$transcription.factor.name.for.df %in% tf_No),]



parejas <- unique(c(parejas_complete,parejas_No))
transcription_factors <- unique(tf_Complete,tf_No)
Comp_df$response <- "Cr-Pr"
No_df$response <- "Nr"


#################
# NO RESPONDERS #
#################
# Both_df <- rbind(Comp_df,No_df)
Both_df <- No_df
# Variable with all the transcription factors of interest

Both_df <- Both_df[which(Both_df$transcription.factor.name.for.df %in% transcription_factors),]
Both_df <- Both_df[which(Both_df$parejas %in% parejas),]

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
  sufixes <- rep(c("_Nr_ligand",
                   "_Nr_receptor"),length(colnames(pivot_df))-4)
  names(new_row)<- paste0(gsub("\\d$", "", names(new_row)),sufixes)
  new_row<- t(as.data.frame(new_row))
  rownames(new_row) <- 1
  molded_df<- rbind(molded_df,new_row)
}
rownames(molded_df) <- unique(pivot_df$tf)
molded_matrix_No <- as.matrix(molded_df)

pairs_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[1])
response_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[2])
fun_ids <-sapply(strsplit(colnames(molded_df), "_"), function(x) x[3])

annotation <- data.frame(fun_ids,pairs_ids)
rownames(annotation) <- colnames(molded_matrix_No) # check out the row names of annotation
pheatmap(molded_matrix_No, annotation_col  = annotation,cluster_rows = F,cluster_cols = F,scale = "none",show_colnames = F)

col_pairs<- colnames(pivot_df[3:(length(colnames(pivot_df))-2)])
split_elements <- strsplit(col_pairs, "-")
col_pairs <- unlist(lapply(split_elements, function(x) rep(paste(x, collapse = "-"), each = 2)))

Heatmap(molded_matrix_No, heatmap_legend_param = list(title = "Correlation TF~LR"),
        cluster_columns = F,
        cluster_rows = F,
        column_split = col_pairs,
        show_column_names = F,
        row_title  = "No Responders | TF ~ Ligand-Receptor Interactions")

#######################
# Complete RESPONDERS #
#######################
# Both_df <- rbind(Comp_df,No_df)
Both_df <- Comp_df
# Variable with all the transcription factors of interest
Both_df <- Both_df[which(Both_df$transcription.factor.name.for.df %in% transcription_factors),]
Both_df <- Both_df[which(Both_df$parejas %in% parejas),]

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
                   "_Cr-Pr_receptor"),length(colnames(pivot_df))-4)
  names(new_row)<- paste0(gsub("\\d$", "", names(new_row)),sufixes)
  new_row<- t(as.data.frame(new_row))
  rownames(new_row) <- 1
  molded_df<- rbind(molded_df,new_row)
}
rownames(molded_df) <- unique(pivot_df$tf)
molded_matrix_Comp <- as.matrix(molded_df)

pairs_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[1])
response_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[2])
fun_ids <-sapply(strsplit(colnames(molded_df), "_"), function(x) x[3])

annotation <- data.frame(fun_ids,pairs_ids)
rownames(annotation) <- colnames(molded_matrix_Comp) # check out the row names of annotation
pheatmap(molded_matrix_Comp, annotation_col  = annotation,cluster_rows = F,cluster_cols = F,scale = "none",show_colnames = F)

col_pairs<- colnames(pivot_df[3:(length(colnames(pivot_df))-2)])
split_elements <- strsplit(col_pairs, "-")
col_pairs <- unlist(lapply(split_elements, function(x) rep(paste(x, collapse = "-"), each = 2)))

Heatmap(molded_matrix_Comp, heatmap_legend_param = list(title = "Correlation TF~LR"),
        cluster_columns = F,
        cluster_rows = F,
        column_split = col_pairs,
        show_column_names = F,
        row_title  = "Complete Responders | TF ~ Ligand-Receptor Interactions")


# BOTH

molded_matrix<- cbind(molded_matrix_No,molded_matrix_Comp)
pairs_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[1])
response_ids <- sapply(strsplit(colnames(molded_df), "_"), function(x) x[2])
fun_ids <-sapply(strsplit(colnames(molded_df), "_"), function(x) x[3])

col_pairs<- colnames(pivot_df[3:(length(colnames(pivot_df))-2)])
split_elements <- strsplit(col_pairs, "-")
col_pairs <- unlist(lapply(split_elements, function(x) rep(paste(x, collapse = "-"), each = 2)))
col_pairs <- c(col_pairs,col_pairs)

ha = HeatmapAnnotation(response = c(rep("Nr",length(pairs_ids)),rep("Cr-Pr",length(pairs_ids))), 
                       annotation_legend_param = list(
                         response = list(
                           title = "Response",
                           at = c("Nr", "Cr-Pr"),
                           labels = c("No Response", "Complete-Partial Response")
                         )))
ht<- Heatmap(molded_matrix, heatmap_legend_param = list(title = "Correlation TF~LR"),
        top_annotation = ha,
        cluster_columns = F,
        cluster_rows = F,
        column_split = col_pairs,
        show_column_names = F,
        row_title  = "Complete Responders | TF ~ Ligand-Receptor Interactions")

draw(ht, merge_legend = TRUE)
