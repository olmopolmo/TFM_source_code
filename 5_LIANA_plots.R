library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(liana)
library(tidyr)
setwd("../Desktop/IJC/DUTRENEO/")
arms <- c("HotDUTRENEO","HotStandar","ColdStandar")
for(i in arms){
  du <- readRDS(paste0("Results_LIANA/RDS/Arms/liana_test_",i,"_Cancer_Comp-No.RDS"))
  n1 <- "Complete Response"
  n2 <- "No Response"
  du$source<- gsub(pattern ="Cancer_No",replacement = "Cancer Border",x = du$source)
  du$source<- gsub(pattern ="Cancer_Comp",replacement = "Cancer Border",x = du$source)
  du$target<- gsub(pattern ="TME_Comp",replacement = "Complete Response",x = du$target)
  du$target<- gsub(pattern ="TME_No",replacement = "No Response",x = du$target)
  
  plot1 <- du %>%
    liana_dotplot(source_groups = "Cancer Border", 
                  target_groups = c(n1,n2))
  
  plot1
  filtered_df <- plot1$data
  
  specificity_90 <-quantile(filtered_df$specificity,0.95) 
  filtered_df <- filtered_df[which(filtered_df$specificity>=specificity_90),]
  magnitude_90 <-quantile(filtered_df$magnitude,0.95)
  filtered_df <- filtered_df[which(filtered_df$magnitude>=magnitude_90),]
  
  plot1$data <- filtered_df
  plot1
  ggsave(paste0("RESULTADOS_PARA_TFM/LIANA/Arms/",i,".png"),plot = plot1, height = 10, width = 12)
  
}
