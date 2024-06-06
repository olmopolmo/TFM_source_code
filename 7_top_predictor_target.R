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
setwd("../Desktop/IJC/DUTRENEO/")
# The objective of this script is to get the pairs with most significance in each distance from misty 
misty.results <- collect_results( paste0("Results_MistyR/vignette_model_footprints_DU1/"))
du$compartment
# Method 1 based from RMSE and R2
m1 <- misty.results$improvements.stats
m1$mean[which(is.na(m1$mean))] <- 0
top_m1 <- c()
measure <- unique(m1$measure)
for(i in measure){
  tmp <- m1[which(m1$measure==i),]
  top<- tmp %>% arrange(desc(mean)) %>% head(10)
  top_m1 <- c(top$target,top_m1)
}
top_m1<- unique(top_m1)
top_m1


# Method 2 by the contribution
m2 <- misty.results$contributions.stats
m2$fraction[which(is.na(m2$fraction))] <- 0
top_m2 <- c()
views <- unique(m2$view)
for(i in views){
  tmp <- m2[which(m2$view==i),]
  top<- tmp %>% arrange(desc(fraction)) %>% head(10)
  top_m2 <- c(paste0(top$target),top_m2)
}
top_m2<- unique(top_m2)
top_m2

# Method 3 by the interaction
m3 <- misty.results$importances.aggregated
m3$Importance[which(is.na(m3$Importance))] <- 0

top_10 <- c()
views <- unique(m3$view)
for(i in views){
  tmp <- m3[which(m3$view==i),]
  top<- tmp %>% arrange(desc(Importance)) %>% head(10)
  top_10 <- c(paste0(top$Predictor,"-",top$Target),top_10)
}
top_10<- unique(top_10)
