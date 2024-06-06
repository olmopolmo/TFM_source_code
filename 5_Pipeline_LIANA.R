library(Seurat)
library(tidyverse)
library(magrittr)
library(liana)
library(tidyr)
setwd(dir = "/Users/ASUS/Desktop/IJC/DUTRENEO/")
# CREATE THE LIANA OBJECTS
n <- list.files(path = "RDS/seurat_v5_SCT/")
for (i in n){
  du <- readRDS(paste0("RDS/seurat_v5_SCT/",i))
  Idents(du) <-du$compartment
  du
  # Run liana
  liana_test <- liana_wrap(sce = du,resource = ,assay = "SCT")
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% dplyr::glimpse()

  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate()
  saveRDS(file = paste0("Results_LIANA/RDS/LIANA_2_",i),object = liana_test)
}


pm <- readRDS("RDS/patients_metatdata.RDS")

HOT <- unique(pm$sample[which(pm$type=="HOT")])
COLD <- unique(pm$sample[which(pm$type=="COLD")])
STD<- unique(pm$sample[which(pm$treatment=="STANDARD")])
DUTRE<- unique(pm$sample[which(pm$treatment!="STANDARD")])
Cr <- unique(pm$sample[which(pm$response=="COMPLETE")])
Pr <- unique(pm$sample[which(pm$response=="PARTIAL")])
Nr <- unique(pm$sample[which(pm$response=="NO")])

COLD<- gsub(pattern = "DU",x = COLD,replacement = "")
HOT<- gsub(pattern = "DU",x = HOT,replacement = "")
STD<- gsub(pattern = "DU",x = STD,replacement = "")
DUTRE<- gsub(pattern = "DU",x = DUTRE,replacement = "")
Cr<- gsub(pattern = "DU",x = Cr,replacement = "")
Pr<- gsub(pattern = "DU",x = Pr,replacement = "")
Nr<- gsub(pattern = "DU",x = Nr,replacement = "")

classes <- list(COLD,HOT,STD,DUTRE,Cr,Pr,Nr) 
names <- c("Cold","Hot","Standar","DUTRENEO","Complete","Partial","No")
n <- list.files(path = "RDS/seurat_v5_SCT/")
genes_du <- c()
for(i in n){  
  du <- readRDS(paste0("/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/seurat_v5_SCT/",i))
  genes <- as.data.frame(AverageExpression(du,assays = "SCT"))
  genes$genes <- rownames(genes)
  colnames(genes) <- c("mean","genes")
  genes_du <- c(genes_du,genes$genes[which(genes$mean>5)])
}
genes_du <- unique(genes_du)
for (c in 1:7){
  big_df <- data.frame()
  for(i in classes[[c]]){
    du <- readRDS(paste0("Results_LIANA/RDS/LIANA_DU",i,".rds"))
    plot1 <- du %>%
      liana_dotplot(source_groups = c("Cancer"), # CHANGE WHEN THE NEW DATA ARRIVES FOR THE NEW CLUSTERS
                    target_groups = c("TME")) 
    # In this state the plot is not readable we have to filter
    df <- plot1$data
    df$sample <- paste0("DU",i)
    
    # Change the symbol for clarity
    df$interaction <- sub(x = df$interaction,pattern = "-",replacement = "_")
    df$interaction <- sub(x = df$interaction,pattern = " _> ",replacement = "-")
    
    #df<- df[which(df$interaction %in% unique(CH$parejas)),]
    plot1$data <- df
    big_df <- rbind(big_df,df)
    
  }
  # Remove genes with low expression in the samples
  # big_df$interaction
  # big_df <- separate(big_df, interaction, into = c("G1", "G2"), sep = "-",remove = F)
  # big_df <- big_df[which(big_df$G1 %in% genes_du & big_df$G2 %in% genes_du),]
  # For loop in order to calculate the mean per target
  gp <- unique(big_df$interaction)
  targets <- unique(big_df$target)
  for (j in gp) {
    small_df <- big_df[which(big_df$interaction==j),]
    for(k in targets){
      big_df$magnitude[which(big_df$interaction==j & big_df$target == k)] <- mean(small_df$magnitude[which(small_df$target == k)])
    }
  }
  big_df$interaction<- as.factor(big_df$interaction)
  g_pairs<- unique(names(which(table(big_df$interaction)>=length(classes[[c]])-4))) # Modify this number based on what classification are you comparing 4*(number of samples - 1)
  df<- big_df[which(big_df$interaction %in% g_pairs),]
  thresh <- sd(big_df$magnitude)*1.5+mean(big_df$magnitude)

  th<- mean(big_df$magnitude)
  high_mean <- function(points, mean_threshold = thresh) { # Should be modified to take the top 5%
    mean_value <- mean(points)
    return(mean_value > mean_threshold)
  }
  # FILTER BASED ON DIFFERENCES, TAKE ONLY TE MOST DIFFERENT PAIRS 

  # Filter sets based on standard deviation so I get those who vary the most.
  mean_df <- big_df %>%
    group_by(interaction) %>%
    filter(high_mean(magnitude))
  filtered_df <-mean_df %>%
    group_by(interaction) 
  plot1 <- du %>%
    liana_dotplot(source_groups = c("Cancer"), # CHANGE WHEN THE NEW DATA ARRIVES FOR THE NEW CLUSTERS
                  target_groups =  c("TME"))
  #plot1$data <- mean_df
  plot1$data <- filtered_df
  ggsave(paste0("Results_LIANA/Plots/plot_",names[c],".png"),plot = plot1, height = 30, width = 15)
}

big_df <- readRDS("Results_LIANA/Borders/Arms/HotDUTRENEO/filtered_df_HotDUTRENEO_Cancer_Comp-No.RDS")
# Filtering: 
#   Take the log2FC that are higher than 0.5 and smaller than -0.5
#   Take the specificity higher than 5%
df <- big_df[which(big_df$specificity>=0.01),]
df <- df[which(abs(df$logfc.logfc_comb)>=0.5),]

pvp <- readRDS("RDS/frontier/Cancer_DU15.rds")
pvp2 <- readRDS("RDS/frontier/Cancer_DU19.rds")

SpatialFeaturePlot(pvp,features = c("TFF2","CXCR4"))
SpatialFeaturePlot(pvp2,features = c("TFF2","CXCR4"))

# SpatialFeaturePlot(pvp,features = c("APP","CD74"))
# SpatialFeaturePlot(pvp2,features = c("APP","CD74"))
