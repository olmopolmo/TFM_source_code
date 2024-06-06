library(Seurat)
library(dplyr)
library(liana)
samples_id <- list.files(path = "/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/seurat_v5_SCT/")
gene_pairs <- data.frame()
pairs <- c()
for(i in samples_id){
  pve <- readRDS(paste0("/Users/ASUS/Desktop/IJC/DUTRENEO/Results_LIANA/RDS/LIANA_",i))
  pairs <- c(pairs,paste0(pve$ligand.complex,"-",pve$receptor.complex))

}

pairs <- pairs[!duplicated(pairs)]
splt<- strsplit(x = pairs,split = "-")
gene1 <- sapply(splt, `[`, 1)
gene2 <- sapply(splt, `[`, 2)
df <- data.frame(pairs,gene1,gene2)
write.csv(x = df,file = "/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/LIANA_genes.csv",row.names = T)
