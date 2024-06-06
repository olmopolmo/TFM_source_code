library(OmnipathR)
library(Seurat)
library(dplyr)
# Get the Ligand-Receptor Data
lr <- curated_ligand_receptor_interactions()
lr
lr <- lr[which(lr$curation_effort>1),]
lr <- lr[which(lr$n_references>1),]
lr <- lr[which(lr$n_resources>1),]

# Substract the ligands and receptors present in our samples
# Load the list of ligandreceptors from the data
# ligand_list <- unique(lr$source_genesymbol)
# receptor_list <- unique(lr$target_genesymbol)
# # Create two empty vectors to fill with the ligands and receptors present in both the database and the samples
# common_Ligands <- c()
# common_receptors <- c()
samples_id <- list.files(path = "/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/seurat_v5_SCT/")
for (i in samples_id){
  du <- readRDS(paste0("/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/seurat_v5_SCT/",i))
  genes <- as.data.frame(AverageExpression(du,assays = "SCT"))
  genes$genes <- rownames(genes)
  colnames(genes) <- c("mean","genes")
  genes_du <- genes$genes[which(genes$mean>0.5)] 
  lr <- lr[which(lr$source_genesymbol %in% genes_du & lr$target_genesymbol %in% genes_du),]

  #ligand_list <- ligand_list[ligand_list %in% genes_du]
  #common_Ligands <- c(common_l,common_Ligands)
  #receptor_list <- receptor_list[receptor_list %in% genes_du]
  #common_receptors <- c(common_r,common_receptors)
}
#common_Ligands <- unique(ligand_list)
#common_receptors <- unique(receptor_list)
#gene_pairs <- unique(common_Ligands,common_receptors)
#saveRDS(file = paste0("/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/Ligan_list.RDS"),object = common_Ligands)
#saveRDS(file = "/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/Receptor_list.RDS",object = common_receptors)
saveRDS(file = "/Users/ASUS/Desktop/IJC/DUTRENEO/RDS/genes_LR.RDS",object = lr)
