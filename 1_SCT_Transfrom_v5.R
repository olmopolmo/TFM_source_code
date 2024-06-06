# Transform to SCT the samples
samples_id <- list.files(path = "RDS/seurat/")
for (i in samples_id){
  gc()
  du <-  readRDS(paste0("RDS/seurat/",i))
  du<- SCTransform(du,vst.flavor="v2",return.only.var.genes = F)
  DefaultAssay(du) <- "SCT"
  saveRDS(file = paste0("RDS/seurat_v5_SCT/",i),object = du)
}
