#Creating patient metadata
sample <- c()
type <- c()
treatment <- c()
response <- c()
samples_id <- list.files(path = "RDS/seurat_v5_SCT/")
for(i in samples_id){
  gc()
  du <-  readRDS(paste0("RDS/seurat_v5_SCT/",i))
  sample <- c(sample,unique(du$sample))
  type <- c(type,unique(du$type))
  treatment <- c(treatment,unique(du$treatment))
  response <- c(response,unique(du$response))
}
patients <- data.frame(sample,type,treatment,response)
colnames(patients) <- c("sample","type","treatment","response")
saveRDS(object = patients,file = "RDS/patients_metatdata.RDS")

