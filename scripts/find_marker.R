library(SeuratObject)
library(Seurat)

# annotations <- readRDS("/datacommons/ochoalab/ssns_gwas/seurat/Gbadegesin_annotations_1.rds")
# pbmc_norm <- NormalizeData(annotations, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc <- FindVariableFeatures(pbmc_norm, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# 
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# DimPlot(pbmc, reduction = "pca")
# DimHeatmap(pbmc, dims = 7:12, cells = 500, balanced = TRUE)
# DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")
# #ElbowPlot(pbmc)
# pbmc <- FindNeighbors(pbmc, dims = 1:15)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# 
# pbmc <- RunUMAP(pbmc, dims = 1:15)
# #DimPlot(pbmc, reduction = "umap")
# #DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
# Markers <- FindAllMarkers(pbmc)
# save(Markers, file = "/datacommons/ochoalab/ssns_gwas/seurat/umap_markers_onset_remission.RData")
# 

# new onset vs remission
load("/datacommons/ochoalab/ssns_gwas/seurat/harmonized_seurat_3_4_5_6.RData")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:30)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.4)

Markers <- FindAllMarkers(harmonized_seurat)
save(Markers, file = "/datacommons/ochoalab/ssns_gwas/seurat/umap_markers_onset_remission_res04.RData")

# harmonized_seurat$new_identity <- ifelse(harmonized_seurat$orig.ident %in% c("sample3", "sample5"), "Group1",
#                                          ifelse(harmonized_seurat$orig.ident %in% c("sample4", "sample6"), "Group2", "Other"))
# 
# # Set the new identity class
# Idents(harmonized_seurat) <- "new_identity"
# table(Idents(harmonized_seurat))

# paired_samples_markers <- FindMarkers(harmonized_seurat, ident.1 = "Group1", ident.2 = "Group2", min.pct = 0, logfc.threshold = 0)
# save(paired_samples_markers, file = "/datacommons/ochoalab/ssns_gwas/seurat/paired_samples_markers_3_4_5_6.RData")

# new onset vs control
# load("/datacommons/ochoalab/ssns_gwas/seurat/harmonized_seurat_remission_control.RData")
# 
# harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:30)
# harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
# harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
# 
# harmonized_seurat$new_identity <- ifelse(harmonized_seurat$orig.ident %in% c("sample4", "sample6"), "Group1",
#                                          ifelse(harmonized_seurat$orig.ident %in% c("sample1", "sample8"), "Group2", "Other"))
# 
# # Set the new identity class
# Idents(harmonized_seurat) <- "new_identity"
# table(Idents(harmonized_seurat))
# DefaultAssay(harmonized_seurat) <- "RNA"
# paired_samples_markers <- FindMarkers(harmonized_seurat, ident.1 = "Group1", ident.2 = "Group2", min.pct = 0, logfc.threshold = 0)
# save(paired_samples_markers, file = "/datacommons/ochoalab/ssns_gwas/seurat/markers_remission_control_RNA.RData")
