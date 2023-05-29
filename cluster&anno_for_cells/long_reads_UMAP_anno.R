library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(harmony)

sce <- list()
sce[[1]] <- readRDS("donor0.rds")
sce[[2]] <- readRDS("donor1.rds")
sce[[3]] <- readRDS("donor2.rds")
sce[[4]] <- readRDS("donor3.rds")

sample <- rep("donor0",nrow(sce[[1]]@meta.data))
sce[[1]] <- AddMetaData(
  object = sce[[1]],
  metadata = sample,
  col.name = 'sample')

sample <- rep("donor1",nrow(sce[[2]]@meta.data))
sce[[2]] <- AddMetaData(
  object = sce[[2]],
  metadata = sample,
  col.name = 'sample')

sample <- rep("donor2",nrow(sce[[3]]@meta.data))
sce[[3]] <- AddMetaData(
  object = sce[[3]],
  metadata = sample,
  col.name = 'sample')

sample <- rep("donor3",nrow(sce[[4]]@meta.data))
sce[[4]] <- AddMetaData(
  object = sce[[4]],
  metadata = sample,
  col.name = 'sample')

merged_seurat_filtered <- merge(sce[[1]], 
                                y = c(sce[2:length(sce)]),
                                project = "sce")


# QC

seurat.list <- SplitObject(object = merged_seurat_filtered, split.by = "sample")


gc()
sceList_filtered <- lapply(X = seurat.list, FUN = function(x) {
  
  table(x@meta.data$orig.ident) 
  
  # position <- substr(x@meta.data[["orig.ident"]],5,6)
  # x <- AddMetaData(
  #   object = x,
  #   metadata = position,
  #   col.name = 'position')
  # 
  # age <- substr(x@meta.data[["orig.ident"]],7,8)
  # x <- AddMetaData(
  #   object = x,
  #   metadata = age,
  #   col.name = 'age')
  
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  metadata <- x@meta.data
  metadata$cells <- rownames(metadata)
  
  rownames(x)[grepl('^mt-',rownames(x),ignore.case = T)]
  rownames(x)[grepl('^RP[SL]',rownames(x),ignore.case = T)]
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  
  rb.genes <- rownames(x)[grep("^RP[SL]",rownames(x),ignore.case = T)]
  C<-GetAssayData(object = x, slot = "counts")
  percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
  x<-AddMetaData(x, percent.ribo, col.name = "percent.ribo")
  
  hist.genes <- rownames(x)[grep("^HIST",rownames(x),ignore.case = T)]
  C<-GetAssayData(object = x, slot = "counts")
  percent.hist <- Matrix::colSums(C[hist.genes,])/Matrix::colSums(C)*100
  x<-AddMetaData(x, percent.hist, col.name = "percent.hist")
  
  counts <- GetAssayData(x, assay = 'RNA')
  counts <- counts[-(which(rownames(counts) %in% c('Hbb','Hba1','Hba2','Hbg2','Hbg1','Hbd'))),]
  x <- subset(x, features = rownames(counts))
  x <- subset(x, nFeature_RNA < 10000  & percent.mt < 50)
  x <- subset(x, nFeature_RNA >200 & nFeature_RNA < 6000 & nCount_RNA > 1000  & percent.mt < 20)
  return(x)
})


gc()

sceList_filtered_with_out_qc <- lapply(X = seurat.list, FUN = function(x) {
  
  table(x@meta.data$orig.ident) 
  
  # position <- substr(x@meta.data[["orig.ident"]],5,6)
  # x <- AddMetaData(
  #   object = x,
  #   metadata = position,
  #   col.name = 'position')
  # 
  # age <- substr(x@meta.data[["orig.ident"]],7,8)
  # x <- AddMetaData(
  #   object = x,
  #   metadata = age,
  #   col.name = 'age')
  
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  metadata <- x@meta.data
  metadata$cells <- rownames(metadata)
  
  rownames(x)[grepl('^mt-',rownames(x),ignore.case = T)]
  rownames(x)[grepl('^RP[SL]',rownames(x),ignore.case = T)]
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  
  rb.genes <- rownames(x)[grep("^RP[SL]",rownames(x),ignore.case = T)]
  C<-GetAssayData(object = x, slot = "counts")
  percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
  x<-AddMetaData(x, percent.ribo, col.name = "percent.ribo")
  
  hist.genes <- rownames(x)[grep("^HIST",rownames(x),ignore.case = T)]
  C<-GetAssayData(object = x, slot = "counts")
  percent.hist <- Matrix::colSums(C[hist.genes,])/Matrix::colSums(C)*100
  x<-AddMetaData(x, percent.hist, col.name = "percent.hist")
  
  counts <- GetAssayData(x, assay = 'RNA')
  counts <- counts[-(which(rownames(counts) %in% c('Hbb','Hba1','Hba2','Hbg2','Hbg1','Hbd'))),]
  # x <- subset(x, features = rownames(counts))
  # 
  # x <- subset(x, nFeature_RNA < quantile(c(x@meta.data[["nFeature_RNA"]]), 0.9)  & percent.mt < 30)
  return(x)
})



seurat_merged_with_out_qc <- merge(sceList_filtered_with_out_qc[[1]], 
                                   y = c(sceList_filtered_with_out_qc[2:length(sceList_filtered)]), project = 'gse', merge.data = TRUE)
# Results before QC

VlnPlot(seurat_merged_with_out_qc, features = c("percent.ribo", "percent.mt",'percent.hist',"nCount_RNA", "nFeature_RNA", "log10GenesPerUMI"), 
        pt.size=0,group.by = "sample")

ggsave("QC_before.pdf",width = 20,height = 10)

VlnPlot(seurat_merged_with_out_qc, features = c("percent.ribo", "percent.mt",'percent.hist',"nCount_RNA", "nFeature_RNA", "log10GenesPerUMI"), 
        pt.size=0,group.by = "sample")

ggsave("QC_before.pdf",width = 20,height = 10)


gc()

seurat_merged <- merge(sceList_filtered[[1]], 
                       y = c(sceList_filtered[2:length(sceList_filtered)]), project = 'gse', merge.data = TRUE)

# Results after QC

VlnPlot(seurat_merged, features = c("percent.ribo", "percent.mt",'percent.hist',"nCount_RNA", "nFeature_RNA", "log10GenesPerUMI"), pt.size=0,group.by = "sample")

ggsave("QC_after.pdf",width = 20,height = 10)

qs::qsave(seurat_merged,"seurat_merged.qs")
qs::qsave(seurat_merged_with_out_qc,"seurat_merged_with_out_qc.qs")
qs::qsave(sceList_filtered,"sceList_filtered.qs")
rm(seurat.list,GSE137804,sceList_filtered,seurat_merged_with_out_qc,sceList_filtered_with_out_qc,merged_seurat_filtered)
gc()


# DoubletFinder

rm(list=ls())
gc()
sceList_filtered <- qs::qread("sceList_filtered.qs")
library(DoubletFinder)
# library(future)
# plan(multisession, workers=1)

idx  <-  1:length(sceList_filtered)
# options(future.globals.maxSize=50000 * 1024^2)
# future.apply::future_lapply(sceList_filtered,function(idx){
## preprocessing

for (data in 1:length(sceList_filtered)){
  sceList_filtered[[data]] <- NormalizeData(sceList_filtered[[data]])
  sceList_filtered[[data]] <- ScaleData(sceList_filtered[[data]], verbose = FALSE)
  sceList_filtered[[data]] <- FindVariableFeatures(sceList_filtered[[data]], verbose = FALSE)
  sceList_filtered[[data]] <- RunPCA(sceList_filtered[[data]], npcs = 30, verbose = FALSE)
  sceList_filtered[[data]] <- RunUMAP(sceList_filtered[[data]], reduction = "pca", dims = 1:30)
  sceList_filtered[[data]] <- FindNeighbors(sceList_filtered[[data]], reduction = "pca", dims = 1:30)
  sceList_filtered[[data]] <- FindClusters(sceList_filtered[[data]], resolution = 0.5)
  sweep.res.list <- paramSweep_v3(sceList_filtered[[data]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sceList_filtered[[data]]@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample14@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(sceList_filtered[[data]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  sceList_filtered[[data]] <- doubletFinder_v3(sceList_filtered[[data]], PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  ## save results
  sceList_filtered[[data]]$doubFind_res = sceList_filtered[[data]]@meta.data %>% dplyr::select(contains('DF.classifications'))
  sceList_filtered[[data]]$doubFind_score = sceList_filtered[[data]]@meta.data %>% dplyr::select(contains('pANN'))
  # return(idx)
  
}
sceList_filtered_DF <- merge(sceList_filtered[[1]], 
                             y = c(sceList_filtered[2:length(sceList_filtered)]), project = "scRNA", merge.data = TRUE)

qs::qsave(sceList_filtered_DF,"sceList_filtered_df.qs")


# harmony_integrated

rm(list=ls())
gc()

source("hanshu.R")
sceList_filtered <- qs::qread("sceList_filtered_df.qs")

merged_seurat_harmony <- harmony_integ(sceList_filtered,group.by="sample",
                                       mt.pattern="^MT-",dim.use=15,mt.cutoff=5,
                                       nf.low=300,nf.high=6000,nfeatures=2000,
                                       res=0.3)

qs::qsave(merged_seurat_harmony,"merged_seurat_harmony.qs")


if (F) {

a = readRDS("integrated_donor.rds")

a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)

class(a)

all.genes <- colnames(a)

a <- ScaleData(a, verbose = FALSE)
a <- RunPCA(a, npcs = 30, verbose = FALSE,features = VariableFeatures(object = a))
a <- RunUMAP(a, reduction = "pca", dims = 1:30)
a <- FindNeighbors(a, reduction = "pca", dims = 1:30)
a <- FindClusters(a, resolution = 0.5)
DimPlot(a)

b <- ElbowPlot(a,ndims = 50,reduction="pca") 

}

merged_seurat_harmony <- qs::qread("merged_seurat_harmony.qs")

# QC after clustering


VlnPlot(merged_seurat_harmony, features = c("percent.ribo", "percent.mt",'percent.hist',"nCount_RNA", "nFeature_RNA", "log10GenesPerUMI"), pt.size=0,
        group.by = "RNA_snn_res.0.3")

ggsave("QC_after2.pdf",width = 20,height = 10)

a <- merged_seurat_harmony

diff.wilcox = FindAllMarkers(a,only.pos = T,logfc.threshold = 0.5)
write.csv(diff.wilcox,file = "diff.wilcox.csv",quote = F)

Idents(a)
ALL_MARKER = FindAllMarkers(a,only.pos = F,logfc.threshold = 0)
write.csv(ALL_MARKER,file = "ALL_MARKER.csv",quote = F)

diff.wilcox <- read.csv(file = "diff.wilcox.csv",row.names = 1)
markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.01,)
top100 = markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


marker.genes.list <- list(
  "CD14+_Monocyte" = c("CD14","CD16","CD36","CCR2","CD64","CD62L","PU1","JUN","FOS","IRF8","KLF4","VCAN","CD163","CD63","S100A12","S100A8"),
  "CD16+_Monocyte" = c("CD14","CD16","CX3CR1","CXCR4","HLA-DR","PU1","NR4A1","KLF2","FCGR3A","IFITM-3","CDKN1C","MTSS1"),
  "pDC" = c("IL3RA", "CLEC4C"),
  "Treg" = c("RTKN2","FOXP3","AC133644.2","CD4,IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2"),
  "CD16+_NK" = c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1",'GNLY','NKG7','KLRF1','TRDC'),
  "CD56+_NK" = c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A",'GNLY','NKG7','KLRF1','TRDC'),
  "naive_B" = c("IGHM","IGHD","IL4R","CXCR4","BTG1","TCL1A","YBX3"),
  "Memory_B" = c("COCH","AIM2","BANK1","SSPN","TEX9","RALGPS2","TNFRSF17","LINC01781"),
  "Plasma"=c("TNFRSF17","SDC1"),
  "CD4+_T_naive" = c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL,NOSIP","LDHB,PIK3IP1","CD3","CD3G","CD8B","IL7R","TRAC","CD3D","IL7R","MAL","LTB"),
  "CD8+_T_naive"= c("CD8B","S100B","CCR7","RGS10","NOSIP","LINC02446","LEF1","CRTAM","CD8A","OXNAD1","CD3","CD3G","IL7R","TRAC","CD3D","IL7R","MAL","LTB"),
  "cDC2" = c('FCER1A','HLA-DQA1','CLEC10A','CD1C','ENHO','PLD4','GSN','SLC38A1','NDRG2','AFF3')
)


library(dplyr)

top100$cluster <- as.factor(top100$cluster)

cluster_markers <- split(top100$gene, top100$cluster)


cluster_annotations <- list()


for (cell_type in names(marker.genes.list)) {
  
  cell_type_markers <- marker.genes.list[[cell_type]]
  
  
  for (cluster in names(cluster_markers)) {
    
    cluster_top_genes <- cluster_markers[[cluster]]
    
   
    overlap <- length(intersect(cell_type_markers, cluster_top_genes))
    
   
    if (is.null(cluster_annotations[[cluster]]) || overlap > cluster_annotations[[cluster]]$overlap) {
      
      cluster_annotations[[cluster]] <- list(cell_type = cell_type, overlap = overlap)
    }
  }
}


cluster_cell_types <- sapply(cluster_annotations, function(x) x$cell_type)


cluster_annotation_top100 <- data.frame(cluster = names(cluster_cell_types), cell_type = cluster_cell_types, stringsAsFactors = FALSE)
cluster_annotation_top100$cell_type

# cluster_annotation_top100$cell_type <- c("CD4+_T_naive","CD14+_Monocyte","CD16+_NK","CD8+_T_naive","CD16+_NK","CD16+_NK"
                                         # ,"CD8+_T_naive","CD16+_Monocyte","Memory_B","CD16+_NK","naive_B","cDC","Treg","pDC","CD56+_NK","CD14+_Monocyte","CD14+_Monocyte","CD56+_NK")

# cluster_cell_types <- c("CD4+_T_naive","CD14+_Monocyte","CD16+_NK","CD8+_T_naive","CD16+_NK","CD16+_NK","CD16+_NK","CD16+_Monocyte","Memory_B","Memory_B","naive_B","cDC","Treg","pDC","CD56+_NK","CD14+_Monocyte","CD14+_Monocyte","CD56+_NK")
cluster_cell_types <- cluster_annotation_top100$cell_type

# Change Cluster name to the most possible cell type 
cluster_cell_types[10] <- "Plasma"

a$cell_type <- cluster_cell_types[a$seurat_clusters]
table(a$cell_type)

# split the barcode list for each cell type
barcode_per_celltype <- split(names(a), a)

for (celltype in names(barcode_per_celltype)) {
  file_name <- paste0(celltype, ".txt")
  file_conn <- file(file_name)
  writeLines(barcode_per_celltype[[celltype]], con = file_conn)
  close(file_conn)
}


# cluster9 --------------------------------------------------------------
# `%notin%` <- Negate(`%in%`)
# a <- subset(a,cell_type%notin%"Plasma")

cell.type.colors <- c("CD16+_NK" = "#bf812d", 
                      "CD4+_T" = "#4a1486", 
                      "CD4+_T_naive" = "#bcbddc",
                      "CD8+_T_naive" = "#a6bddb",
                      "CD8+_T_GZMK+" = "#74a9cf",
                      "cDC2" = "#f46d43",
                      "naive_B" = "#c7e9c0",
                      # "Memory_B" = "#a1d99b",
                      "Plasma" = "#a1d99b",
                      "pDC" = "#d73027",
                      "Treg" = "#CCCC4D",
                      "CD14+_Monocyte" = "#fee090",
                      "CD16+_Monocyte" = "#fdae61",
                      "atypical_B" = "#238b45",
                      "CD16+_NK" = "#bf812d",
                      "CD56+_NK" = "#8c510a"
 
                      )
DimPlot(a,label = T)
ggsave("UMAP_DimPlot_res=0.3.pdf",width = 15,height = 10)

plot1 <- DimPlot(a, reduction="umap",group.by = "cell_type",cols = cell.type.colors,label = TRUE,pt.size = 2.5,label.size = 5)
ggsave("UMAP.pdf",plot1,width = 15,height = 10)

plot1 <- DimPlot(a, reduction="tsne",group.by = "cell_type",cols = cell.type.colors,label = TRUE,pt.size = 2.5,label.size = 5)
ggsave("TSNE.pdf",plot1,width = 15,height = 10)

barcode_per_cluster <- seurat_obj@meta.data$barcode[seurat_obj@meta.data$seurat_clusters == 1]


# bubble plot ---------------------------------------------------------------------

library(ggplot2)
library(RColorBrewer) 
library(viridis)
library(wesanderson)
pal <- wes_palette("Zissou1", 10, type = "continuous")
marker_top3 <- marker.genes.list %>% unlist()
##  RotatedAxis() scale_colour_gradientn are both important

Idents(a) <- a$cell_type

markers <- c("HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN",    "SLC38A1" , "NDRG2",  "AFF3",
  "FCGR3A", "IFITM-3","CDKN1C","FCER1G",
  "TNFRSF17",
  "CD14","VCAN","S100A12","S100A8", 
  "IRF8","IL3RA",  "CLEC4C", "GZMB","CD4",
  "IGHD","IL4R","BTG1","TCL1A",  "YBX3","COCH","AIM2","BANK1", "SSPN","TEX9","RALGPS2","LINC01781",
  "IL7R",
  "GNLY","TYROBP", "NKG7","TRDC","PRF1","FGFBP2", "SPON2", "KLRF1",  "XCL2",
  "TCF7", "CCR7",  
  "FHIT","LEF1","MAL,NOSIP","LDHB,PIK3IP1", "CD3", "CD3G","CD8B","TRAC",  
  "CD3D","MAL", "LTB", "S100B",  "RGS10",  "NOSIP",  "LINC02446","CRTAM",  "CD8A",  
  "OXNAD1")
  
#   "CD16","CD36","CCR2","CD64","CD62L",  "PU1", "JUN", "FOS",
# "KLF4","CD163",  "CD63","CX3CR1", "CXCR4", 
#  "HLA-DR", "NR4A1",  "KLF2", "MTSS1",  
#  "RTKN2",  "FOXP3",  "AC133644.2","IL2RA","TIGIT",  "CTLA4",  "FCRL3",  "LAIR2",  "IKZF2", 
# "SPINK2", "KLRC1",  "XCL1","SPTSSB", "PPP1R9A","NCAM1",  "TNFRSF11A",   
# "IGHM","SDC1", "FCER1A",  )       

DotPlot(a, features = unique(rev(markers)) , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "Cell markers", y = "", x="")

## bubble blot
ggsave("markers_bubble_plot.pdf",width = 25,height = 8)

# plot2 <- FeaturePlot(a, 
#                      reduction = "umap", 
#                      features = marker.genes.list$`CD14+_Monocyte`, 
#                      order = TRUE, 
#                      label = TRUE,
#                      min.cutoff = 'q10',
#                      label.size = 1)
# plot3 <- FeaturePlot(a, 
#             reduction = "umap", 
#             features = marker.genes.list$pDC, 
#             order = TRUE, 
#             label = TRUE,
#             min.cutoff = 'q10',
#             label.size = 1)
# plot4 <- FeaturePlot(a, 
#                      reduction = "umap", 
#                      features = c("CCR7","LEF1"), 
#                      order = TRUE, 
#                      label = TRUE,
#                      min.cutoff = 'q10',
#                      label.size = 1)
# plot5 <- FeaturePlot(a, 
#                      reduction = "umap", 
#                      features = marker.genes.list$naive_B, 
#                      order = TRUE, 
#                      label = TRUE,
#                      min.cutoff = 'q10',
#                      label.size = 1)
# plot2+plot3+plot4+plot5
gene_int<- intersect(unique(marker_top3),rownames(a))
gene_int <- c("CCR7","NKG7","IL7R","BANK1","GZMB","S100A8","TNFRSF17","FCGR3A","HLA-DQA1")
pdf("FeaturePlot.pdf",width = 12,height = 10)
# for (p in gene_int){
# p <- FeaturePlot(a, 
#            reduction = "umap",
#            features = p,
#            order = TRUE,
#            label = TRUE,
#            min.cutoff = 'q10',
#            label.size = 1)
# plot(p)
# }

FeaturePlot(a, 
                       reduction = "umap",
                       features = gene_int,
                       order = TRUE,
                       label = TRUE,
                       min.cutoff = 'q10',ncol = 3,
                       label.size = 1)
dev.off()

qs::qsave(a,file = "named.qs")

saveRDS(a,file = "named.Rds")


# ACTIONet auto-anno --------------------------------------------------------------------
rm(list = ls())
seurat_merged <- qs::qread('named.qs')


library(scater)
library(loomR)
library(Seurat)
library(patchwork)
library(ACTIONet)

# Run ACTIONet
ace <- Seurat::as.SingleCellExperiment(seurat_merged)

gc()
ace = reduce.ace(ace = ace)
ace = normalize.ace(ace)
ace = reduce.ace(ace)
depth = 20
ACTIONet_results = runACTIONet(ace = ace, k_max = depth)

# Annotate cell-types
data("curatedMarkers_human")

markers = curatedMarkers_human$Blood$PBMC$Ding2019$marker.genes
annot.out = annotate.cells.using.markers(ace = ACTIONet_results, markers = markers)
ace$celltypes = annot.out$Label

# Visualize output
# plot.ACTIONet(ace, "celltypes")
# ggsave("plot.ACTIONet.pdf",width=10,heigth=9)
# Export results as AnnData
# ACE2AnnData(ace, fname = "GSE174332.h5ad")
seurat_merged$celltypes <- ace$celltypes
qs::qsave(seurat_merged,file = "seurat_merged_ACTIONet.qs")

seurat_merged <- qs::qread(file = "seurat_merged_ACTIONet.qs")
plot1 <- DimPlot(seurat_merged, reduction="umap",group.by = "celltypes",
                 # cols = cell.type.colors,seurat_merged
                 label = TRUE,pt.size = 2,label.size = 5)
ggsave("UMAP_ACTIONet.pdf",plot1,width = 15,height = 10)


################## GSVA go ##################


sce <- qs::qread("named.qs")

sce$cell_type

av <- AverageExpression(sce,group.by = "cell_type",
                        assays = "RNA")
av <- av[[1]]
dim(av)
library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens",
                         category = "C5",subcategory = "BP")
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))

library(gplots)
library(GSVA)
library(GSEABase)
gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs <- lapply(gs, unique)

gsc <- GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
},gs,names(gs)))
gsc

geneset <- gsc
X <- av
es.max <- GSVA::gsva(X,geneset,
                     mx.diff=F,verbose=F,parallel.sz=10)
write.csv(es.max,file = "all_GSVA.csv",quote = F,row.names = T)

ROWNAMES(es.max)
pdf("GSVA.pdf",width = 20,height = 20)
pheatmap::pheatmap(es.max)
dev.off()

# Select the differences
cg <- names(tail(sort(apply(es.max, 1, sd)),50))

# pheatmap::pheatmap(es.max[cg,],show_rownames = F,cluster_cols = F)
pdf("top50_GSVA.pdf",width = 12,height = 10)
pheatmap::pheatmap(es.max[cg,],show_rownames = T,show_colnames = T,luster_cols = F,cluster_rows = T,cluster_cols = F)
dev.off()

df <- do.call(rbind,
              lapply(1:ncol(es.max), function(i){
                dat=data.frame(
                  path = rownames(es.max),
                  cluster=colnames(es.max)[i],
                  sd.1=es.max[,i],
                  sd.2=apply(es.max[,-i], 1, sd)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top <- df %>% group_by(cluster)  %>% top_n(10,fc)

gsav_pheatmap <- pheatmap::pheatmap(es.max[top$path,],show_rownames = T,cellheight = 12,cellwidth = 20,cluster_cols = F,main = "GSVA top 10"
)
ggsave(plot = gsav_pheatmap,"GSVA_TOP10.pdf",width = 15,height = 22)


if (F) {

# Only select the genes related to metabolism

# es.max_neu <- es.max[grep(("NEUR|NERV|AXON|SYNAP|CYNAP|MYLI"),rownames(es.max)),]
es.max_neu <- es.max[grep(toupper(c("metab|acid|fat|lipid|phos|biosy")),rownames(es.max)),]

# A matrix of metabolism-related genes was derived
write.csv(es.max_neu,file = "GSVA_gobp_metab.csv",quote = F)

# Derived gene concentration
all_gene_sets_int <- subset(all_gene_sets,gs_name%in%rownames(es.max_neu))
write.csv(all_gene_sets_int[,1:4],file = "GSVA_gobp_metab_all_genes.csv",quote = F)

# Select the differentials
cg <- names(tail(sort(apply(es.max_neu, 1, sd)),20))

pheatmap::pheatmap(es.max_neu[cg,],show_rownames = F)


df <- do.call(rbind,
              lapply(1:ncol(es.max_neu), function(i){
                dat=data.frame(
                  path = rownames(es.max_neu),
                  cluster=colnames(es.max_neu)[i],
                  sd.1=es.max_neu[,i],
                  sd.2=apply(es.max_neu[,-i], 1, sd)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top <- df %>% group_by(cluster)  %>% top_n(10,fc)

gsav_pheatmap <- pheatmap::pheatmap(es.max_neu[top$path,],show_rownames = T,cellheight = 12,cellwidth = 20,cluster_cols = F,
                                    main = "GSVA top 10 metabolism"
)
ggsave(plot = gsav_pheatmap,"GSVA_TOP10_metab.pdf",width = 12,height = 15)

}


### GSVA kegg 

av <- AverageExpression(sce,group.by = "cell_type",
                        assays = "RNA")
av <- av[[1]]
dim(av)
library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens",
                         category = "C2")
all_gene_sets <- filter(all_gene_sets,grepl("^KEGG",all_gene_sets$gs_name))
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))

library(gplots)
library(GSVA)
library(GSEABase)
gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs <- lapply(gs, unique)

gsc <- GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
},gs,names(gs)))
gsc

geneset <- gsc
X <- av
es.max <- GSVA::gsva(X,geneset,
                     mx.diff=F,verbose=F,parallel.sz=10)
write.csv(es.max,file = "all_GSVA_kegg.csv",quote = F,row.names = T)

ROWNAMES(es.max)
pdf("GSVA_kegg.pdf",width = 20,height = 20)
pheatmap::pheatmap(es.max)
dev.off()

# Select the differences
cg <- names(tail(sort(apply(es.max, 1, sd)),50))

# pheatmap::pheatmap(es.max[cg,],show_rownames = F,cluster_cols = F)
pdf("top50_GSVA_kegg.pdf",width = 12,height = 10)
pheatmap::pheatmap(es.max[cg,],show_rownames = T,show_colnames = T,luster_cols = F,cluster_rows = T,cluster_cols = F)
dev.off()

df <- do.call(rbind,
              lapply(1:ncol(es.max), function(i){
                dat=data.frame(
                  path = rownames(es.max),
                  cluster=colnames(es.max)[i],
                  sd.1=es.max[,i],
                  sd.2=apply(es.max[,-i], 1, sd)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top <- df %>% group_by(cluster)  %>% top_n(10,fc)

gsav_pheatmap <- pheatmap::pheatmap(es.max[top$path,],show_rownames = T,cellheight = 12,cellwidth = 20,cluster_cols = F,main = "GSVA top 10"
)
ggsave(plot = gsav_pheatmap,"GSVA_TOP10_kegg.pdf",width = 15,height = 25)



# Metabolism-related genes
if (F) {
 
# es.max_neu <- es.max[grep(("NEUR|NERV|AXON|SYNAP|CYNAP|MYLI"),rownames(es.max)),]
es.max_neu <- es.max[grep(toupper(c("metab|acid|fat|lipid|phos|biosy")),rownames(es.max)),]
# Derive a matrix of metabolism-related genes 
write.csv(es.max_neu,file = "GSVA_kegg_metab.csv",quote = F)
# Derived gene concentration
all_gene_sets_int <- subset(all_gene_sets,gs_name%in%rownames(es.max_neu))
write.csv(all_gene_sets_int[,1:4],file = "GSVA_kegg_metab_all_genes.csv",quote = F)
# Difference selection
cg <- names(tail(sort(apply(es.max_neu, 1, sd)),20))

pheatmap::pheatmap(es.max_neu[cg,],show_rownames = F)


df <- do.call(rbind,
              lapply(1:ncol(es.max_neu), function(i){
                dat=data.frame(
                  path = rownames(es.max_neu),
                  cluster=colnames(es.max_neu)[i],
                  sd.1=es.max_neu[,i],
                  sd.2=apply(es.max_neu[,-i], 1, sd)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top <- df %>% group_by(cluster)  %>% top_n(10,fc)

gsav_pheatmap <- pheatmap::pheatmap(es.max_neu[top$path,],show_rownames = T,cellheight = 12,cellwidth = 20,cluster_cols = F,
                                    main = "GSVA top 10 metabolism"
)
ggsave(plot = gsav_pheatmap,"GSVA_TOP10_metab_kegg.pdf",width = 12,height = 15)

}


################## GSVA C8 ##################

av <- AverageExpression(sce,group.by = "cell_type",
                        assays = "RNA")
av <- av[[1]]
dim(av)
library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens",
                         category = "C8")
# all_gene_sets <- filter(all_gene_sets,grepl("^KEGG",all_gene_sets$gs_name))
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))

library(gplots)
library(GSVA)
library(GSEABase)
gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs <- lapply(gs, unique)

gsc <- GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
},gs,names(gs)))
gsc

geneset <- gsc
X <- av
es.max <- GSVA::gsva(X,geneset,
                     mx.diff=F,verbose=F,parallel.sz=10)
write.csv(es.max,file = "all_GSVA_C8.csv",quote = F,row.names = T)

ROWNAMES(es.max)
pdf("GSVA_C8.pdf",width = 20,height = 20)
pheatmap::pheatmap(es.max)
dev.off()

# Difference selection
cg <- names(tail(sort(apply(es.max, 1, sd)),50))

# pheatmap::pheatmap(es.max[cg,],show_rownames = F,cluster_cols = F)
pdf("top50_GSVA_C8.pdf",width = 12,height = 10)
pheatmap::pheatmap(es.max[cg,],show_rownames = T,show_colnames = T,luster_cols = F,cluster_rows = T,cluster_cols = F)
dev.off()

df <- do.call(rbind,
              lapply(1:ncol(es.max), function(i){
                dat=data.frame(
                  path = rownames(es.max),
                  cluster=colnames(es.max)[i],
                  sd.1=es.max[,i],
                  sd.2=apply(es.max[,-i], 1, sd)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top <- df %>% group_by(cluster)  %>% top_n(10,fc)

gsav_pheatmap <- pheatmap::pheatmap(es.max[top$path,],show_rownames = T,cellheight = 12,cellwidth = 20,cluster_cols = F,main = "GSVA top 10"
)
ggsave(plot = gsav_pheatmap,"GSVA_TOP10_C8.pdf",width = 15,height = 22)



# GO kegg----------------------------------------------------------------------
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
#showCategoryÖž¶šÕ¹ÊŸµÄGO TermsµÄžöÊý£¬Ä¬ÈÏÎª10£¬ŒŽp.adjust×îÐ¡µÄ10žö
library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)

all_gene <- FindAllMarkers(sce,logfc.threshold = 0,only.pos = F)


# all_gene <- read.csv('ALL_MARKER.csv')
# all_gene <- subset(all_gene,avg_log2FC>0.585&p_val_adj>0.05)
for (i in unique(all_gene$cluster)){
  geneList <- subset(all_gene,cluster==i&(avg_log2FC>0.25&p_val_adj<0.05))
  gene.df <- clusterProfiler::bitr(geneList$gene, fromType = "SYMBOL", #fromTypeÊÇÖžÄãµÄÊýŸÝIDÀàÐÍÊÇÊôÓÚÄÄÒ»ÀàµÄ
                                   toType = c("ENTREZID",'ENSEMBL'), #toTypeÊÇÖžÄãÒª×ª»»³ÉÄÄÖÖIDÀàÐÍ£¬¿ÉÒÔÐŽ¶àÖÖ£¬Ò²¿ÉÒÔÖ»ÐŽÒ»ÖÖ
                                   OrgDb = org.Hs.eg.db)#OrgdbÊÇÖž¶ÔÓŠµÄ×¢ÊÍ°üÊÇÄÄžö
  
  ENSEMBL <- as.character(gene.df$ENSEMBL)
  gene <- as.character(geneList$gene)
  
  ###   GO 
  
  if (T) {
    
    
    
    
    GO <- enrichGO(gene          = ENSEMBL, #²îÒì»ùÒò vector
                   keyType       = "ENSEMBL",  #²îÒì»ùÒòµÄ ID ÀàÐÍ£¬ÐèÒªÊÇ OrgDb Ö§³ÖµÄ
                   OrgDb         = org.Hs.eg.db, #¶ÔÓŠµÄOrgDb
                   ont           = "ALL", #GO ·ÖÀàÃû³Æ£¬CC BP MF
                   pvalueCutoff  = 0.05, #Pvalue ãÐÖµ   £špvalue=1ÖžÊä³öËùÓÐœá¹û£¬pvalue=0.05ÖžÊä³ö·ûºÏÒªÇóµÄœá¹û£©
                   qvalueCutoff  = 0.05, #qvalue ãÐÖµ
                   pAdjustMethod = "BH", #Pvalue œÃÕý·œ·š
                   readable      = FALSE) #TRUE ÔòÕ¹ÊŸSYMBOL£¬FALSE ÔòÕ¹ÊŸÔ­ÀŽµÄID£šÑ¡falseÊÇÒòÎª²»ÊÇËùÓÐgene¶ŒÓÐsymbolµÄ£©
    
    GO_b <- GO
    # GO@result <- GO@result[grep(("neur|nerv|axon|synap|cynap|myli"),GO@result[["Description"]]),]
    
    
    BP <- enrichGO(gene          = ENSEMBL, #²îÒì»ùÒò vector
                   keyType       = "ENSEMBL",  #²îÒì»ùÒòµÄ ID ÀàÐÍ£¬ÐèÒªÊÇ OrgDb Ö§³ÖµÄ
                   OrgDb         = org.Hs.eg.db, #¶ÔÓŠµÄOrgDb
                   ont           = "BP", #GO ·ÖÀàÃû³Æ£¬CC BP MF
                   pvalueCutoff  = 0.05, #Pvalue ãÐÖµ   £špvalue=1ÖžÊä³öËùÓÐœá¹û£¬pvalue=0.05ÖžÊä³ö·ûºÏÒªÇóµÄœá¹û£©
                   qvalueCutoff  = 0.05, #qvalue ãÐÖµ
                   pAdjustMethod = "BH", #Pvalue œÃÕý·œ·š
                   readable      = FALSE) #TRUE ÔòÕ¹ÊŸSYMBOL£¬FALSE ÔòÕ¹ÊŸÔ­ÀŽµÄID£šÑ¡falseÊÇÒòÎª²»ÊÇËùÓÐgene¶ŒÓÐsymbolµÄ£©
    
    
    BP
    dir.create("GO")
    pdf(paste0("GO/cluster_",i,"_top10 GO terms of BP.pdf"),height = 6,width = 6)
    p <- dotplot(BP,title=paste0("cluster_",i,"_top10 GO terms of BP"),
                 showCategory=10)
    # p <- barplot(GO,showCategory = 15) 
    plot(p)
    dev.off()
    
    
    BP_b <- BP
    # BP@result <- BP@result[grep(("neur|nerv|axon|synap|cynap|myli"),BP@result[["Description"]]),]
    pdf(paste0("GO/cluster ",i,"_top10 GO terms of BP.pdf"),height = 6,width = 6)
    p <- dotplot(BP,title=paste0("cluster ",i,"_top10 GO terms of BP"),
                 showCategory=10)
    # p <- barplot(GO,showCategory = 15) 
    plot(p)
    dev.off()
    
    
    # GO_BP <- clusterProfiler::plotGOgraph(BP)
    
    # write.csv(BP,paste0("BP.csv"))
    
    # GO_filtered <- gofilter(GO, level=10)
    
    # GO_simp <- simplify(GO,cutoff=0.6, by="p.adjust", select_fun=min)
    
    # barplot(BP_simp, showCategory=15, color = 'p.adjust')
    
    # dotplot(GO_simp,showCategory=15,title=paste0("cluster_",i,'_top55 GO terms of each sub-class'))
    # write.csv(ego,'1w.go.bp.csv')
    
    pdf(paste0("GO/cluster_",i,"_top10 GO terms of each sub-class.pdf"),height = 6,width = 6)
    p <- dotplot(GO,title=paste0("cluster_",i,"_top10 GO terms of each sub-class"),
                 showCategory=10,split='ONTOLOGY')+
      facet_grid(ONTOLOGY~.,scale="free")
    # p <- barplot(GO,showCategory = 15) 
    plot(p)
    dev.off()
  }
  # browseVignettes("DOSE")
  
  # groupGO
  # In clusterProfiler, groupGO is genetically classified according to the distribution of GO at a particular level.
  if (T) {
    
    dir.create("groupGO")
    ggo <- groupGO(gene     = ENSEMBL,
                   OrgDb    = org.Hs.eg.db,
                   keyType  = "ENSEMBL",
                   ont      = "BP",
                   level    = 3,
                   readable = TRUE)
    
    ggodf <- data.frame(GeneRatio=ggo@result[["GeneRatio"]],description=ggo@result[["Description"]])
    ggodf <- subset(ggodf,GeneRatio!=paste0("0/",nrow(gene.df)))
    ggodf <- ggodf[order(ggodf$GeneRatio,decreasing = T),]
    write_csv(ggodf,file = paste0("groupGO/cluster_",i,"_BP_groupGO.csv"))
    
    ggo <- groupGO(gene     = ENSEMBL,
                   OrgDb    = org.Hs.eg.db,
                   keyType  = "ENSEMBL",
                   ont      = "MF",
                   level    = 3,
                   readable = TRUE)
    ggodf <- data.frame(GeneRatio=ggo@result[["GeneRatio"]],description=ggo@result[["Description"]])
    ggodf <- subset(ggodf,GeneRatio!=paste0("0/",nrow(gene.df)))
    ggodf <- ggodf[order(ggodf$GeneRatio,decreasing = T),]
    write_csv(ggodf,file = paste0("groupGO/cluster_",i,"_MF_groupGO.csv"))
    
    ggo <- groupGO(gene     = ENSEMBL,
                   OrgDb    = org.Hs.eg.db,
                   keyType  = "ENSEMBL",
                   ont      = "CC",
                   level    = 3,
                   readable = TRUE)
    
    ggodf <- data.frame(GeneRatio=ggo@result[["GeneRatio"]],description=ggo@result[["Description"]])
    ggodf <- subset(ggodf,GeneRatio!=paste0("0/",nrow(gene.df)))
    ggodf <- ggodf[order(ggodf$GeneRatio,decreasing = T),]
    write_csv(ggodf,file = paste0("groupGO/cluster_",i,"_CC_groupGO.csv"))
    
  }
  
  
  # head(ggo,3)
  # clusterProfiler::groupGO(ggo)
  # ggo
  # 
  # kegg --------------------------------------------------------------------
  
  if (T) {
    # gene <- read.csv('esc_res1.0.markers.csv')
    # geneList <- gene[251:792,]
    ## feature 2: named vector
    # geneList <- as.character(geneList$gene)
    
    # gene.df <- bitr(geneList, fromType = "SYMBOL", #fromTypeÊÇÖžÄãµÄÊýŸÝIDÀàÐÍÊÇÊôÓÚÄÄÒ»ÀàµÄ
    #                 toType = c("ENTREZID",'ENSEMBL'), #toTypeÊÇÖžÄãÒª×ª»»³ÉÄÄÖÖIDÀàÐÍ£¬¿ÉÒÔÐŽ¶àÖÖ£¬Ò²¿ÉÒÔÖ»ÐŽÒ»ÖÖ
    #                 OrgDb = org.Mm.eg.db)#OrgdbÊÇÖž¶ÔÓŠµÄ×¢ÊÍ°üÊÇÄÄžö
    # gene <- as.character(gene.df$ENTREZID)
    
    # Every time R is opened, it automatically links to kegg's website for the most recent species annotation information, so the database must be up to date
    kegg <- enrichKEGG(
      gene = gene.df$ENTREZID, 
      keyType = 'kegg',  
      organism = 'hsa',  
      pAdjustMethod = 'fdr', 
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.75 
      
    )
    # kegg_b <- kegg
    # kegg@result <- kegg@result[grep(("neur|nerv|axon|synap|cynap|myli"),kegg@result[["Description"]]),]
    
    dir.create("kegg")
    write.table(kegg, paste0("kegg/cluster_",i,'_kegg.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
    # if (!kegg@result[["qvalue"]]>"0.05") {
    #   GO@organism
    #       }
    pdf(paste0("kegg/cluster_",i,"_top10 kegg terms_dot.pdf"),height = 10,width = 8)
    # dotplot(kegg,title=paste0("cluster_",i,"_top55 kegg terms"),
    #         showCategory=15)
    p <-  clusterProfiler::dotplot(kegg,showCategory = 10)
    # p
    plot(p)
    dev.off()

    # }
  }
  
  
  # MSigDb analysis ---------------------------------------------------------
  if (T) {
    
    
    dim(geneList)
    
    # Ms_msigdbr <- msigdbr(species="Mus musculus")
    library(msigdbr)
    MsGO <- msigdbr(category="C8"
                    # ,subcategory = "BP"
    ) %>% 
      dplyr::select(gs_name, entrez_gene, gene_symbol)
    head(MsGO)
    
    em <- enricher(gene,TERM2GENE=MsGO[,c(1,3)])
    
    head(em,5)
    dir.create("MSigDb_C8")
    write.csv(em, paste0("MSigDb_C8/cluster_",i,'_msigdbr.csv'))
    pdf(paste0("MSigDb_C8/cluster_",i,'_msigdbr.pdf'))
    p1<- barplot(em,showCategory=10)
    
    p2 <- cnetplot(em, showCategory = 5)
    plot(p1)
    plot(p1)
    dev.off()
    # goplot(em, showCategory = 5)
    
    
    
    # em1 <- GSEA(gene,TERM2GENE=MsGO[,c(1,3)])
    # head(em1,1)
    # dim(em)
    
    # dotplot(em,showCategory=10)
    # goplot(em)
    
  }
  
}


# Cell cycle 
merged_seurat_harmony <- qs::qread("merged_seurat_harmony.qs")


# mouse.cc.gene=readRDS("../mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
merged_seurat_harmony=CellCycleScoring(object = merged_seurat_harmony, 
                              g2m.features = cc.genes$g2m.genes, 
                              s.features = cc.genes$s.genes)
p4=VlnPlot(merged_seurat_harmony, features = c("S.Score", "G2M.Score"), group.by = "RNA_snn_res.0.3", 
           ncol = 2, pt.size = 0.1)

pdf("Vlnplot_CellCycleScoring.pdf",width = 10,height = 10)
p4
dev.off()


pdf("11.cycle_details.pdf")
sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
dev.off()
# Those with high S.core are S stage, G2M stage. Those with high Score are G2M stage, and those with relatively low score are G1 stage

# cluster 9 ---------------------------------------------------------------

merged_seurat_harmony <- qs::qread("merged_seurat_harmony.qs")

# To demonstrate whether cluster9 is specific for its high expression genes
pdf("cluster9_marker.pdf",width = 10,height = 10)
FeaturePlot(merged_seurat_harmony,features = c("JCHAIN","MZB1","TXNDC5","STMN1"),ncol = 2)
VlnPlot(merged_seurat_harmony,features = c("JCHAIN","MZB1","TXNDC5","STMN1"),ncol = 2)
DotPlot(merged_seurat_harmony,features = c("JCHAIN","MZB1","TXNDC5","STMN1"))
dev.off()
# If it's Plasma Markers there are BCMA(TNFRSF17) and CD138 (SDC1)

pdf("cluster9_plasma_marker.pdf",width = 10,height = 5)
FeaturePlot(merged_seurat_harmony,features = c("TNFRSF17","SDC1"),reduction = "umap")
dev.off()

# If it's Memory B
pdf("cluster9_memory_B_marker.pdf",width = 10,height = 5)
FeaturePlot(merged_seurat_harmony,features = c("COCH","AIM2","BANK1","SSPN","TEX9","RALGPS2","TNFRSF17","LINC01781"),reduction = "umap")
dev.off()
"Memory_B" = c("COCH","AIM2","BANK1","SSPN","TEX9","RALGPS2","TNFRSF17","LINC01781")