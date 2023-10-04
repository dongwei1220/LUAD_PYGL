### LUAD单细胞数据分析 ###
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)

setwd("/data/LUAD_PYGL/GSE131907")

# 读取count数据
data <- readRDS("/data/LUAD_PYGL/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
dim(data)
#[1]  29634 208506
data[1:5,1:5]
#         AAACCTGAGAAACCGC_LN_05 AAACCTGAGAAACGCC_NS_13
#A1BG                          0                      0
#A1BG-AS1                      0                      0
#A1CF                          0                      0

# 读取meta样本信息
meta <- read.table("/data/LUAD_PYGL/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",check.names=F,header=T,row.names=1)
head(meta)
#                                   Barcode   Sample Sample_Origin     Cell_type
#AAACCTGCAAGGTGTG_LUNG_N01 AAACCTGCAAGGTGTG LUNG_N01         nLung Myeloid cells
#AACTCCCGTTCACCTC_LUNG_N01 AACTCCCGTTCACCTC LUNG_N01         nLung Myeloid cells
#AACTCCCTCACGCGGT_LUNG_N01 AACTCCCTCACGCGGT LUNG_N01         nLung Myeloid cells
dim(meta)
#[1] 208506      6

# 构建seurat对象
data_obj <- CreateSeuratObject(counts = data, 
                               meta.data=meta, project = "LUAD",
                               min.cells = 3, min.features = 200)
data_obj
#An object of class Seurat
#27578 features across 208506 samples within 1 assay
#Active assay: RNA (27578 features, 0 variable features)

table(data_obj$Sample_Origin)
#mBrain    mLN    nLN  nLung     PE   tL/B  tLung
# 29060  21479  37446  42995  20304  12073  45149
table(data_obj$Cell_type)
#    B lymphocytes Endothelial cells  Epithelial cells       Fibroblasts
#            27657              1996             36467              4172
#       MAST cells     Myeloid cells          NK cells  Oligodendrocytes
#             3396             42245             11551               716
#    T lymphocytes      Undetermined
#            79676               630

# 计算线粒体基因含量``
data_obj <- PercentageFeatureSet(data_obj, "^MT-", col.name = "percent_mito")

# subset tumor samples("tLung","tL/B","mLN","mBrain")
Idents(data_obj) <- "Sample_Origin"
data_obj <- subset(data_obj,idents=c("tLung","tL/B","mLN","mBrain"))

### integration with seurat SCTransform pipeline ###
data_list <- SplitObject(data_obj, split.by = "Sample_Origin")
data_list

# 分别对每个对象进行SCTransform标准化处理
for (i in names(data_list)) {
    data_list[[i]] <- SCTransform(data_list[[i]], vars.to.regress=c("percent_mito"), verbose = FALSE)
}

data_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000, verbose = FALSE)
head(data_features)
#[1] "IGHG1" "IGKC"  "SPP1"  "IGHG4" "GNLY"  "IGLC2"
length(data_features)
#[1] 3000

data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = data_features, verbose = FALSE)
data_list

# integrate data with cca method
data_anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = data_features, verbose = FALSE, reduction = "cca")
data_anchors
#An AnchorSet object containing 582720 anchors between 6 Seurat objects
# This can be used as input to IntegrateData.

#options(future.rng.onMisuse="ignore")
data_integrated <- IntegrateData(anchorset = data_anchors, normalization.method = "SCT", verbose = FALSE)
data_integrated

DefaultAssay(data_integrated)
#[1] "integrated"

# pca dimensional reduction
data_integrated <- RunPCA(object = data_integrated, verbose = FALSE, npcs = 50)

pdf(file='All_sample_seurat_integrated_ElbowPlot.pdf',height=5,width=6.5)
ElbowPlot(data_integrated, ndims=50)
dev.off()
pdf(file='All_sample_seurat_integrated_PCAPlot.pdf',height=5,width=6)
DimPlot(data_integrated, reduction = "pca", group.by="Sample_Origin")
dev.off()

# Run non-linear dimensional reduction
#The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells
#together in low-dimensional space.
#sample_colors <- rev(c("#a6cee3","#1f78b4","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#b15928"))
sample_colors <- as.vector(ArchR::ArchRPalettes$circus)

data_integrated <- RunUMAP(object = data_integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pdf(file='All_sample_seurat_integrated_UMAP_dimplot_by_sample.pdf',height=5,width=6)
DimPlot(data_integrated, reduction = "umap", group.by = "Sample_Origin") + scale_color_manual(values=sample_colors)
dev.off()

# cell clustering
data_integrated <- FindNeighbors(data_integrated, reduction = "pca", dims = 1:30)
data_integrated <- FindClusters(data_integrated, resolution = 0.3)

#cluster_colors <- c("#5A2955", "#D51F26", "#753A80", "#3E5D8D", "#2B7F4A", "#C0545B", "#DCABCF", "#FDDD03", "#F69421", "#9C82BA", "#AA875A", "#3EB3A7", "#245359", "#711E21", "#B36B45", "#DB7D8D", "#BF5D58", "#8ED2D0", "#9A7456", "#216069", "#3D3D3D", "#B4B883", "#C39C6F", "#9B8EC4", "#6264A0", "#BFCADF", "#8AC972", "#D51F26", "#272E6A", "#208A42")
cluster_colors <- c(as.vector(ArchR::ArchRPalettes$stallion),"#BFCADF", "#8AC972", "#D51F26", "#272E6A", "#208A42","#B4B883", "#C39C6F", "#9B8EC4", "#6264A0","#BF5D58", "#8ED2D0", "#9A7456")

pdf(file='All_sample_seurat_integrated_UMAP_dimplot_by_cluster.pdf',height=5,width=6.2)
DimPlot(data_integrated, reduction = "umap", label = T, label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()
pdf(file='All_sample_seurat_integrated_UMAP_dimplot_by_celltype1.pdf',height=5,width=6.5)
DimPlot(data_integrated, reduction = "umap", label = F, group.by= "Cell_type", label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()
pdf(file='All_sample_seurat_integrated_UMAP_dimplot_by_celltype2.pdf',height=5,width=6.5)
DimPlot(data_integrated, reduction = "umap", label = F, group.by = "Cell_type.refined", label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()
pdf(file='All_sample_seurat_integrated_UMAP_dimplot_by_celltype3.pdf',height=5,width=11)
DimPlot(data_integrated, reduction = "umap", label = F, group.by= "Cell_subtype", label.size = 5)
dev.off()
pdf(file='All_sample_seurat_integrated_UMAP_dimplot_split_by_sample.pdf',height=5.5,width=20)
DimPlot(data_integrated, reduction = "umap", label = F, split.by="Sample_Origin",label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()

data_integrated
#An object of class Seurat
#55880 features across 208506 samples within 3 assays
#Active assay: integrated (3000 features, 3000 variable features)
# 2 other assays present: RNA, SCT
# 3 dimensional reductions calculated: pca, umap, tsne

DefaultAssay(data_integrated) <- "SCT"

pdf(file='FeaturePlot_by_PYGL.pdf',height=5,width=5.5)
FeaturePlot(data_integrated, reduction = "umap", features = "PYGL", order = T, min.cutoff="q5",max.cutoff="q95",cols = c("grey", "red"))
dev.off()

pdf(file='DotPlot_by_PYGL_MHCII.pdf',height=3.5,width=4.5)
DotPlot(subset(data_integrated,idents=c(2,5,6,9,13)), features = c("PYGL","HLA-DPA1","HLA-DPB1","HLA-DOA","HLA-DMA","HLA-DMB","HLA-DQA1","CIITA","CD276","CD274","PDCD1","HAVCR2"), col.min=0, col.max=1, scale.max = 40, cluster.idents=T) + coord_flip() + scale_color_viridis(direction=1)
dev.off()

#saveRDS(data_integrated,"All_sample_seurat_integrated_LUAD.rds")

###############################################################################
# sample-cluster stats
sample_stats <- table(data_integrated$Sample_Origin, data_integrated$Cell_type)
library(reshape2)
sample_stats <- melt(sample_stats)
colnames(sample_stats) <- c("Sample","Cell_type","Value")
head(sample_stats)
#  Sample     Cell_type Value
#1 mBrain B lymphocytes  1311
#2    mLN B lymphocytes  6062
#3    nLN B lymphocytes 10584
#4  nLung B lymphocytes   634
#5     PE B lymphocytes  3285
#6   tL/B B lymphocytes   469

sample_colors <- as.vector(ArchR::ArchRPalettes$circus)

pdf("All_sample_celltype_stats_by_celltype.pdf",height=7,width=6)
ggplot(sample_stats,aes(Cell_type,Value,fill=Sample)) + geom_bar(stat="identity",position="fill") + geom_hline(yintercept=0.5,color="black",linetype=2) + scale_fill_manual(values=sample_colors) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

pdf("All_sample_celltype_stats_by_sample.pdf",height=6.5,width=6)
ggplot(sample_stats,aes(Sample,Value,fill=Cell_type)) + geom_bar(stat="identity",position="fill") + scale_fill_manual(values=cluster_colors) + theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

### DE analysis by celltype ###
DefaultAssay(data_integrated) <- "RNA"
data_integrated <- NormalizeData(data_integrated)

dir_DE_All_Cluster <- "DE_All_CellType"
dir.create(dir_DE_All_Cluster)

Idents(data_integrated) <- "Cell_type"

# DE analysis by celltype
FindDE_all_cluster <- FindAllMarkers(data_integrated, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)

# filtering
# p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
# filter RPLXX (Ribosomal protein) and MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

output=paste0(dir_DE_All_Cluster, "/", "DE_all.celltype.pValAdj_0.05.csv")
write.csv(FindDE_all_cluster, file=output, quote=F)

DefaultAssay(data_integrated) <- "SCT"

library(dplyr)
library(viridis)
library(RColorBrewer)

top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(paste0(dir_DE_All_Cluster, "/", "DE_all.celltype.top10.heatmap.pdf"), width=24, height=15)
#p <- DoHeatmap(data_integrated, features = top5$gene) + NoLegend()
p <- DoHeatmap(data_integrated, features = as.vector(top10$gene), group.colors=cluster_colors) + scale_fill_viridis()
print(p)
dev.off()

top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf(paste0(dir_DE_All_Cluster, "/","DotPlot_all_cluster_top5_markers.pdf"),height=5,width=15)
DotPlot(object = data_integrated, features=unique(top5$gene), assay="RNA") + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

### DE analysis between C6-vs-C13
dir_DE_All_Cluster <- "DE_All_Cluster"
dir.create(dir_DE_All_Cluster)
DefaultAssay(data_integrated) <- "RNA"
Idents(data_integrated) <- "seurat_clusters"

FindDE_cluster <- FindMarkers(data_integrated, ident.1 = "6", ident.2 = "13", logfc.threshold = 0, min.pct = 0)
output=paste0(dir_DE_All_Cluster, "/", "DE_C6-vs-C13.all.csv")
write.csv(FindDE_cluster, file=output, quote=F)

# gene filtering
# logfc >= 0.25 and p val adjust < 0.05
FindDE_cluster <- FindDE_cluster[FindDE_cluster$avg_log2FC >= 0.25 & FindDE_cluster$p_val_adj < 0.05, ]
# filter RPLXX (Ribosomal protein) and MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_cluster), perl=T)
if(length(tmp) > 0){FindDE_cluster <- FindDE_cluster[-tmp, ]}

output=paste0(dir_DE_All_Cluster, "/", "DE_C6-vs-C13.log2fc_0.25.pValAdj_0.05.csv")
write.csv(FindDE_cluster, file=output, quote=F)

### subset Epithelial cells for reclustering analysis ###
Idents(data_integrated) <- "seurat_clusters"
table(Idents(data_integrated))
#    0     1     2     3     4     5     6     7     8     9    10    11    12
#19472 14826 12651 11404  9930  9484  5756  5566  3208  3058  2577  2488  2309
#   13    14    15    16    17    18    19    20
# 1269   896   746   648   561   443   303   166

#'Epithelial cells': 2,5,6,9,13
data_epi <- subset(data_integrated,idents=c(2,5,6,9,13))
data_epi
#An object of class Seurat
#55482 features across 32218 samples within 3 assays
#Active assay: RNA (27578 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

pdf(file='Epithelial_reclustering_UMAP_dimplot_by_sample.pdf',height=5,width=5.5)
DimPlot(data_epi, reduction = "umap", group.by = "Sample_Origin") + scale_color_manual(values=cluster_colors)
dev.off()
pdf(file='Epithelial_reclustering_UMAP_dimplot_by_cluster.pdf',height=5,width=5.5)
DimPlot(data_epi, reduction = "umap", label = F, label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()
pdf(file='Epithelial_reclustering_UMAP_dimplot_split_by_sample.pdf',height=5.5,width=20)
DimPlot(data_epi, reduction = "umap", label = F, split.by="Sample_Origin",label.size = 5) + scale_color_manual(values=cluster_colors)
dev.off()

# DE analysis by celltype
dir_DE_All_Cluster <- "DE_All_Cluster"
DefaultAssay(data_epi) <- "RNA"

FindDE_all_cluster <- FindAllMarkers(data_epi, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
# filtering
# p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
# filter RPLXX (Ribosomal protein) and MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

output=paste0(dir_DE_All_Cluster, "/", "DE_epithelial.clusters.pValAdj_0.05.csv")
write.csv(FindDE_all_cluster, file=output, quote=F)

DefaultAssay(data_epi) <- "SCT"

library(dplyr)
library(viridis)
library(RColorBrewer)

top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(paste0(dir_DE_All_Cluster, "/", "DE_epithelial.clusters.top10.heatmap.pdf"), width=20, height=12)
p <- DoHeatmap(data_epi, features = as.vector(top10$gene), group.colors=cluster_colors, draw.lines = F) + scale_fill_viridis()
print(p)
dev.off()

top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf(paste0(dir_DE_All_Cluster, "/","DotPlot_epithelial.clusters_top5_markers.pdf"),height=3.5,width=8.5)
DotPlot(object = data_epi, features=unique(top5$gene), assay="RNA") + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file='FeaturePlot_by_PYGL_MHCII.pdf',height=15,width=16.5)
FeaturePlot(data_epi, reduction = "umap", features = c("PYGL","HLA-DPA1","HLA-DPB1","HLA-DOA","HLA-DMA","HLA-DMB","HLA-DQA1","CD74","CIITA"), order = T, min.cutoff="q5",max.cutoff="q95",cols = c("grey", "red"))
dev.off()

data_epi
#An object of class Seurat
#55482 features across 32218 samples within 3 assays
#Active assay: SCT (24904 features, 0 variable features)
# 2 other assays present: RNA, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

saveRDS(data_epi,"data_epi.rds")


