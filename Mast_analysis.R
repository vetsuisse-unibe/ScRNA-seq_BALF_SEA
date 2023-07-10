# INDEPENDENT ANALYSIS OF MAST CELLS

# After subsetting mast cells from the dataset,
# we rerun the following steps: sctransform, integration, PCA and clustering. 

# Set up environment
set.seed(42)
library(Seurat)
library(clustree)
library(celldex)
library(SingleR)
library(cowplot)
library(patchwork)
library(dplyr)
library(Matrix)
library(gdata)
library(ggtext)
library(kableExtra)
library(ggplot2)
require("biomaRt")
library(RColorBrewer)
library(stringr)
library(viridis)
library(BiocParallel)
library(scDblFinder)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(scFeatureFilter)
library(glmGamPoi)
library(scran)
library(PCAtools)
library(scSorter)
library(metap)
library(tibble)
library(speckle)


# Retrieve the Seurat object
sc22 <- readRDS("sc22_seed.rds")
# Set the grouping ident to the 6 major cell types identified
Idents(sc22) <- "MajCellType"

# Subset the mast cells from the complete data set 
sc22_mast <- subset(sc22, idents = "Mast cells") 
DimPlot(sc22_mast) 

# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_mast@meta.data[,grep("integrated_snn_res", colnames(sc22_mast@meta.data))] <- NULL

# NORMALIZATION
sc22_mast <- SCTransform(sc22_mast, method = "glmGamPoi", variable.features.n = 3000)

vf_top10_mast <- head(VariableFeatures(sc22_mast), 10)
vf_plot_mast <- VariableFeaturePlot(sc22_mast)
LabelPoints(plot = vf_plot_mast,
            points = vf_top10_mast, repel = TRUE, xnudge = 0, ynudge = 0)

#Perform dimensionality reduction by PCA and UMAP embedding
sc22_mast <- RunPCA(sc22_mast, verbose = FALSE)
sc22_mast  <- RunUMAP(sc22_mast, dims = 1:30, verbose = FALSE)
sc22_mast <- FindNeighbors(sc22_mast, dims = 1:30, verbose = FALSE)
sc22_mast  <- FindClusters(sc22_mast, verbose = FALSE)
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_mast <- CellCycleScoring(sc22_mast, 
                            s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                            set.ident = TRUE)

DimPlot(sc22_mast, label = TRUE) + NoLegend()
DimPlot(sc22_mast, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22_mast, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
  ggtitle('Sample ID') 
DimPlot(sc22_mast, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase')
DimPlot(sc22_mast, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22_mast, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22_mast, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_mast, features=c("nCount_RNA","nFeature_RNA"))
## => based on visualization, data should be integrated to account for the effect of disease group and individual horse

# INTEGRATION AND NORMALIZATION
sc22_mast.list <- SplitObject(sc22_mast, split.by = "disease_state")
sc22_mast.list <- lapply(X = sc22_mast.list, FUN = SCTransform)
features_mast <- SelectIntegrationFeatures(object.list = sc22_mast.list, nfeatures = 3000) 
sc22_mast.list <- PrepSCTIntegration(object.list = sc22_mast.list, anchor.features = features_mast)

immune.anchors_mast <- FindIntegrationAnchors(object.list = sc22_mast.list, normalization.method = "SCT",
                                            anchor.features = features_mast)
sc22_mast <- IntegrateData(anchorset = immune.anchors_mast, normalization.method = "SCT")
sc22_mast <- RunPCA(sc22_mast, verbose = FALSE)
sc22_mast <- RunUMAP(sc22_mast, reduction = "pca", dims = 1:30)


# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(sc22_mast, ndims = 30)

# Determine percent of variation associated with each PC
pct_mast <- sc22_mast[["pca"]]@stdev / sum(sc22_mast[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_mast <- cumsum(pct_mast)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_mast  <- which(cumu_mast  > 90 & pct_mast  < 5)[1]
co1_mast

# Determine the difference between variation of PC and subsequent PC
co2_mast  <- sort(which((pct_mast [1:length(pct_mast) - 1] - pct_mast[2:length(pct_mast)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_mast

# Minimum of the two calculation
pcs_mast  <- min(co1_mast, co2_mast)
pcs_mast

# Create a dataframe with values
plot_mast <- data.frame(pct = pct_mast, 
                      cumu = cumu_mast, 
                      rank = 1:length(pct_mast))

# Elbow plot to visualize 
ggplot(plot_mast, aes(cumu_mast, pct_mast, label = rank, color = rank > pcs_mast)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_mast[pct_mast > 5]), color = "grey") +
  theme_bw()
print(x = sc22_mast[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
## we pick 7 PCs

# CLUSTERING 
sc22_mast <- FindNeighbors(sc22_mast, reduction = "pca", dims = 1:7)
sc22_mast <- FindClusters(sc22_mast, resolution = seq(0.1, 1, by=0.1))
head(sc22_mast@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_mast@meta.data[,grep("integrated_snn_res", colnames(sc22_mast@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_mast, group.by=grep("integrated_snn_res",colnames(sc22_mast@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_mast, group.by=grep("integrated_snn_res",colnames(sc22_mast@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_mast, reduction = "umap", group.by = "integrated_snn_res.0.1", label=T)
# We choose resolution 0.1 (1 cluster)

# Visualization of the dataset at clustering resolution = 0.1
p1_mast <-DimPlot(sc22_mast, reduction = "umap", group.by = "integrated_snn_res.0.1", label=T) + ggtitle('snn_res.0.4')
p2_mast <-DimPlot(sc22_mast, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_mast + p2_mast

DimPlot(sc22_mast, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(sc22_mast, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_mast, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_mast, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_mast, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_mast, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_mast, features=c("nCount_RNA","nFeature_RNA"))

# Gene expression patterns
sc22_mast <- AddModuleScore(sc22_mast,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
sc22_mast <- AddModuleScore(sc22_mast,features = list(eq.s.genes), name = "feat_eq.s.genes")
VlnPlot(sc22_mast, "feat_eq.g2m.genes1") +  FeaturePlot(sc22_mast, "feat_eq.g2m.genes1") 

# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_mastcell_info<-data.frame(sc22_mast@meta.data$CellSubtype,sc22_mast@meta.data$sample_id, sc22_mast@meta.data$disease_state)
colnames(sc_mastcell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_mastcell_info)
# run Propeller
x<-propeller(clusters = sc_mastcell_info$cellType, sample = sc_mastcell_info$sampleId, 
             group = sc_mastcell_info$diseaseGroup)
write.table(x, file = "cellProp_mast.txt",quote = FALSE, sep = "\t")
 

# SAVE THE FINAL OBJECT
saveRDS(sc22_mast, "sc22_mast.rds")
