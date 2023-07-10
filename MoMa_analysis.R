# INDEPENDENT ANALYSIS OF MONOCYTES-MACROPHAGES

# After subsetting monocytes-macrophages from the dataset,
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
library(dittoSeq)
library(data.table)
library(ggpubr)
library(tibble)
library(speckle)

# Retrieve the Seurat object
sc22 <- readRDS("sc22_seed.rds")
# Set the grouping ident to the 6 major cell types identified
Idents(sc22) <- "MajCellType"

# Subset the Mo/Ma cell group from the complete data set 
sc22_moma <- subset(sc22, idents = "Mo/Ma") 
DimPlot(sc22_moma) 
sc22_moma

# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_moma@meta.data[,grep("integrated_snn_res", colnames(sc22_moma@meta.data))] <- NULL

# NORMALIZATION
sc22_moma <- SCTransform(sc22_moma, method = "glmGamPoi", variable.features.n = 3000)

vf_top10_moma <- head(VariableFeatures(sc22_moma), 10)
vf_plot_moma <- VariableFeaturePlot(sc22_moma)
LabelPoints(plot = vf_plot_moma,
            points = vf_top10_moma, repel = TRUE, xnudge = 0, ynudge = 0)

#Perform dimensionality reduction by PCA and UMAP embedding
sc22_moma <- RunPCA(sc22_moma, verbose = FALSE)
sc22_moma <- RunUMAP(sc22_moma, dims = 1:30, verbose = FALSE)
sc22_moma <- FindNeighbors(sc22_moma, dims = 1:30, verbose = FALSE)
sc22_moma <- FindClusters(sc22_moma, verbose = FALSE)
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_moma <- CellCycleScoring(sc22_moma, 
                              s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                              set.ident = TRUE)

DimPlot(sc22_moma, label = TRUE) + NoLegend()
DimPlot(sc22_moma, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22_moma, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
  ggtitle('Sample ID') 
DimPlot(sc22_moma, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase')
DimPlot(sc22_moma, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22_moma, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22_moma, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_moma, features=c("nCount_RNA","nFeature_RNA"))

# INTEGRATION AND NORMALIZATION
sc22_moma.list <- SplitObject(sc22_moma, split.by = "disease_state")
sc22_moma.list <- lapply(X = sc22_moma.list, FUN = SCTransform)
features_moma <- SelectIntegrationFeatures(object.list = sc22_moma.list, nfeatures = 3000) 
sc22_moma.list <- PrepSCTIntegration(object.list = sc22_moma.list, anchor.features = features_moma)

immune.anchors_moma <- FindIntegrationAnchors(object.list = sc22_moma.list, normalization.method = "SCT",
                                              anchor.features = features_moma)
sc22_moma <- IntegrateData(anchorset = immune.anchors_moma, normalization.method = "SCT")
sc22_moma <- RunPCA(sc22_moma, verbose = FALSE)
sc22_moma <- RunUMAP(sc22_moma, reduction = "pca", dims = 1:30)


# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(sc22_moma, ndims = 30)

# Determine percent of variation associated with each PC
pct_moma <- sc22_moma[["pca"]]@stdev / sum(sc22_moma[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_moma <- cumsum(pct_moma)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_moma <- which(cumu_moma > 90 & pct_moma < 5)[1]
co1_moma 

# Determine the difference between variation of PC and subsequent PC
co2_moma <- sort(which((pct_moma[1:length(pct_moma) - 1] - pct_moma[2:length(pct_moma)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_moma 

# Minimum of the two calculation
pcs_moma <- min(co1_moma, co2_moma)
pcs_moma 

# Create a dataframe with values
plot_moma <- data.frame(pct = pct_moma, 
                        cumu = cumu_moma, 
                        rank = 1:length(pct_moma))

# Elbow plot to visualize 
ggplot(plot_moma, aes(cumu_moma, pct_moma, label = rank, color = rank > pcs_moma)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_moma[pct_moma > 5]), color = "grey") +
  theme_bw()
print(x = sc22_moma[["pca"]], 
      dims = 1:25, 
      nfeatures = 5)
## we pick 18 PCs

# CLUSTERING 
sc22_moma <- FindNeighbors(sc22_moma, reduction = "pca", dims = 1:18)
sc22_moma <- FindClusters(sc22_moma, resolution = seq(0.05, 0.4, by=0.05))
head(sc22_moma@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_moma@meta.data[,grep("integrated_snn_res", colnames(sc22_moma@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_moma, group.by=grep("integrated_snn_res",colnames(sc22_moma@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_moma, group.by=grep("integrated_snn_res",colnames(sc22_moma@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_moma, reduction = "umap", group.by = "integrated_snn_res.0.2", label=T)
# We choose resolution 0.2 (6 clusters)

# Visualization of the dataset at clustering resolution = 0.2
p1_moma <-DimPlot(sc22_moma, reduction = "umap", group.by = "integrated_snn_res.0.2", label=T) + ggtitle('snn_res.0.2')
p2_moma <-DimPlot(sc22_moma, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_moma + p2_moma

DimPlot(sc22_moma, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(sc22_moma, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_moma, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_moma, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_moma, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_moma, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_moma, features=c("nCount_RNA","nFeature_RNA"))


# CELL CLUSTER ANNOTATION
# Rename the Mo/Ma subtypes 
Idents(sc22_moma) <- "integrated_snn_res.0.2"
sc22_moma <- RenameIdents(sc22_moma, "0"="Mo/Ma 0", "1"="Mo/Ma 1", "2"="Mo/Ma 2", "3"="Mo/Ma 3", 
                          "4"="Mo/Ma 4", "5"="Mo/Ma 5")
sc22_moma$CellSubtype <- Idents(object = sc22_moma)
table(Idents(sc22_moma))

# Create vectors of marker genes (features) for cell subtypes 
feat_AM <- c("MARCO", "APOE", "MSR1", "CD163" ) # alveolar macrophages markers
feat_HLA <- c("DQB", "DQB.1", "DQA", "DQA.1","DRA", "DRB") # MHCII-associated genes
feat_ISG <- c("OASL", "IFI6", "IFI44", "IRF7") # interferon stimulated genes
feat_T <- c("CD2", "CD3D", "CD3E", "CD3G") # T cell markers
eq.s.genes <- readRDS("eq.s.genes.rds") # S phase cell cycle markers
eq.g2m.genes <- readRDS("eq.g2m.genes.rds") # G2M phase cell cycle markers
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc22_moma <- AddModuleScore(sc22_moma,features = feat_AM, name = "feat_AM")
sc22_moma <- AddModuleScore(sc22_moma,features = feat_HLA, name = "feat_HLA")
sc22_moma <- AddModuleScore(sc22_moma,features = feat_ISG, name = "feat_ISG")
sc22_moma <- AddModuleScore(sc22_moma,features = feat_T, name = "feat_T")
sc22_moma <- AddModuleScore(sc22_moma,features = list(eq.s.genes), name = "feat_eq.s.genes")
sc22_moma <- AddModuleScore(sc22_moma,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")


# MARKER GENES FOR THE CELL CLUSTERS
# Find the markers for the Mo/Ma cell subtypes 
Idents(sc22_moma) <- "CellSubtype" 
markers_all_moma <- FindAllMarkers(sc22_moma, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_moma <- subset(markers_all_moma, markers_all_moma$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_moma, "markers_all_moma.rds")
# Supplementary table 7
write.csv(markers_all_moma, 'markers_all_moma.csv')


# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_momacell_info<-data.frame(sc22_moma@meta.data$CellSubtype,sc22_moma@meta.data$sample_id, sc22_moma@meta.data$disease_state)
colnames(sc_momacell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_momacell_info)
# run Propeller
x<-propeller(clusters = sc_momacell_info$cellType, sample = sc_momacell_info$sampleId, 
             group = sc_momacell_info$diseaseGroup)
write.table(x, file = "cellProp_moma.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22_moma, "sc22_moma.rds")
