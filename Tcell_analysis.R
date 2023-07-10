# INDEPENDENT ANALYSIS OF T CELLS

# After subsetting T cells from the dataset,
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
library(reshape2)
library(tibble)
library(speckle)

# Retrieve the Seurat object
sc22 <- readRDS("sc22_seed.rds")
# Set the grouping ident to the 6 major cell types identified
Idents(sc22) <- "MajCellType"

# Subset the T cell group from the complete data set 
sc22_T <- subset(sc22, idents = "T cells") 
DimPlot(sc22_T) 

# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_T@meta.data[,grep("integrated_snn_res", colnames(sc22_T@meta.data))] <- NULL

# NORMALIZATION
sc22_T <- SCTransform(sc22_T, method = "glmGamPoi", variable.features.n = 3000)

vf_top10_T <- head(VariableFeatures(sc22_T), 10)
vf_plot_T <- VariableFeaturePlot(sc22_T)
LabelPoints(plot = vf_plot_T,
            points = vf_top10_T, repel = TRUE, xnudge = 0, ynudge = 0)

#Perform dimensionality reduction by PCA and UMAP embedding
sc22_T <- RunPCA(sc22_T, verbose = FALSE)
sc22_T <- RunUMAP(sc22_T, dims = 1:30, verbose = FALSE)
sc22_T <- FindNeighbors(sc22_T, dims = 1:30, verbose = FALSE)
sc22_T <- FindClusters(sc22_T, verbose = FALSE)
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_T <- CellCycleScoring(sc22_T, 
                              s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                              set.ident = TRUE)

DimPlot(sc22_T, label = TRUE) + NoLegend()
DimPlot(sc22_T, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22_T, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
  ggtitle('Sample ID') #
DimPlot(sc22_T, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase')
DimPlot(sc22_T, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22_T, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22_T, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_T, features=c("nCount_RNA","nFeature_RNA"))

# INTEGRATION AND NORMALIZATION
sc22_T.list <- SplitObject(sc22_T, split.by = "disease_state")
sc22_T.list <- lapply(X = sc22_T.list, FUN = SCTransform)
features_T <- SelectIntegrationFeatures(object.list = sc22_T.list, nfeatures = 3000) 
sc22_T.list <- PrepSCTIntegration(object.list = sc22_T.list, anchor.features = features_T)

immune.anchors_T <- FindIntegrationAnchors(object.list = sc22_T.list, normalization.method = "SCT",
                                              anchor.features = features_T)
sc22_T <- IntegrateData(anchorset = immune.anchors_T, normalization.method = "SCT")
sc22_T <- RunPCA(sc22_T, verbose = FALSE)
sc22_T <- RunUMAP(sc22_T, reduction = "pca", dims = 1:30)

# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(sc22_T, ndims = 30)

# Determine percent of variation associated with each PC
pct_T <- sc22_T[["pca"]]@stdev / sum(sc22_T[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_T <- cumsum(pct_T)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_T <- which(cumu_T > 90 & pct_T < 5)[1]
co1_T 

# Determine the difference between variation of PC and subsequent PC
co2_T <- sort(which((pct_T[1:length(pct_T) - 1] - pct_T[2:length(pct_T)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_T 

# Minimum of the two calculation
pcs_T <- min(co1_T, co2_T)
pcs_T 

# Create a dataframe with values
plot_T <- data.frame(pct = pct_T, 
                        cumu = cumu_T, 
                        rank = 1:length(pct_T))

# Elbow plot to visualize 
ggplot(plot_T, aes(cumu_T, pct_T, label = rank, color = rank > pcs_T)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_T[pct_T > 5]), color = "grey") +
  theme_bw()
print(x = sc22_T[["pca"]], 
      dims = 1:25, 
      nfeatures = 5)
## we pick 18 PCs

# CLUSTERING 
sc22_T <- FindNeighbors(sc22_T, reduction = "pca", dims = 1:18)
sc22_T <- FindClusters(sc22_T, resolution = seq(0.1, 1, by=0.1))
head(sc22_T@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_T@meta.data[,grep("integrated_snn_res", colnames(sc22_T@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_T, group.by=grep("integrated_snn_res",colnames(sc22_T@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_T, group.by=grep("integrated_snn_res",colnames(sc22_T@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_T, reduction = "umap", group.by = "integrated_snn_res.0.3", label=T)
# We choose resolution 0.3 (7 clusters)

# Visualization of the dataset at clustering resolution = 0.3
p1_T <-DimPlot(sc22_T, reduction = "umap", group.by = "integrated_snn_res.0.3", label=T) + ggtitle('snn_res.0.3')
p2_T <-DimPlot(sc22_T, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_T + p2_T

DimPlot(sc22_T, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(sc22_T, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_T, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_T, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_T, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_T, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_T, features=c("nCount_RNA","nFeature_RNA"))


# CELL CLUSTER ANNOTATION
# Rename the T cells subtypes 
Idents(sc22_T) <- "integrated_snn_res.0.3"
sc22_T <- RenameIdents(sc22_T, "0"="T0", "1"="T1", "2"="T2", "3"="T3", 
                          "4"="T4", "5"="T5", "6"="T6")
sc22_T$CellSubtype <- Idents(object = sc22_T)

# Create vectors of marker genes (features) for cell subtypes 
eq.s.genes <- readRDS("eq.s.genes.rds") # S phase cell cycle markers
eq.g2m.genes <- readRDS("eq.g2m.genes.rds") # G2M phase cell cycle markers
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc22_T <- AddModuleScore(sc22_T,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
sc22_T <- AddModuleScore(sc22_T,features = list(eq.s.genes), name = "feat_eq.s.genes")


# MARKER GENES FOR THE CELL CLUSTERS
# Find the markers for the T cell subtypes 
Idents(sc22_T) <- "CellSubtype" 
markers_all_T <- FindAllMarkers(sc22_T, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_T <- subset(markers_all_T, markers_all_T$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_T, "markers_all_T.rds")
# Supplementary Table 6
write.csv(markers_all_T, 'markers_all_T.csv')


# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_Tcell_info<-data.frame(sc22_T@meta.data$CellSubtype,sc22_T@meta.data$sample_id, sc22_T@meta.data$disease_state)
colnames(sc_Tcell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_Tcell_info)
# run Propeller
x<-propeller(clusters = sc_Tcell_info$cellType, sample = sc_Tcell_info$sampleId, 
             group = sc_Tcell_info$diseaseGroup)
write.table(x, file = "cellProp_T.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22_T, "sc22_T_seed.rds")

