# INDEPENDENT ANALYSIS OF NEUTROPHILS

# After subsetting neutrophils from the dataset,
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

# Subset the neutrophils from the complete data set 
sc22_ne <- subset(sc22, idents = "Neutrophils") 
DimPlot(sc22_ne) 

# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_ne@meta.data[,grep("integrated_snn_res", colnames(sc22_ne@meta.data))] <- NULL

# NORMALIZATION
sc22_ne <- SCTransform(sc22_ne, method = "glmGamPoi", variable.features.n = 3000)

vf_top10_ne <- head(VariableFeatures(sc22_ne), 10)
vf_plot_ne <- VariableFeaturePlot(sc22_ne)
LabelPoints(plot = vf_plot_ne,
            points = vf_top10_ne, repel = TRUE, xnudge = 0, ynudge = 0)

#Perform dimensionality reduction by PCA and UMAP embedding
sc22_ne <- RunPCA(sc22_ne, verbose = FALSE)
sc22_ne  <- RunUMAP(sc22_ne, dims = 1:30, verbose = FALSE)
sc22_ne  <- FindNeighbors(sc22_ne, dims = 1:30, verbose = FALSE)
sc22_ne  <- FindClusters(sc22_ne, verbose = FALSE)
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_ne <- CellCycleScoring(sc22_ne, 
                           s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                           set.ident = TRUE)

DimPlot(sc22_ne, label = TRUE) + NoLegend()
DimPlot(sc22_ne, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22_ne, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
  ggtitle('Sample ID') #
DimPlot(sc22_ne, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase')
DimPlot(sc22_ne, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22_ne, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22_ne, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_ne, features=c("nCount_RNA","nFeature_RNA"))

# INTEGRATION AND NORMALIZATION
sc22_ne.list <- SplitObject(sc22_ne, split.by = "disease_state")
sc22_ne.list <- lapply(X = sc22_ne.list, FUN = SCTransform)
features_ne <- SelectIntegrationFeatures(object.list = sc22_ne.list, nfeatures = 3000) 
sc22_ne.list <- PrepSCTIntegration(object.list = sc22_ne.list, anchor.features = features_ne)

immune.anchors_ne <- FindIntegrationAnchors(object.list = sc22_ne.list, normalization.method = "SCT",
                                           anchor.features = features_ne)
sc22_ne <- IntegrateData(anchorset = immune.anchors_ne, normalization.method = "SCT")
sc22_ne <- RunPCA(sc22_ne, verbose = FALSE)
sc22_ne <- RunUMAP(sc22_ne, reduction = "pca", dims = 1:30)

# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(sc22_ne, ndims = 30)

# Determine percent of variation associated with each PC
pct_ne <- sc22_ne[["pca"]]@stdev / sum(sc22_ne[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_ne <- cumsum(pct_ne)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_ne  <- which(cumu_ne  > 90 & pct_ne  < 5)[1]
co1_ne  

# Determine the difference between variation of PC and subsequent PC
co2_ne  <- sort(which((pct_ne [1:length(pct_ne ) - 1] - pct_ne[2:length(pct_ne )]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_ne  

# Minimum of the two calculation
pcs_ne  <- min(co1_ne, co2_ne)
pcs_ne  

# Create a dataframe with values
plot_ne <- data.frame(pct = pct_ne, 
                     cumu = cumu_ne, 
                     rank = 1:length(pct_ne))

# Elbow plot to visualize 
ggplot(plot_ne, aes(cumu_ne, pct_ne, label = rank, color = rank > pcs_ne)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_ne[pct_ne > 5]), color = "grey") +
  theme_bw()
print(x = sc22_ne[["pca"]], 
      dims = 1:20, 
      nfeatures = 5)
## we pick 11 PCs

# CLUSTERING 
sc22_ne <- FindNeighbors(sc22_ne, reduction = "pca", dims = 1:11)
sc22_ne <- FindClusters(sc22_ne, resolution = seq(0.05, 0.5, by=0.05))
head(sc22_ne@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_ne@meta.data[,grep("integrated_snn_res", colnames(sc22_ne@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_ne, group.by=grep("integrated_snn_res",colnames(sc22_ne@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_ne, group.by=grep("integrated_snn_res",colnames(sc22_ne@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_ne, group.by=grep("integrated_snn_res",colnames(sc22_ne@meta.data),value = TRUE)[9:12], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_ne, reduction = "umap", group.by = "integrated_snn_res.0.15", label=T)
# We choose resolution 0.15 (3 clusters)

# Visualization of the dataset at clustering resolution = 0.15
p1_ne <-DimPlot(sc22_ne, reduction = "umap", group.by = "integrated_snn_res.0.15", label=T) + ggtitle('snn_res.0.15')
p2_ne <-DimPlot(sc22_ne, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_ne + p2_ne

DimPlot(sc22_ne, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(sc22_ne, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_ne, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_ne, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_ne, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_ne, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_ne, features=c("nCount_RNA","nFeature_RNA"))


# CELL CLUSTER ANNOTATION
# Rename the neutrophils clusters
Idents(sc22_ne) <- "integrated_snn_res.0.15"
sc22_ne <- RenameIdents(sc22_ne, "0"="Neu 0", "1"="Neu 1", "2"="Neu 2")
sc22_ne$CellSubtype <- Idents(object = sc22_ne)

# Create vectors of marker genes (features) for cell subtypes 
feat_infl <- c("IL1A","IL1B","CXCL2","CXCL8") # pro-inflammatory genes
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc22_ne <- AddModuleScore(sc22_ne,features = feat_infl, name = "feat_infl")


# MARKER GENES FOR THE CELL CLUSTERS
# Find the markers for the neutrophil clusters 
Idents(sc22_ne) <- "CellSubtype" 
markers_all_ne <- FindAllMarkers(sc22_ne, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_ne <- subset(markers_all_ne, markers_all_ne$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_ne, "markers_all_ne.rds")
# Supplementary table 4
write.csv(markers_all_ne, 'markers_all_ne.csv')

# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_necell_info<-data.frame(sc22_ne@meta.data$CellSubtype,sc22_ne@meta.data$sample_id, sc22_ne@meta.data$disease_state)
colnames(sc_necell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_necell_info)
# run Propeller
x<-propeller(clusters = sc_necell_info$cellType, sample = sc_necell_info$sampleId, 
             group = sc_necell_info$diseaseGroup)
write.table(x, file = "cellProp_ne.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22_ne, "sc22_ne.rds")

