# INDEPENDENT ANALYSIS OF DENDRITIC CELLS

# After subsetting dendritic cells from the dataset,
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

# Subset the dendritic cells from the complete data set 
sc22_DC <- subset(sc22, idents = "Dendritic cells") 
DimPlot(sc22_DC) 

# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_DC@meta.data[,grep("integrated_snn_res", colnames(sc22_DC@meta.data))] <- NULL

# NORMALIZATION
sc22_DC <- SCTransform(sc22_DC, method = "glmGamPoi", variable.features.n = 3000)

vf_top10_DC <- head(VariableFeatures(sc22_DC), 10)
vf_plot_DC <- VariableFeaturePlot(sc22_DC)
LabelPoints(plot = vf_plot_DC,
            points = vf_top10_DC, repel = TRUE, xnudge = 0, ynudge = 0)

# Perform dimensionality reduction by PCA and UMAP embedding
sc22_DC <- RunPCA(sc22_DC, verbose = FALSE)
sc22_DC <- RunUMAP(sc22_DC, dims = 1:30, verbose = FALSE)
sc22_DC<- FindNeighbors(sc22_DC, dims = 1:30, verbose = FALSE)
sc22_DC <- FindClusters(sc22_DC, verbose = FALSE)
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_DC<- CellCycleScoring(sc22_DC, 
                              s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                              set.ident = TRUE)

DimPlot(sc22_DC, label = TRUE) + NoLegend()
DimPlot(sc22_DC, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22_DC, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
  ggtitle('Sample ID') 
DimPlot(sc22_DC, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase')
DimPlot(sc22_DC, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22_DC, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22_DC, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_DC, features=c("nCount_RNA","nFeature_RNA"))
## => based on visualization, data should be integrated to account for the effect of disease group and individual horse

# INTEGRATION AND NORMALIZATION
sc22_DC.list <- SplitObject(sc22_DC, split.by = "disease_state")
sc22_DC.list <- lapply(X = sc22_DC.list, FUN = SCTransform)
features_DC <- SelectIntegrationFeatures(object.list = sc22_DC.list, nfeatures = 3000) 
sc22_DC.list <- PrepSCTIntegration(object.list = sc22_DC.list, anchor.features = features_DC)

immune.anchors_DC <- FindIntegrationAnchors(object.list = sc22_DC.list, normalization.method = "SCT",
                                              anchor.features = features_DC)
sc22_DC<- IntegrateData(anchorset = immune.anchors_DC, normalization.method = "SCT")
sc22_DC<- RunPCA(sc22_DC, verbose = FALSE)
sc22_DC<- RunUMAP(sc22_DC, reduction = "pca", dims = 1:30)

# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(sc22_DC, ndims = 30)

# Determine percent of variation associated with each PC
pct_DC <- sc22_DC[["pca"]]@stdev / sum(sc22_DC[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_DC <- cumsum(pct_DC)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_DC  <- which(cumu_DC  > 90 & pct_DC  < 5)[1]
co1_DC  

# Determine the difference between variation of PC and subsequent PC
co2_DC  <- sort(which((pct_DC [1:length(pct_DC) - 1] - pct_DC[2:length(pct_DC)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_DC  

# Minimum of the two calculation
pcs_DC  <- min(co1_DC, co2_DC)
pcs_DC  

# Create a dataframe with values
plot_DC <- data.frame(pct = pct_DC, 
                        cumu = cumu_DC, 
                        rank = 1:length(pct_DC))

# Elbow plot to visualize 
ggplot(plot_DC, aes(cumu_DC, pct_DC, label = rank, color = rank > pcs_DC)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_DC[pct_DC > 5]), color = "grey") +
  theme_bw()
print(x = sc22_DC[["pca"]], 
      dims = 1:15, 
      nfeatures = 5)
## we pick 10 PCs

# CLUSTERING 
sc22_DC<- FindNeighbors(sc22_DC, reduction = "pca", dims = 1:10)
sc22_DC<- FindClusters(sc22_DC, resolution = seq(0.1, 1, by=0.1))
head(sc22_DC@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_DC@meta.data[,grep("integrated_snn_res", colnames(sc22_DC@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_DC, group.by=grep("integrated_snn_res",colnames(sc22_DC@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_DC, group.by=grep("integrated_snn_res",colnames(sc22_DC@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_DC, reduction = "umap", group.by = "integrated_snn_res.0.1", label=T)
# We choose resolution 0.1 (4 clusters)

# Visualization of the dataset at the chosen clustering resolution
p1_DC <-DimPlot(sc22_DC, reduction = "umap", group.by = "integrated_snn_res.0.1", label=T) + ggtitle('snn_res.0.1')
p2_DC <-DimPlot(sc22_DC, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_DC + p2_DC

DimPlot(sc22_DC, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(sc22_DC, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_DC, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_DC, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_DC, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_DC, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_DC, features=c("nCount_RNA","nFeature_RNA"))

# CELL CLUSTER ANNOTATION
# Rename the clusters
Idents(sc22_DC) <- "integrated_snn_res.0.1"
sc22_DC<- RenameIdents(sc22_DC, "0"="DC 0", "1"="DC 1", "2"="DC 2", "3" = "DC 3")
sc22_DC$CellSubtype <- Idents(object = sc22_DC)

# Create vectors of marker genes (features) for cell subtypes 
feat_cDC2 <- c("LOC100072936","LOC100072933","CD1C","FCER1A")
feat_HLA <- c("DQB", "DQB.1", "DQA", "DQA.1","DRA", "DRB") # MHCII-associated genes
feat_act <- c("CCR7","LAMP3","IDO1","CD83") # genes associated with DC activation
feat_AM <- c("MARCO","APOE","CD163","MSR1") # markers for alveolar macrophages
feat_cDC1 <- c("XCR1","CLEC9A","CADM1")
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc22_DC <- AddModuleScore(sc22_DC,features = list(feat_cDC2), name = "feat_cDC2") 
sc22_DC <- AddModuleScore(sc22_DC, features = feat_HLA, name = "feat_HLA")
sc22_DC <- AddModuleScore(sc22_DC,features = list(feat_act), name = "feat_act") 
sc22_DC <- AddModuleScore(sc22_DC,features = list(feat_AM), name = "feat_AM") 
sc22_DC <- AddModuleScore(sc22_DC,features = list(feat_cDC1), name = "feat_cDC1")
                           

# MARKER GENES FOR THE CELL CLUSTERS
# Find the markers for the DC cell clusters 
Idents(sc22_DC) <- "CellSubtype" 
markers_all_DC <- FindAllMarkers(sc22_DC, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_DC <- subset(markers_all_DC, markers_all_DC$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_DC, "markers_all_DC.rds")
# Supplementary table 8
write.csv(markers_all_DC, 'markers_all_DC.csv')


# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_dccell_info<-data.frame(sc22_dc@meta.data$CellSubtype,sc22_dc@meta.data$sample_id, sc22_dc@meta.data$disease_state)
colnames(sc_dccell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_dccell_info)
# run Propeller
x<-propeller(clusters = sc_dccell_info$cellType, sample = sc_dccell_info$sampleId, 
             group = sc_dccell_info$diseaseGroup)
write.table(x, file = "cellProp_DC.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22_DC, "sc22_DC.rds")
