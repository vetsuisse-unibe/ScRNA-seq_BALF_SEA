# INDEPENDENT ANALYSIS OF B CELLS
# After subsetting B cells from the dataset,
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

# Retrieve the Seurat object with all cells
sc22 <- readRDS("sc22_seed.rds")
# Set the grouping ident to the 6 major cell types identified
Idents(sc22) <- "MajCellType"

# Subset the B cell group from the complete data set 
sc22_B <- subset(sc22, idents = "B/Plasma cells") 
DimPlot(sc22_B) 
# Remove the existing metadata column with cluster resolution to avoid confusion with later subclustering
sc22_B@meta.data[,grep("SCT_snn_res", colnames(sc22_B@meta.data))] <- NULL

# NORMALIZATION
sc22_B <- SCTransform(sc22_B, method = "glmGamPoi", variable.features.n = 3000)
DefaultAssay(sc22_B)

vf_top10_B <- head(VariableFeatures(sc22_B), 10)
vf_plot_B <- VariableFeaturePlot(sc22_B)
LabelPoints(plot = vf_plot_B,
            points = vf_top10_B, repel = TRUE, xnudge = 0, ynudge = 0)

# Cell cycle scoring
eq.s.genes <- readRDS("eq.s.genes.rds")
eq.g2m.genes <- readRDS("eq.g2m.genes.rds")
sc22_B <- CellCycleScoring(sc22_B, 
                           s.features = eq.s.genes, g2m.features = eq.g2m.genes, 
                           set.ident = TRUE)

#Perform dimensionality reduction by PCA and UMAP embedding
sc22_B <- RunPCA(sc22_B, verbose = FALSE)

# INTEGRATION 
# Cannot be performed due to the low number of cells (gives an error) 

# PCA
# Choice of the number of dimensions
ElbowPlot(sc22_B, ndims = 30)

# Determine percent of variation associated with each PC
pct_B <- sc22_B[["pca"]]@stdev / sum(sc22_B[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu_B <- cumsum(pct_B)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1_B <- which(cumu_B > 90 & pct_B < 5)[1]
co1_B 

# Determine the difference between variation of PC and subsequent PC
co2_B <- sort(which((pct_B[1:length(pct_B) - 1] - pct_B[2:length(pct_B)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2_B

# Minimum of the two calculation
pcs_B <- min(co1_B, co2_B)
pcs_B 

# Create a dataframe with values
plot_B <- data.frame(pct = pct_B, 
                     cumu = cumu_B, 
                     rank = 1:length(pct_B))

# Elbow plot to visualize 
ggplot(plot_B, aes(cumu_B, pct_B, label = rank, color = rank > pcs_B)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_B[pct_B > 5]), color = "grey") +
  theme_bw()
print(x = sc22_B[["pca"]], 
      dims = 1:15, 
      nfeatures = 5)
## we pick 11 PCs
sc22_B <- RunUMAP(sc22_B, dims = 1:11, verbose = FALSE)


# CELL CLUSTERING
sc22_B <- FindNeighbors(sc22_B, reduction = "pca", dims = 1:11)
sc22_B <- FindClusters(sc22_B, resolution = seq(0.1, 1.2, by=0.1))
head(sc22_B@meta.data)

# CHOOSE THE RESOLUTION 
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(sc22_B@meta.data[,grep("SCT_snn_res", colnames(sc22_B@meta.data))],
         prefix = "SCT_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = sc22_B, group.by=grep("SCT_snn_res",colnames(sc22_B@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc22_B, group.by=grep("SCT_snn_res",colnames(sc22_B@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc22_B, reduction = "umap", group.by = "SCT_snn_res.0.3", label=T)
# We choose resolution 0.3 (3 clusters)

# Visualization of the dataset at clustering resolution = 0.3
p1_B <-DimPlot(sc22_B, reduction = "umap", group.by = "SCT_snn_res.0.3", label=T) + ggtitle('snn_res.0.3')
p2_B <-DimPlot(sc22_B, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1_B + p2_B

DimPlot(sc22_B, group.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_B, group.by='breed', split.by='breed', reduction='umap') + ggtitle('Breed') 
DimPlot(sc22_B, group.by='sex', reduction='umap') + ggtitle('Sex')
DimPlot(sc22_B, group.by='sex', split.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(sc22_B, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22_B, features=c("nCount_RNA","nFeature_RNA"))

# CELL CLUSTER ANNOTATION
# Rename the B cells subtypes 
Idents(sc22_B) <- "SCT_snn_res.0.3"
sc22_B <- RenameIdents(sc22_B, "0"="B0", "1"="B1", "2"="B2")
sc22_B$CellSubtype <- Idents(object = sc22_B)

# Create vectors of marker genes (features) for cell subtypes 
feat_HLA <- c("DQB", "DQB.1", "DQA", "DQA.1","DRA", "DRB") # MHCII-associated genes
feat_ab <- c("TXNDC5", "HSP90B1", "FAM46C") #genes associated with antibody production (FAM46C=TENT5C)
feat_igm <- c("JCHAIN", "MZB1") #genes associated with IgM production
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc22_B <- AddModuleScore(sc22_B,features = list(feat_HLA), name = "feat_HLA") 
sc22_B <- AddModuleScore(sc22_B,features = list(feat_ab), name = "feat_ab") 
sc22_B <- AddModuleScore(sc22_B,features = list(feat_igm), name = "feat_igm") 

# Find the markers for the B cell subtypes 
sc22_B <- SetIdent(sc22_B, value = sc22_B$CellSubtype)
markers_all_B <- FindAllMarkers(sc22_B, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_B <- subset(markers_all_B, markers_all_B$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_B, "markers_all_B.rds")
# Supplementary Table 5 
write.csv(markers_all_B, 'markers_B.csv')


# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_Bcell_info<-data.frame(sc22_B@meta.data$CellSubtype,sc22_B@meta.data$sample_id, sc22_B@meta.data$disease_state)
colnames(sc_Bcell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_Bcell_info)
# run Propeller
propeller(clusters = sc_Bcell_info$cellType, sample = sc_Bcell_info$sampleId, 
          group = sc_Bcell_info$diseaseGroup)
write.table(x, file = "cellProp_B.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22_B, "sc22_B_seed.rds")

