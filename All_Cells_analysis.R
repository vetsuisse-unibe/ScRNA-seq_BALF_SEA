# Analysis of the scRNA-seq dataset obtained from the cryopreserved bronchoalveolar cells of 
# 6 horses with severe equine asthma (SEA) and 5 control horses

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

sea_sc2021 <- readRDS("initial_seurat_object.rds")

# Quality control and filtering 
# Add mitochondrial reads to the Seurat object's meta.data
mito.genes<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB","MT-tRNA-Phe",  
              "MT-s-rRNA","MT-tRNA-Val","MT-l-rRNA","MT-tRNA-Leu","MT-tRNA-Ile","MT-tRNA-Gln","MT-tRNA-Met","MT-tRNA-Trp","MT-tRNA-Ala",
              "MT-tRNA-Cys","MT-tRNA-Tyr","MT-tRNA-Ser","MT-tRNA-Asp","MT-tRNA-Lys","MT-tRNA-Gly","MT-tRNA-Arg","MT-tRNA-His","MT-tRNA-Ser.1",
              "MT-tRNA-Leu.1","MT-tRNA-Glu","MT-tRNA-Thr","MT-tRNA-Pro")#36/37 MT genes because tRNA-Asn is absent from this dataset

sea_sc2021[["percent.mito"]] <- PercentageFeatureSet(sea_sc2021, features=mito.genes)
summary(sea_sc2021$percent.mito)

sea_sc2021<- PercentageFeatureSet(sea_sc2021, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")

# Initial gene metrics 
kable(do.call("cbind", tapply(sea_sc2021$nFeature_RNA,
                              Idents(sea_sc2021),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of Genes/Cell by Sample") %>% kable_styling()

kable(round(do.call("cbind", tapply(sea_sc2021$percent.mito, Idents(sea_sc2021),quantile,probs=seq(0,1,0.05))), digits = 3),
      caption = "5% Quantiles of Percent Mitochondria by Sample") %>% kable_styling()

# Sample wise Quality metrics 
VlnPlot(
  sea_sc2021,
  features = c("nFeature_RNA", "nCount_RNA","percent.mito"),
  ncol = 1, pt.size = 0.3)

RidgePlot(sea_sc2021, features=c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol = 2)
plot(sort(Matrix::rowSums(GetAssayData(sea_sc2021) >= 3), decreasing = TRUE) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito", shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired")) + geom_vline(xintercept = c(1000,12000)) + geom_hline(yintercept = 8)
FeatureScatter(sea_sc2021, feature1 = "nFeature_RNA", feature2 = "percent.mito", shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired")) + geom_vline(xintercept = 700) + geom_hline(yintercept = 8)
FeatureScatter(sea_sc2021, "nCount_RNA", "nFeature_RNA",pt.size = 0.5, shuffle = TRUE,cols=brewer.pal(n = 11, name = "Paired"))  + geom_vline(xintercept = c(1000,12000)) + geom_hline(yintercept = 700)

# Looking for the best thresholds to filter for nFeature_RNA
(VlnPlot(sea_sc2021, features = "nFeature_RNA", pt.size = 0) + 
    geom_hline(yintercept=200, linetype="dashed", color = "red") +
    geom_hline(yintercept=8000, linetype="dashed", color = "red") +
    labs(title= NULL, x=NULL, y="Feature count")) &
  (scale_fill_discrete(limits = c("30", "33", "44","46","51","52","56","68","89","90","91"),
                       labels = c("Horse 30", "Horse 33", "Horse 44","Horse 46","Horse 51","Horse 52","Horse 56","Horse 68","Horse 89","Horse 90","Horse 91")))
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept=200, linetype="dashed", color = "red") +
  geom_hline(yintercept=8000, linetype="dashed", color = "red") + 
  labs(title= NULL, x="RNA count", y="Feature count")

# Figure justifying choice of nFeature_RNA and percent.mito thresholds for cell filtering
Plot_FilterVln <- ((VlnPlot(sea_sc2021, features = "nFeature_RNA", pt.size = 0) + 
                      geom_hline(yintercept=200, linetype="dashed", color = "red") +
                      geom_hline(yintercept=8000, linetype="dashed", color = "red") +
                      labs(title= NULL, x=NULL, y="Feature count")) &
                     (VlnPlot(sea_sc2021, features = "percent.mito", pt.size = 0) + 
                        geom_hline(yintercept=15, linetype="dashed", color = "red") +
                        labs(title= NULL, x=NULL, y="Mitochondrial reads %")) +
                     plot_layout(ncol = 1)) &
  (theme(axis.title = element_text(size=10), axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(), axis.text.y = element_text(size=8))) &
  scale_fill_discrete(limits = c("30", "33", "44","46","51","52","56","68","89","90","91"),
                      labels = c("Horse 30", "Horse 33", "Horse 44","Horse 46","Horse 51","Horse 52","Horse 56","Horse 68","Horse 89","Horse 90","Horse 91"))

Plot_FilterFeat <- ((FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
                       geom_hline(yintercept=200, linetype="dashed", color = "red") +
                       geom_hline(yintercept=8000, linetype="dashed", color = "red") + 
                       labs(title= NULL, x="RNA count", y="Feature count")) +
                      (FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
                         geom_hline(yintercept=15, linetype="dashed", color = "red") 
                       + labs(title= NULL, x="RNA count", y="Mitochondrial reads %")) +
                      plot_layout(ncol = 1)) &
  theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) &
  NoLegend()

(Plot_FilterFeat | Plot_FilterVln) + plot_annotation(tag_levels = 'A')

# QC metrics grouped by disease state
VlnPlot(sea_sc2021, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 3, group.by = "disease_state",
        pt.size = 0.001, combine = T) & theme(legend.position = 'none',
                                              axis.title.x = element_blank())
VlnPlot(sea_sc2021, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 3, group.by = "disease_state",
        pt.size = 0, combine = T) & theme(legend.position = 'none', 
                                          axis.title.x = element_blank())

FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by = "disease_state")
FeatureScatter(sea_sc2021, feature1 = "nFeature_RNA", feature2 = "percent.ribo", group.by = "disease_state")

# Doublet detection 
## Convert Seurat object into a SingleCellExperiment object
sea.sce <- as.SingleCellExperiment(sea_sc2021)
## Find Doublets 
sea.sce <- scDblFinder(sea.sce, samples="sample_id", clusters=NULL)
## Convert sce object back to Seurat object
sea_seurat <- as.Seurat(sea.sce, counts = "counts", data = "logcounts")

head(sea_seurat@meta.data) %>%
  kable(format = "html", col.names = colnames(head(sea_seurat@meta.data))) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

## Count the number of doublets detected
table(sea_seurat$scDblFinder.class)
prop.table(table(sea_seurat$scDblFinder.class))

## Filter the doublets
sea_singlet <- subset(sea_seurat, subset = scDblFinder.class  == "singlet")
dim(sea_seurat)
dim(sea_singlet)

## Apply additional QC filtering using the threshols previously defined 
sea_singlet_flt <- subset(sea_singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 8000
                          & percent.mito < 15 )
dim(sea_singlet_flt)


# Update Metadata with MT and Rb genes information
sea_singlet_flt[["percent.mito"]] <- PercentageFeatureSet(sea_singlet_flt, features=mito.genes)
summary(sea_singlet_flt$percent.mito)

sea_singlet_flt<- PercentageFeatureSet(sea_singlet_flt, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
summary(sea_singlet_flt$percent.ribo)

VlnPlot(sea_singlet_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, group.by = "disease_state",
        pt.size = 0.001, combine = T) & theme(legend.position = 'none',
                                              axis.title.x = element_blank())

VlnPlot(sea_singlet_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, group.by = "disease_state",
        pt.size = 0, combine = T) & theme(legend.position = 'none', 
                                          axis.title.x = element_blank())

FeatureScatter(sea_singlet_flt, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "disease_state")
FeatureScatter(sea_singlet_flt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "disease_state") 

# NORMALIZATION AND QC
sc22 <- sea_singlet_flt 

# Normalization using sctransform
sc22<- SCTransform(sc22, method = "glmGamPoi", variable.features.n = 3000)
DefaultAssay(sc22) 

vf_top10 <- head(VariableFeatures(sc22), 10)
vf_plot <- VariableFeaturePlot(sc22)
LabelPoints(plot = vf_plot,
            points = vf_top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Calculate a cell cycle score for each cell
convertHumanGeneList <- function(x){
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  horse = useEnsembl("ensembl", dataset = "ecaballus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("external_gene_name"), martL = horse, uniqueRows=T)
  humanx <- unique(genes[, 2])
  print(head(humanx))
  return(humanx)
}
eq.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
eq.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

sc22 <- CellCycleScoring(sc22, s.features = eq.s.genes, g2m.features = eq.g2m.genes, set.ident = TRUE)
head(sc22)
table(sc22@meta.data$Phase) %>% kable(caption = "Number of Cells in each Cell Cycle Stage", col.names = c("Stage", "Count"), align = "c") %>% kable_styling()
table(Idents(sc22))
Idents(sc22) <- "orig.ident"

#Perform dimensionality reduction by PCA and UMAP embedding
sc22 <- RunPCA(sc22, verbose = FALSE)
sc22 <- RunUMAP(sc22, dims = 1:30, verbose = FALSE)
sc22 <- FindNeighbors(sc22, dims = 1:30, verbose = FALSE)
sc22 <- FindClusters(sc22, verbose = FALSE)
DimPlot(sc22, label = TRUE) + NoLegend()
DimPlot(sc22, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
DimPlot(sc22, group.by='orig.ident', reduction='umap',
        cols=c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD","#5E4FA2")) +
        ggtitle('Sample ID') 
DimPlot(sc22, group.by='Phase', reduction='umap') + ggtitle('Cell Cycle Phase') 
DimPlot(sc22, group.by='breed', reduction='umap') + ggtitle('Breed')  
DimPlot(sc22, group.by='sex', reduction='umap') + ggtitle('Sex') 
FeaturePlot(sc22, features=c("percent.mito","percent.ribo"))
FeaturePlot(sc22, features=c("nCount_RNA","nFeature_RNA"))
## => based on visualization, data should be integrated to account for the effect of disease group and individual horse

# INTEGRATION AND NORMALIZATION
sc22.list <- SplitObject(sc22, split.by = "disease_state")
sc22.list <- lapply(X = sc22.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = sc22.list, nfeatures = 3000) #default (pilot: 2000)
sc22.list <- PrepSCTIntegration(object.list = sc22.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = sc22.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct

# PCA
# Choice of the number of dimensions for PCA 
ElbowPlot(immune.combined.sct, ndims = 30)

# Determine percent of variation associated with each PC
pct <- immune.combined.sct[["pca"]]@stdev / sum(immune.combined.sct[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2 

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs 

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
print(x = immune.combined.sct[["pca"]], 
      dims = 1:25, 
      nfeatures = 5)
## we pick 16 PCs


# CLUSTERING 
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:16)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = seq(0.1, 1.2, by=0.1))
head(immune.combined.sct@meta.data)
# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)

# Visualize how clusters sub-divide at increasing resolution:
clustree(immune.combined.sct@meta.data[,grep("integrated_snn_res", colnames(immune.combined.sct@meta.data))],
         prefix = "integrated_snn_res.")

# Visualize the UMAP at different resolutions:
(DimPlot(object = immune.combined.sct, group.by=grep("res",colnames(immune.combined.sct@meta.data),value = TRUE)[2:5], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

(DimPlot(object = immune.combined.sct, group.by=grep("res",colnames(immune.combined.sct@meta.data),value = TRUE)[6:9], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

(DimPlot(object = immune.combined.sct, group.by=grep("res",colnames(immune.combined.sct@meta.data),value = TRUE)[10:13], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(immune.combined.sct, reduction = "umap", group.by = "integrated_snn_res.0.6", label=T)
# We choose a resolution of 0.6 (19 clusters)

# Visualization of the dataset at clustering resolution = 0.6
p1<-DimPlot(immune.combined.sct, reduction = "umap", group.by = "integrated_snn_res.0.6", label=T) + ggtitle('SCT_snn_res.0.6')
p2<-DimPlot(immune.combined.sct, group.by='disease_state', reduction='umap') + ggtitle('Disease type')
p1 + p2

DimPlot(immune.combined.sct, group.by='orig.ident', reduction='umap') + ggtitle('Sample origin')
DimPlot(immune.combined.sct, group.by='breed', reduction='umap') + ggtitle('Breed')
DimPlot(immune.combined.sct, group.by='sex', reduction='umap') + ggtitle('Sex')
FeaturePlot(immune.combined.sct, features=c("percent.mito","percent.ribo"))
FeaturePlot(immune.combined.sct, features=c("nCount_RNA","nFeature_RNA"))


# CELL CLUSTER ANNOTATION
# Retrieve annotations from Pilot Study (Sage et al. 2022)
anno <- read.csv("cellMarkersFromPilot_top10.csv")
anno %>%
  kable(format = "html", col.names = colnames(anno)) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "300px")

#The scSorter method is based on a semi-supervised learning algorithm. It uses pre-processed data as input. 
#Therefore, we used the previously normalized, scaled and transformed expression matrix generated by SCT method as input.
#The top 3000 highly variable genes selected by SCT method are also given as input to scSorter.
exprMat <- immune.combined.sct@assays$SCT@data
topgenes <- features
topgene_filter = rowSums(exprMat[topgenes, ]!=0) > ncol(exprMat)*.1
topgenes = topgenes[topgene_filter]

## At last, we subset the preprocessed expression data and run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
exprMat = exprMat[rownames(exprMat) %in% picked_genes, ]

# Application of scSorter 
## The marker genes are fit on the expression matrix. 
## The predicted cell type assignments are stored in the Pred_Type vector.
rts <- scSorter(exprMat, anno)
immune.combined.sct$cell_type <- rts$Pred_Type
## Cell type distribution
table(rts$Pred_Type)
## Add cell type information to the metadata in the original Seurat object
cell_type.info <- data.frame(cell_type = immune.combined.sct$cell_type, row.names= colnames(immune.combined.sct))
immune.combined.sct <- AddMetaData(object = immune.combined.sct, metadata = cell_type.info)
## Convert cell type to factor
immune.combined.sct$cell_type <- as.factor(immune.combined.sct$cell_type) 
Idents(immune.combined.sct) <- "cell_type"
immune.combined.sct$ident <- NULL

#Stacked bar plot
#Analyze the cell type composition across different disease groups and store results in the data frame.
df <- data.frame(immune.combined.sct@meta.data$cell_type, immune.combined.sct@meta.data$disease_state) 
colnames(df) <- c("Cell_Type", "Condition")
df <- df %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb/C*100) 

df<-df[df$Cell_Type!="Unknown",]

xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10)) 
ggplot(df, aes(fill=Condition, y=Percent, x=Cell_Type)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_viridis(discrete = T) + xlab("") + xtheme + 
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Confirm automated annotation by visualizing the expression of canonical marker genes
FeaturePlot(immune.combined.sct, features = c("JCHAIN","MS4A1","TCF4","FSCN1",
                                            "CCR7","CD83","RET.1","LTC4S",
                                            "MS4A2","APOE","C1QA","C1QB",
                                            "IL1RN","CXCL2","SNX10","CCL5",
                                            "GZMA","CD3E"), pt.size = 0.2, ncol = 5) + theme(legend.position = 'none', 
                                            axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- DimPlot(immune.combined.sct, group.by = "cell_type", reduction = "umap", label = FALSE, cols= c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7", "#00b0e6"))
p1 + p4

# Create vectors of marker genes(features) for cell types 
feat_mac <- c("CD68","CD163") # markers for monocytes-macrophages (Mo/Ma)
feat_T <- c("CD2","CD3D","CD3E","CD3G") # markers for T cells
feat_B <- c("MS4A1","CD79A","CD79B") # markers for B/Plasma cells
feat_neut <- c("RGS2","ILT11B","CSF3R") # markers for neutrophils
feat_mast <- c("LTC4S","HPGDS","GCSAML","MS4A2") # markers for mast cells
feat_DC <- c("CD83","CCR7")# markers for dendritic cells 
feat_eos <- c("CEBPE","ALOX15","LRRC17","LOC100068723","THBS4","TMIGD3","LOC100066165","MYCT1","LOC100066165") # markers for eosinophils
# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_mac), name = "feat_mac") 
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_T), name = "feat_T") 
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_B), name = "feat_B") 
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_neut), name = "feat_neut")
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_mast), name = "feat_mast") 
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(feat_DC), name = "feat_DC")
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
immune.combined.sct <- AddModuleScore(immune.combined.sct,features = list(eq.s.genes), name = "feat_eq.s.genes")


# Merge the cell clusters into the major cell types previously identified
sc22 <- immune.combined.sct
Idents(sc22) <- "integrated_snn_res.0.6"
sc22 <- RenameIdents(object=sc22,'0'="Mo/Ma",'4'="Mo/Ma", '6'="Mo/Ma",'11'="Mo/Ma",
                     '12'="Mo/Ma",'14'="Mo/Ma",'15'="Mo/Ma",
                     '1'="T cells",'2'="T cells",'5'="T cells",'7'="T cells",
                     '8'="T cells", '9'="T cells",'10'="T cells",'18'="T cells",
                     '3'="Neutrophils", 
                     '13'="Mast cells",
                     '16'="B cells",
                     '17'="Dendritic cells")

# Store the current identities in a new column of meta.data called MajCellType
sc22$MajCellType <- Idents(object = sc22)
DimPlot(sc22, reduction = "umap", group.by = "MajCellType", label=T)
table(Idents(sc22))

# MARKER GENES FOR CELL CLUSTERS
# 1 - Markers for the major cell types
Idents(sc22) <- "MajCellType"
markers_all_major <- FindAllMarkers(sc22, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_major <- subset(markers_all_major, markers_all_major$p_val_adj < 0.05) #filtering the non significant genes
saveRDS(markers_all_major, "markers_all_major.rds")
# Supplementary table 3:
write.csv(markers_all_major, 'markers_major.csv')

# 2 - Markers for each of the 19 cell clusters (0-18)
Idents(sc22) <- "integrated_snn_res.0.6"
markers_all <- FindAllMarkers(sc22, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all <- subset(markers_all, markers_all$p_val_adj < 0.05) #filtering the non significant genes
# Supplementary table 2:
write.csv(markers_all, 'markers_all_19.csv')

# CELL TYPE DISTRIBUTION ACROSS SAMPLES
sc_cell_info<-data.frame(sc22@meta.data$MajCellType,sc22@meta.data$sample_id, sc22@meta.data$disease_state)
colnames(sc_cell_info)<-c("cellType","sampleId","diseaseGroup")
head(sc_cell_info)
# run Propeller
propeller(clusters = sc_cell_info$cellType, sample = sc_cell_info$sampleId, 
          group = sc_cell_info$diseaseGroup)
write.table(x, file = "cellPropMajor.txt",quote = FALSE, sep = "\t")


# SAVE THE FINAL OBJECT
saveRDS(sc22, "sc22.rds")

