# Differential gene expression (DGE) analysis between the severe equine asthma (SEA) group and the control group
# using the package Nebula
# DGE analysis is conducted on major cell types, B cell subtypes, monocyte-macrophage subtypes, 
# neutrophil subtypes and T cell subtypes

# Set up environment
set.seed(42)
library(Seurat)
library(celldex)
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
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(scFeatureFilter)
library(glmGamPoi)
library(scran)
library(PCAtools)
library(metap)
library(UpSetR)
library(purrr)
library(nebula)

#********************MAJOR CELL TYPE DGE*****************************************

sc22 <- readRDS("sc22.rds")
DimPlot(sc22, reduction = "umap", group.by = "CellSubtype", label=T)

#B_Plasma 
sce_B_Plasma<-sce_no_unknown[,sce_no_unknown$cluster_id=='B_Plasma']
model.df = model.matrix(~sce_B_Plasma$group_id, data=sce_B_Plasma$group_id)
re = nebula(counts(sce_B_Plasma),sce_B_Plasma$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_B_Plasma$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_B_Plasma$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_B_Plasma_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_B_Plasma_nebula_res$gene <- genes

colnames(sce_B_Plasma_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_B_Plasma_nebula_res[sce_B_Plasma_nebula_res$logFC > abs(1.0) & sce_B_Plasma_nebula_res$padj <=0.05,])

#Dendritic 
sce_Dendritic<-sce_no_unknown[,sce_no_unknown$cluster_id=='Dendritic']
model.df = model.matrix(~sce_Dendritic$group_id, data=sce_Dendritic$group_id)
re = nebula(counts(sce_Dendritic),sce_Dendritic$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_Dendritic$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_Dendritic$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_Dendritic_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_Dendritic_nebula_res$gene <- genes

colnames(sce_Dendritic_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_Dendritic_nebula_res[sce_Dendritic_nebula_res$logFC > abs(1.0) & sce_Dendritic_nebula_res$padj <=0.05,])

#Mast 
sce_Mast<-sce_no_unknown[,sce_no_unknown$cluster_id=='Mast']
model.df = model.matrix(~sce_Mast$group_id, data=sce_Mast$group_id)
re = nebula(counts(sce_Mast),sce_Mast$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_Mast$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_Mast$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_Mast_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_Mast_nebula_res$gene <- genes

colnames(sce_Mast_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_Mast_nebula_res[sce_Mast_nebula_res$logFC > abs(1.0) & sce_Mast_nebula_res$padj <=0.05,])

# MoMa
sce_MoMa<-sce_no_unknown[,sce_no_unknown$cluster_id=='MoMa']
model.df = model.matrix(~sce_MoMa$group_id, data=sce_MoMa$group_id)
re = nebula(counts(sce_MoMa),sce_MoMa$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_MoMa$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_MoMa$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_MoMa_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_MoMa_nebula_res$gene <- genes

colnames(sce_MoMa_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_MoMa_nebula_res[sce_MoMa_nebula_res$logFC > abs(1.0) & sce_MoMa_nebula_res$padj <=0.05,])

#Neutrophils
sce_Neutrophils<-sce_no_unknown[,sce_no_unknown$cluster_id=='Neutrophils']
model.df = model.matrix(~sce_Neutrophils$group_id, data=sce_Neutrophils$group_id)
re = nebula(counts(sce_Neutrophils),sce_Neutrophils$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_Neutrophils$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_Neutrophils$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_Neutrophils_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_Neutrophils_nebula_res$gene <- genes

colnames(sce_Neutrophils_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_Neutrophils_nebula_res[sce_Neutrophils_nebula_res$logFC > abs(1.0) & sce_Neutrophils_nebula_res$padj <=0.05,])

# T cells
sce_Tcells<-sce_no_unknown[,sce_no_unknown$cluster_id=='T cells']
model.df = model.matrix(~sce_Tcells$group_id, data=sce_Tcells$group_id)
re = nebula(counts(sce_Tcells),sce_Tcells$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_Tcells$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_Tcells$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_Tcells_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_Tcells_nebula_res$gene <- genes

colnames(sce_Tcells_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_Tcells_nebula_res[sce_Tcells_nebula_res$logFC > abs(1.0) & sce_Tcells_nebula_res$padj <=0.05,])


#*************************B CELL SUBTYPES DGE*****************************************

sc22_B <- readRDS("sc22_B.rds")
DimPlot(sc22_B, reduction = "umap", group.by = "CellSubtype", label=T)

# Convert to SingleCellExperiment object
sce_B <- as.SingleCellExperiment(sc22_B,assay = "SCT")

# remove undetected genes
sce_B <- sce_B[rowSums(counts(sce_B) > 0) > 0, ]
dim(sce_B)

# remove cells with few or many detected genes
qc <- perCellQCMetrics(sce_B)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_B <- sce_B[, !ol]
dim(sce_B)

# remove lowly expressed genes
sce_B <- sce_B[rowSums(counts(sce_B) > 1) >= 10, ]
dim(sce_B)

assays(sce_B)$vstresiduals <- vst(counts(sce_B), verbosity = FALSE)$y

(sce_B <- prepSCE(sce_B, 
                  kid = "CellSubtype", # subpopulation assignments
                  gid = "disease_state",  # group IDs
                  sid = "sample_id",   # sample IDs 
                  drop = FALSE))  # don't drop all other colData columns

#*DGE analysis using Nebula
#B0
sce_B0<-sce_B[,sce_B$cluster_id=='B0']
model.df = model.matrix(~sce_B0$group_id, data=sce_B0$group_id)
re = nebula(counts(sce_B0),sce_B0$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_B0$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_B0$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_B0_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_B0_nebula_res$gene <- genes

colnames(sce_B0_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_B0_nebula_res[sce_B0_nebula_res$logFC > abs(1.0) & sce_B0_nebula_res$padj <=0.05,])

#B1
sce_B1<-sce_B[,sce_B$cluster_id=='B1']
model.df = model.matrix(~sce_B1$group_id, data=sce_B1$group_id)
re = nebula(counts(sce_B1),sce_B1$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_B1$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_B1$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_B1_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_B1_nebula_res$gene <- genes

colnames(sce_B1_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_B1_nebula_res[sce_B1_nebula_res$logFC > abs(1.0) & sce_B1_nebula_res$padj <=0.05,])


#B2
sce_B2<-sce_B[,sce_B$cluster_id=='B2']
model.df = model.matrix(~sce_B2$group_id, data=sce_B2$group_id)
re = nebula(counts(sce_B2),sce_B2$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_B2$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_B2$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_B2_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_B2_nebula_res$gene <- genes

colnames(sce_B2_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_B2_nebula_res[sce_B2_nebula_res$logFC > abs(1.0) & sce_B2_nebula_res$padj <=0.05,])


# Print tables
# Only DEGs with an adjusted P-value <0.05 and an absolute average log2 fold change > 1 are considered
x<-sce_B0_nebula_res[sce_B0_nebula_res$padj < 0.05 & abs(sce_B0_nebula_res$logFC) > 1.0,]
file="B0_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_B1_nebula_res[sce_B1_nebula_res$padj < 0.05 & abs(sce_B1_nebula_res$logFC) > 1.0,]
file="B1_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_B2_nebula_res[sce_B2_nebula_res$padj < 0.05 & abs(sce_B2_nebula_res$logFC) > 1.0,]
file="B2_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)


#*#*************************DC SUBTYPES DGE*****************************************

sc22_DC <- readRDS("sc22_DC.rds")
DimPlot(sc22_DC, reduction = "umap", group.by = "CellSubtype", label=T)

# Convert to SingleCellExperiment object
sce_DC<-as.SingleCellExperiment(sc22_DC,assay = "SCT")

# remove undetected genes
sce_DC <- sce_DC[rowSums(counts(sce_DC) > 0) > 0, ]
dim(sce_DC)

# remove cells with few or many detected genes
qc <- perCellQCMetrics(sce_DC)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_DC <- sce_DC[, !ol]
dim(sce_DC)

# remove lowly expressed genes
sce_DC <- sce_DC[rowSums(counts(sce_DC) > 1) >= 10, ]
dim(sce_DC)

assays(sce_DC)$vstresiduals <- vst(counts(sce_DC), verbosity = FALSE)$y

(sce_DC <- prepSCE(sce_DC, 
                   kid = "CellSubtype", # subpopulation assignments
                   gid = "disease_state",  # group IDs
                   sid = "sample_id",   # sample IDs 
                   drop = FALSE)) 

# DGE analysis using Nebula
#DC 0
sce_DC0<-sce_DC[,sce_DC$cluster_id=='DC 0']
model.df = model.matrix(~sce_DC0$group_id, data=sce_DC0$group_id)
re = nebula(counts(sce_DC0),sce_DC0$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_DC0$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_DC0$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_DC0_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_DC0_nebula_res$gene <- genes

colnames(sce_DC0_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_DC0_nebula_res[sce_DC0_nebula_res$logFC > abs(1.0) & sce_DC0_nebula_res$padj <=0.05,])
write.table(sce_DC0_nebula_res,file="DC0_nebula_SEA_CTL.csv",sep="\t",row.names=F,quote=F)

#DC 1
sce_DC1<-sce_DC[,sce_DC$cluster_id=='DC 1']
model.df = model.matrix(~sce_DC1$group_id, data=sce_DC1$group_id)
re = nebula(counts(sce_DC1),sce_DC1$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_DC1$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_DC1$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_DC1_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_DC1_nebula_res$gene <- genes

colnames(sce_DC1_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_DC1_nebula_res[sce_DC1_nebula_res$logFC > abs(1.0) & sce_DC1_nebula_res$padj <=0.05,])


#DC 2
sce_DC2<-sce_DC[,sce_DC$cluster_id=='DC 2']
model.df = model.matrix(~sce_DC2$group_id, data=sce_DC2$group_id)
re = nebula(counts(sce_DC2),sce_DC2$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_DC2$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_DC2$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_DC2_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_DC2_nebula_res$gene <- genes

colnames(sce_DC2_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_DC2_nebula_res[sce_DC2_nebula_res$logFC > abs(1.0) & sce_DC2_nebula_res$padj <=0.05,])


#DC 3
sce_DC3<-sce_DC[,sce_DC$cluster_id=='DC 3']
model.df = model.matrix(~sce_DC3$group_id, data=sce_DC3$group_id)
re = nebula(counts(sce_DC3),sce_DC3$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_DC3$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_DC3$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_DC3_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_DC3_nebula_res$gene <- genes

colnames(sce_DC3_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_DC3_nebula_res[sce_DC3_nebula_res$logFC > abs(1.0) & sce_DC3_nebula_res$padj <=0.05,])

# Print tables
# Only DEGs with an adjusted P-value <0.05 and an absolute average log2 fold change > 1 are considered
x<-sce_DC0_nebula_res[sce_DC0_nebula_res$padj < 0.05 & abs(sce_DC0_nebula_res$logFC) > 1.0,]
file="DC0_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)



#*************************MO/MA SUBTYPES DGE*****************************************

sc22_moma <- readRDS("sc22_moma.rds")
DimPlot(sc22_moma, reduction = "umap", group.by = "CellSubtype", label=T)

# Convert to SingleCellExperiment object
sc22_moma<- RunTSNE(sc22_moma, reduction = "pca", dims = 1:30)
sce_moma<-as.SingleCellExperiment(sc22_moma,assay = "SCT")

# remove undetected genes
sce_moma <- sce_moma[rowSums(counts(sce_moma) > 0) > 0, ]
dim(sce_moma)

# remove cells with few or many detected genes
qc <- perCellQCMetrics(sce_moma)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_moma <- sce_moma[, !ol]
dim(sce_moma)

# remove lowly expressed genes
sce_moma <- sce_moma[rowSums(counts(sce_moma) > 1) >= 10, ]
dim(sce_moma)
assays(sce_moma)$vstresiduals <- vst(counts(sce_moma), verbosity = FALSE)$y
(sce_moma <- prepSCE(sce_moma, 
                     kid = "CellSubtype", # subpopulation assignments
                     gid = "disease_state",  # group IDs
                     sid = "sample_id",   # sample IDs 
                     drop = FALSE))  # don't drop all other colData columns

# DGE analysis using Nebula
## Mo/Ma 0
sce_moma0<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 0']
model.df = model.matrix(~sce_moma0$group_id, data=sce_moma0$group_id)
re = nebula(counts(sce_moma0),sce_moma0$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma0$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma0$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma0_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma0_nebula_res$gene <- genes

colnames(sce_moma0_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma0_nebula_res[sce_moma0_nebula_res$logFC > abs(1.0) & sce_moma0_nebula_res$padj <=0.05,])

## Mo/Ma 1
sce_moma1<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 1']
model.df = model.matrix(~sce_moma1$group_id, data=sce_moma1$group_id)
re = nebula(counts(sce_moma1),sce_moma1$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma1$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma1$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma1_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma1_nebula_res$gene <- genes

colnames(sce_moma1_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma1_nebula_res[sce_moma1_nebula_res$logFC > abs(1.0) & sce_moma1_nebula_res$padj <=0.05,])

## Mo/Ma 2
sce_moma2<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 2']
model.df = model.matrix(~sce_moma2$group_id, data=sce_moma2$group_id)
re = nebula(counts(sce_moma2),sce_moma2$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma2$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma2$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma2_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma2_nebula_res$gene <- genes

colnames(sce_moma2_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma2_nebula_res[sce_moma2_nebula_res$logFC > abs(1.0) & sce_moma2_nebula_res$padj <=0.05,])

## Mo/Ma 3
sce_moma3<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 3']
model.df = model.matrix(~sce_moma3$group_id, data=sce_moma3$group_id)
re = nebula(counts(sce_moma3),sce_moma3$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma3$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma3$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma3_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma3_nebula_res$gene <- genes

colnames(sce_moma3_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma3_nebula_res[sce_moma3_nebula_res$logFC > abs(1.0) & sce_moma3_nebula_res$padj <=0.05,])

## Mo/Ma 4
sce_moma4<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 4']
model.df = model.matrix(~sce_moma4$group_id, data=sce_moma4$group_id)
re = nebula(counts(sce_moma4),sce_moma4$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma4$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma4$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma4_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma4_nebula_res$gene <- genes

colnames(sce_moma4_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma4_nebula_res[sce_moma4_nebula_res$logFC > abs(1.0) & sce_moma4_nebula_res$padj <=0.05,])

## Mo/Ma 5
sce_moma5<-sce_moma[,sce_moma$cluster_id=='Mo/Ma 5']
model.df = model.matrix(~sce_moma5$group_id, data=sce_moma5$group_id)
re = nebula(counts(sce_moma5),sce_moma5$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_moma5$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_moma5$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_moma5_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_moma5_nebula_res$gene <- genes

colnames(sce_moma5_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_moma5_nebula_res[sce_moma5_nebula_res$logFC > abs(1.0) & sce_moma5_nebula_res$padj <=0.05,])


# Print tables
# Only DEGs with an adjusted P-value <0.05 and an absolute average log2 fold change > 1 are considered

x<-sce_moma0_nebula_res[sce_moma0_nebula_res$padj < 0.05 & abs(sce_moma0_nebula_res$logFC) > 1.0,]
file="moma0_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_moma1_nebula_res[sce_moma1_nebula_res$padj < 0.05 & abs(sce_moma1_nebula_res$logFC) > 1.0,]
file="moma1_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_moma2_nebula_res[sce_moma2_nebula_res$padj < 0.05 & abs(sce_moma2_nebula_res$logFC) > 1.0,]
file="moma2_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_moma3_nebula_res[sce_moma3_nebula_res$padj < 0.05 & abs(sce_moma3_nebula_res$logFC) > 1.0,]
file="moma3_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_moma4_nebula_res[sce_moma4_nebula_res$padj < 0.05 & abs(sce_moma4_nebula_res$logFC) > 1.0,]
file="moma4_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_moma5_nebula_res[sce_moma5_nebula_res$padj < 0.05 & abs(sce_moma5_nebula_res$logFC) > 1.0,]
file="moma5_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)



#*************************NEUTROPHIL SUBTYPES DGE*****************************************

sc22_ne <- readRDS("sc22_ne.rds")
DimPlot(sc22_ne, reduction = "umap", group.by = "CellSubtype", label=T)

# Convert to SingleCellExperiment object
sce_ne<-as.SingleCellExperiment(sc22_ne,assay = "SCT")

# remove undetected genes
sce_ne <- sce_ne[rowSums(counts(sce_ne) > 0) > 0, ]
dim(sce_ne)

# remove cells with few or many detected genes
qc <- perCellQCMetrics(sce_ne)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_ne <- sce_ne[, !ol]
dim(sce_ne)

# remove lowly expressed genes
sce_ne <- sce_ne[rowSums(counts(sce_ne) > 1) >= 10, ]
dim(sce_ne)

assays(sce_ne)$vstresiduals <- vst(counts(sce_ne), verbosity = FALSE)$y

(sce_ne <- prepSCE(sce_ne, 
                   kid = "CellSubtype", # subpopulation assignments
                   gid = "disease_state",  # group IDs
                   sid = "sample_id",   # sample IDs 
                   drop = FALSE))  # don't drop all other colData columns

# DGE analysis using Nebula
## Neu 0
sce_neu0<-sce_ne[,sce_ne$cluster_id=='Neu 0']
model.df = model.matrix(~sce_neu0$group_id, data=sce_neu0$group_id)
re = nebula(counts(sce_neu0),sce_neu0$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_neu0$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_neu0$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_neu0_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_neu0_nebula_res$gene <- genes

colnames(sce_neu0_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_neu0_nebula_res[sce_neu0_nebula_res$logFC > abs(1.0) & sce_neu0_nebula_res$padj <=0.05,])

## Neu 1
sce_neu1<-sce_ne[,sce_ne$cluster_id=='Neu 1']
model.df = model.matrix(~sce_neu1$group_id, data=sce_neu1$group_id)
re = nebula(counts(sce_neu1),sce_neu1$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_neu1$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_neu1$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_neu1_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_neu1_nebula_res$gene <- genes

colnames(sce_neu1_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_neu1_nebula_res[sce_neu1_nebula_res$logFC > abs(1.0) & sce_neu1_nebula_res$padj <=0.05,])

## Neu 2
sce_neu2<-sce_ne[,sce_ne$cluster_id=='Neu 2']
model.df = model.matrix(~sce_neu2$group_id, data=sce_neu2$group_id)
re = nebula(counts(sce_neu2),sce_neu2$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_neu2$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_neu2$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_neu2_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_neu2_nebula_res$gene <- genes

colnames(sce_neu2_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_neu2_nebula_res[sce_neu2_nebula_res$logFC > abs(1.0) & sce_neu2_nebula_res$padj <=0.05,])


# Print tables
# Only DEGs with an adjusted P-value <0.05 and an absolute average log2 fold change > 1 are considered

x<-sce_neu0_nebula_res[sce_neu0_nebula_res$padj < 0.05 & abs(sce_neu0_nebula_res$logFC) > 1.0,]
file="neu0_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_neu1_nebula_res[sce_neu1_nebula_res$padj < 0.05 & abs(sce_neu1_nebula_res$logFC) > 1.0,]
file="neu1_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_neu2_nebula_res[sce_neu2_nebula_res$padj < 0.05 & abs(sce_neu2_nebula_res$logFC) > 1.0,]
file="neu2_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)


#*************************T CELL SUBTYPES DGE*****************************************

sc22_T <- readRDS("sc22_T.rds")
DimPlot(sc22_T, reduction = "umap", group.by = "CellSubtype", label=T)

# Convert to a SingleCellExperiment object
sce_T<-as.SingleCellExperiment(sc22_T,assay = "SCT")

# remove undetected genes
sce_T <- sce_T[rowSums(counts(sce_T) > 0) > 0, ]
dim(sce_T)

# remove cells with few or many detected genes
qc <- perCellQCMetrics(sce_T)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_T <- sce_T[, !ol]
dim(sce_T)

# remove lowly expressed genes
sce_T <- sce_T[rowSums(counts(sce_T) > 1) >= 10, ]
dim(sce_T)

assays(sce_T)$vstresiduals <- vst(counts(sce_T), verbosity = FALSE)$y

(sce_T <- prepSCE(sce_T, 
                  kid = "CellSubtype", # subpopulation assignments
                  gid = "disease_state",  # group IDs
                  sid = "sample_id",   # sample IDs 
                  drop = FALSE))  # don't drop all other colData columns

# DGE analysis using Nebula
## T0
sce_T0<-sce_T[,sce_T$cluster_id=='T0']
model.df = model.matrix(~sce_T0$group_id, data=sce_T0$group_id)
re = nebula(counts(sce_T0),sce_T0$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T0$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T0$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T0_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T0_nebula_res$gene <- genes

colnames(sce_T0_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T0_nebula_res[sce_T0_nebula_res$logFC > abs(1.0) & sce_T0_nebula_res$padj <=0.05,])

## T1
sce_T1<-sce_T[,sce_T$cluster_id=='T1']
model.df = model.matrix(~sce_T1$group_id, data=sce_T1$group_id)
re = nebula(counts(sce_T1),sce_T1$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T1$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T1$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T1_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T1_nebula_res$gene <- genes

colnames(sce_T1_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T1_nebula_res[sce_T1_nebula_res$logFC > abs(1.0) & sce_T1_nebula_res$padj <=0.05,])

## T2
sce_T2<-sce_T[,sce_T$cluster_id=='T2']
model.df = model.matrix(~sce_T2$group_id, data=sce_T2$group_id)
re = nebula(counts(sce_T2),sce_T2$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T2$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T2$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T2_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T2_nebula_res$gene <- genes

colnames(sce_T2_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T2_nebula_res[sce_T2_nebula_res$logFC > abs(1.0) & sce_T2_nebula_res$padj <=0.05,])

## T3
sce_T3<-sce_T[,sce_T$cluster_id=='T3']
model.df = model.matrix(~sce_T3$group_id, data=sce_T3$group_id)
re = nebula(counts(sce_T3),sce_T3$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T3$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T3$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T3_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T3_nebula_res$gene <- genes

colnames(sce_T3_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T3_nebula_res[sce_T3_nebula_res$logFC > abs(1.0) & sce_T3_nebula_res$padj <=0.05,])

## T4
sce_T4<-sce_T[,sce_T$cluster_id=='T4']
model.df = model.matrix(~sce_T4$group_id, data=sce_T4$group_id)
re = nebula(counts(sce_T4),sce_T4$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T4$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T4$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T4_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T4_nebula_res$gene <- genes

colnames(sce_T4_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T4_nebula_res[sce_T4_nebula_res$logFC > abs(1.0) & sce_T4_nebula_res$padj <=0.05,])

## T5
sce_T5<-sce_T[,sce_T$cluster_id=='T5']
model.df = model.matrix(~sce_T5$group_id, data=sce_T5$group_id)
re = nebula(counts(sce_T5),sce_T5$sample_id,pred=model.df)
pvals <- re$summary$`p_sce_T5$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T5$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T5_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T5_nebula_res$gene <- genes

colnames(sce_T5_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T5_nebula_res[sce_T5_nebula_res$logFC > abs(1.0) & sce_T5_nebula_res$padj <=0.05,])

## T6
sce_T6<-sce_T[,sce_T$cluster_id=='T6']
model.df = model.matrix(~sce_T6$group_id, data=sce_T6$group_id)
re = nebula(counts(sce_T6),sce_T6$sample_id,pred=model.df) #raises a waring:In sqrt(re_all[, i]) : NaNs produced
pvals <- re$summary$`p_sce_T6$group_idSEA`
genes <- re$summary$gene
logfc <- re$summary$`logFC_sce_T6$group_idSEA`
padj <- p.adjust(pvals,method = "BH")

sce_T6_nebula_res <- data.frame(cbind(logfc,pvals,padj))
sce_T6_nebula_res$gene <- genes

colnames(sce_T6_nebula_res) <- c("logFC","pvalue","padj","gene")
dim(sce_T6_nebula_res[sce_T6_nebula_res$logFC > abs(1.0) & sce_T6_nebula_res$padj <=0.05,])

# Print tables
# Only DEGs with an adjusted P-value <0.05 and an absolute average log2 fold change > 1 are considered

x<-sce_T0_nebula_res[sce_T0_nebula_res$padj < 0.05 & abs(sce_T0_nebula_res$logFC) > 1.0,]
file="T0_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T1_nebula_res[sce_T1_nebula_res$padj < 0.05 & abs(sce_T1_nebula_res$logFC) > 1.0,]
file="T1_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T2_nebula_res[sce_T2_nebula_res$padj < 0.05 & abs(sce_T2_nebula_res$logFC) > 1.0,]
file="T2_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T3_nebula_res[sce_T3_nebula_res$padj < 0.05 & abs(sce_T3_nebula_res$logFC) > 1.0,]
file="T3_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T4_nebula_res[sce_T4_nebula_res$padj < 0.05 & abs(sce_T4_nebula_res$logFC) > 1.0,]
file="T4_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T5_nebula_res[sce_T5_nebula_res$padj < 0.05 & abs(sce_T5_nebula_res$logFC) > 1.0,]
file="T5_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)

x<-sce_T6_nebula_res[sce_T6_nebula_res$padj < 0.05 & abs(sce_T6_nebula_res$logFC) > 1.0,]
file="T6_nebula_SEA_CTL_filtered.csv"
write.table(x,file=file,sep="\t",row.names=F,quote=F)


