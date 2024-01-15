# Single-cell profiling of bronchoalveolar cells reveals a Th17 signature in neutrophilic severe equine asthma
# Sophie E. Sage, Tosso Leeb, Vidhya Jagannathan, Vinzenz Gerber
# First published: 28 December 2023 in Immunology, Wiley
# https://doi.org/10.1111/imm.13745

# SCRIPT USED TO PRODUCE FIGURES (MANUSCRIPT AND SUPPLEMENTARY FILES)

# Set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set up environment
set.seed(42)
library(Seurat)
library(clustree)
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
library(magrittr)
library(scFeatureFilter)
library(glmGamPoi)
library(scran)
library(metap)
library(reshape2)
library(tidyr)
library(grid)
library(gridExtra)
library(EnhancedVolcano)

# Load the Seurat objects corresponding to each cell population
sc22 <- readRDS('sc22.rds')
sc22_B <- readRDS('sc22_B.rds')
sc22_ne <- readRDS('sc22_ne.rds')
sc22_T <- readRDS('sc22_T.rds')
sc22_moma <- readRDS('sc22_moma.rds')
sc22_DC <- readRDS('sc22_DC.rds')

# Create custom function to label clusters
theme.c <- theme(
  axis.line.x.bottom = element_line(color = 'black'),
  axis.line.y.left   = element_line(color = 'black'),
  aspect.ratio = 1,
  panel.background = element_blank(),
  panel.border = element_blank()
)

GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}
# Custom function to add circles labels to UMAP visualizations
custom.LabelClusters <- function(
    plot, # Use DimPlot to generate base ggplot to apply function
    id,   # The seurat cluster identifier
    clusters = NULL,
    labels = NULL,
    split.by = NULL,
    repel = F,
    colors = colors,
    circle.size = circle.size,
    text.size = text.size,
    ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape=21, fill = "white", stroke = 0.66, col = "black",
    mapping = aes(x = !!sym(xynames['x']), y = !!sym(xynames['y'])),
    ...
  ) + geom_text( 
    size = text.size,
    data = labels.loc, col = "black",
    mapping = aes(x = !!sym(xynames['x']), y = !!sym(xynames['y']), label = !!sym(id)),
    ...
  )
  return(plot)
}


#************************************FIGURE 1**************************************************

# FIGURE 1A

# Define cluster colors (1 color = 1 cell type, using different shades for the clusters) 
colors3<-c("#FF6600","#0080FE","#78BFEA","#9675B4","#F57D20","#00BED3","#FFCC66",
           "#5BC1BF","#678297","#489ECE","#84C8E2","#FFCC33","#FF9933", "#FF0000",
           "#FFCC00","#F05729", "#9ACA3A","#996600","#5BC1BF")

Idents(sc22) <- "integrated_snn_res.0.6"
p_sc22 <-DimPlot(object = sc22, pt.size = 0.01,label = F,group.by = "integrated_snn_res.0.6") +
  scale_color_manual(values = colors3) + scale_fill_manual(values = colors3) + 
  xlab(expression('UMAP'[1])) + ylab(expression('UMAP'[2])) +
  theme(legend.position = "none") + labs(title="")
p_sc22 <- custom.LabelClusters(plot = p_sc22, 
                               id = "integrated_snn_res.0.6", 
                               clusters = levels(sc22@meta.data[,"integrated_snn_res.0.6"]),
                               circle.size = 6, 
                               color = "black", 
                               text.size = 3,
                               shape = 21,
                               fill = colors3,
                               repel = T)

#Print plot
dpi=300
tiff(file='sc22_UMAP_major.tiff', width = dpi*7, height = dpi*6, units = "px",res = dpi)
print(p_sc22)
dev.off()
## Circling of the major cell types performed with Sketch (www.sketch.com)


# FIGURE 1B

Idents(sc22) <- "MajCellType"
pt_major <- table(Idents(sc22), sc22$orig.ident) 
pt_major <- as.data.frame(pt_major)
pt_major$Var1 <- as.character(pt_major$Var1)
colnames(pt_major)<-c("cell_type","horse","Freq")
pt_major$cell_type <- as.factor(pt_major$cell_type)
pt_major$cell_type <- recode(pt_major$cell_type,"Dendritic cells"="DCs")

# Create a new column 'status' based on 'horse_id'
# pt_major$status <- ifelse(pt_major$horse %in% c("30", "33", "44","52","56","68"), "Control", "SEA") 
pt_major$status <- ifelse(pt_major$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
pt_major$status <- as.factor(pt_major$status)

# Create a new column 'horse_id' based on 'status'
pt_major$horse_id <- ifelse(pt_major$status %in% "SEA", paste0("A", pt_major$horse), paste0("C", pt_major$horse))

# Create a new column "prop" with cell type proportion for each horse
pt_major <- pt_major %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

major_boxplot <- ggplot(pt_major, aes(x = cell_type, y = prop, fill = status)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

major_boxplot <- major_boxplot +
  annotate("segment", x = 0.85, xend = 1.15, y = 5, yend = 5, colour = "black") +
  annotate("text", x = 1, y = 6, label = "***") +
  annotate("segment", x = 4.85, xend = 5.15, y = 27, yend = 27, colour = "black") +
  annotate("text", x = 5, y = 28, label = "***") 

# Print plot
tiff(file='major_boxplot.tiff', width = 2300, height = 1500, units = "px",res = 300)
major_boxplot
dev.off()


# FIGURE 1C

# Reorder levels of 'horse_id' to have control horses on the left
pt_major$horse_id <- factor(pt_major$horse_id, levels = c(paste0("C", unique(pt_major$horse)), paste0("A", unique(pt_major$horse))))

# Create plot
major_plot <- ggplot(pt_major, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cell type") +
  scale_fill_manual(values = c("#99CC33", "#FFFF66", "#FF0000", "#FFCC66", "#9675B4", "#99CCFF")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='major_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
major_plot
dev.off()


#************************************FIGURE 2**************************************************

# FIGURE 2A

Idents(sc22_B) <- "SCT_snn_res.0.3"
cellnumber <- dim(sc22_B)[2]
group.ident <-"SCT_snn_res.0.3"

# Define cluster colors 
B_colors<-c("#03C04A","#466D1D","#006633")
pt.size=0.01
label=FALSE

# Create plot
p_B <- 
  DimPlot(object = sc22_B, 
          pt.size = pt.size,
          label = F,
          group.by = group.ident) +
  scale_color_manual(values = B_colors)  + 
  scale_fill_manual(values = B_colors)  + 
  theme.c +
  xlab(expression('UMAP'[1])) +
  ylab(expression('UMAP'[2])) +
  theme(legend.position = "noB",axis.title = element_text(size = 10),axis.text = element_text(size=7)) + labs(title = "")

# Apply custom labels
p_B <- custom.LabelClusters(plot = p_B, 
                            id = group.ident, 
                            clusters = levels(sc22_B@meta.data[,group.ident]),
                            circle.size = 6, 
                            text.size = 4, 
                            colors = B_colors,
                            repel = T)

# Print plot
tiff(file='B_UMAP.tiff', width = 2500, height = 2000, units = "px",res = 300)
print(p_B)
dev.off()


# FIGURE 2B

Idents(sc22_B) <- "CellSubtype"
pt_B <- table(Idents(sc22_B), sc22_B$orig.ident) 
pt_B <- as.data.frame(pt_B)
pt_B$Var1 <- as.character(pt_B$Var1)
colnames(pt_B)<-c("cell_type","horse","Freq")

# Create a new column 'status' based on 'horse_id'
# pt_major$status <- ifelse(pt_major$horse %in% c("30", "33", "44","52","56","68"), "Control", "SEA") 
pt_B$status <- ifelse(pt_B$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
# Create a new column 'horse_id' based on 'status'
pt_B$horse_id <- ifelse(pt_B$status %in% "SEA", paste0("A", pt_B$horse), paste0("C", pt_B$horse))

# Create a new column "prop" with cell type proportion for each horse
pt_B <- pt_B %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

# Create plot
B_boxplot <- ggplot(pt_B, aes(x = cell_type, y = prop, fill = status)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

B_boxplot <- B_boxplot + 
  annotate("segment", x = 2.85, xend = 3.15, y = 30, yend = 30, colour = "black") +
  annotate("text", x = 3, y = 31, label = "*") 

#Print plot
tiff(file='B_boxplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
B_boxplot
dev.off()

# FIGURE 3C

# Reorder levels of 'horse_id' to have control horses on the left
pt_B$horse_id <- factor(pt_B$horse_id, levels = c(paste0("C", unique(pt_B$horse)), paste0("A", unique(pt_B$horse))))

# Create plot
B_plot <- ggplot(pt_B, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values = brewer.pal(3, "Greens")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='B_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
B_plot
dev.off()


#************************************FIGURE 3**************************************************

# FIGURE 3A

Idents(sc22_ne) <- "integrated_snn_res.0.15"
cellnumber <- dim(sc22_ne)[2]
group.ident <-"integrated_snn_res.0.15"

# Define cluster colors 
ne_colors <- c("#9966CC","#663399","#330066")
pt.size=0.01
label=FALSE

# Create plot
p_ne <- 
  DimPlot(object = sc22_ne, 
          pt.size = pt.size,
          label = F,
          group.by = group.ident) +
  scale_color_manual(values = ne_colors)  + 
  scale_fill_manual(values = ne_colors)  + 
  theme.c +
  xlab(expression('UMAP'[1])) +
  ylab(expression('UMAP'[2])) +
  theme(legend.position = "none",axis.title = element_text(size = 10),axis.text = element_text(size=7)) + labs(title = "")

# Apply custom labels
p_ne <- custom.LabelClusters(plot = p_ne, 
                             id = group.ident, 
                             clusters = levels(sc22_ne@meta.data[,group.ident]),
                             circle.size = 6, 
                             text.size = 4, 
                             colors = ne_colors,
                             repel = T)

# Print plot
tiff(file='ne_UMAP.tiff', width = 2500, height = 2000, units = "px",res = 300)
print(p_ne)
dev.off()

# FIGURE 3B

Idents(sc22_ne) <- "CellSubtype"
pt_ne <- table(Idents(sc22_ne), sc22_ne$orig.ident) 
pt_ne <- as.data.frame(pt_ne)
pt_ne$Var1 <- as.character(pt_ne$Var1)
colnames(pt_ne)<-c("cell_type","horse","Freq")

# Create a new column 'status' based on 'horse_id'
pt_ne$status <- ifelse(pt_ne$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
# Create a new column 'horse_id' based on 'status'
pt_ne$horse_id <- ifelse(pt_ne$status %in% "SEA", paste0("A", pt_ne$horse), paste0("C", pt_ne$horse))

# Create a new column "prop" with cell type proportion for each horse
pt_ne <- pt_ne %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

# Create plot
ne_boxplot <- ggplot(pt_ne, aes(x = cell_type, y = prop, fill = status)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='neu_boxplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
ne_boxplot
dev.off()


# FIGURE 3C

# Reorder levels of 'horse_id' to have control horses on the left
pt_ne$horse_id <- factor(pt_ne$horse_id, levels = c(paste0("C", unique(pt_ne$horse)), paste0("A", unique(pt_ne$horse))))

# Create the plot
ne_plot <- ggplot(pt_ne, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values = brewer.pal(3, "Purples")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='neu_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
ne_plot
dev.off()

#************************************FIGURE 4**************************************************

# FIGURE 4A

Idents(sc22_T) <- "integrated_snn_res.0.3"
cellnumber <- dim(sc22_T)[2]
group.ident <-"integrated_snn_res.0.3"

# Define cluster colors 
T_colors <- c("#3399CC","#006699", "#3366CC","#99CCFF","#6699FF","#66CCFF","#6699FF")
pt.size=0.01
label=FALSE

# Create plot
p_T <- 
  DimPlot(object = sc22_T, 
          pt.size = pt.size,
          label = F,
          group.by = group.ident) +
  scale_color_manual(values = T_colors)  + 
  scale_fill_manual(values = T_colors)  + 
  theme.c +
  xlab(expression('UMAP'[1])) +
  ylab(expression('UMAP'[2])) +
  theme(legend.position = "noT",axis.title = element_text(size = 10),axis.text = element_text(size=7)) + labs(title = "")

# Apply custom labels
p_T <- custom.LabelClusters(plot = p_T, 
                            id = group.ident, 
                            clusters = levels(sc22_T@meta.data[,group.ident]),
                            circle.size = 6, 
                            text.size = 4, 
                            colors = T_colors,
                            repel = T)

# Print plot
tiff(file='T_UMAP.tiff', width = 2500, height = 2000, units = "px",res = 300)
print(p_T)
dev.off()


# FIGURE 4B

Idents(sc22_T) <- "CellSubtype"
pt_T <- table(Idents(sc22_T), sc22_T$orig.ident) 
pt_T <- as.data.frame(pt_T)
pt_T$Var1 <- as.character(pt_T$Var1)
colnames(pt_T)<-c("cell_type","horse","Freq")

# Create a new column 'status' based on 'horse_id'
pt_T$status <- ifelse(pt_T$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
# Create a new column 'horse_id' based on 'status'
pt_T$horse_id <- ifelse(pt_T$status %in% "SEA", paste0("A", pt_T$horse), paste0("C", pt_T$horse))

# Creating a new column "prop" with cell type proportion for each horse
pt_T <- pt_T %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

T_boxplot <- ggplot(pt_T, aes(x = cell_type, y = prop, fill = status)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='T_boxplot.tiff', width = 2300, height = 1500, units = "px",res = 300)
T_boxplot
dev.off()

# FIGURE 4C

# Reorder levels of 'horse_id' to have control horses on the left
pt_T$horse_id <- factor(pt_T$horse_id, levels = c(paste0("C", unique(pt_T$horse)), paste0("A", unique(pt_T$horse))))

# Create plot
T_plot <- ggplot(pt_T, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values = brewer.pal(7, "Blues")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='T_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
T_plot
dev.off()


# FIGURE 3D

# Extract the DEGs between asthmatic and control horses using Nebula
T_table<-read.table("Tcells_nebula_SEA_CTL_all.csv",header=T,sep="\t")
T_df<-data.frame(T_table$logFC,T_table$pvalue,T_table$padj)
row.names(T_df)<-T_table$gene
colnames(T_df)<-c("logFC","pvalue","padj")

# Create volcano plot
tiff(file='Volcano_T.tiff', width = 3000, height = 2000, units = "px",res = 300)
(EnhancedVolcano(T_df,
                 lab = rownames(T_df),
                 x = 'logFC',
                 y = 'padj',
                 title = '',
                 subtitle = '',
                 selectLab = c('IL17A', 'IL17F', 'IL21','CCL20'),
                 xlab = bquote(~Log[2]~ 'fold change'),
                 pCutoff = 0.05,#I changed 10e-2 to 0.05 (criterium we used for DEG)
                 FCcutoff = 1.0,
                 pointSize = 3.0,
                 labSize = 4.0,
                 colAlpha = 1,
                 # legendLabSize = 14,
                 legendPosition = 'none',
                 # legendIconSize = 3.0,
                 boxedLabels = TRUE,
                 drawConnectors = TRUE,
                 widthConnectors = 0.5,
                 colConnectors = 'black',
                 gridlines.major = TRUE,
                 gridlines.minor = FALSE,
                 border = 'partial',
                 borderWidth = 1.5,
                 borderColour = 'black')
) +
  theme(axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.line = element_line(linewidth = 0.5))
dev.off()


#************************************FIGURE 5**************************************************

# FIGURE 5A

Idents(sc22_moma) <- "integrated_snn_res.0.2"
cellnumber <- dim(sc22_moma)[2]
group.ident <-"integrated_snn_res.0.2"

# Define cluster colors 
moma_colors<-c("#ED7014","#FAC80A","#80400B","#EC9706","#FCAE1E","#D23B05")
pt.size=0.01
label=FALSE

# Create plot
p_moma <- 
  DimPlot(object = sc22_moma, 
          pt.size = pt.size,
          label = F,
          group.by = group.ident) +
  scale_color_manual(values = moma_colors)  + 
  scale_fill_manual(values = moma_colors)  + 
  theme.c +
  xlab(expression('UMAP'[1])) +
  ylab(expression('UMAP'[2])) +
  theme(legend.position = "nomoma",axis.title = element_text(size = 10),axis.text = element_text(size=7)) + labs(title = "")

# Apply custom labels
p_moma <- custom.LabelClusters(plot = p_moma, 
                               id = group.ident, 
                               clusters = levels(sc22_moma@meta.data[,group.ident]),
                               circle.size = 6, 
                               text.size = 4, 
                               colors = moma_colors,
                               repel = T)

# Print plot
tiff(file='moma_UMAP.tiff', width = 2500, height = 2000, units = "px",res = 300)
print(p_moma)
dev.off()


# FIGURE 5B

Idents(sc22_moma) <- "CellSubtype"
pt_moma <- table(Idents(sc22_moma), sc22_moma$orig.ident) 
pt_moma <- as.data.frame(pt_moma)
pt_moma$Var1 <- as.character(pt_moma$Var1)
colnames(pt_moma)<-c("cell_type","horse","Freq")

# Create a new column 'status' based on 'horse_id'
pt_moma$status <- ifelse(pt_moma$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
# Create a new column 'horse_id' based on 'status'
pt_moma$horse_id <- ifelse(pt_moma$status %in% "SEA", paste0("A", pt_moma$horse), paste0("C", pt_moma$horse))

# Create a new column "prop" with cell type proportion for each horse
pt_moma <- pt_moma %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

# Create plot
moma_boxplot <- ggplot(pt_moma, aes(x = cell_type, y = prop, fill = status)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='moma_boxplot.tiff', width = 2200, height = 1500, units = "px",res = 300)
moma_boxplot
dev.off()


# FIGURE 5C

# Reorder levels of 'horse_id' to have control horses on the left
pt_moma$horse_id <- factor(pt_moma$horse_id, levels = c(paste0("C", unique(pt_moma$horse)), paste0("A", unique(pt_moma$horse))))

# Create plot
moma_plot <- ggplot(pt_moma, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values = brewer.pal(6, "Oranges")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot

tiff(file='moma_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
moma_plot
dev.off()

# FIGURE 5D

# Extract the DEGs between asthmatic and control horses using Nebula
moma_table<-read.table("MoMa_nebula_SEA_CTL_all.csv",header=T,sep="\t")
moma_df<-data.frame(moma_table$logFC,moma_table$pvalue,moma_table$padj)
row.names(moma_df)<-moma_table$gene
colnames(moma_df)<-c("logFC","pvalue","padj")

# Create volcano plot
tiff(file='Volcano_moma.tiff', width = 3000, height = 2000, units = "px",res = 300)
(EnhancedVolcano(moma_df,
                 lab = rownames(moma_df),
                 x = 'logFC',
                 y = 'padj',
                 title = '',
                 subtitle = '',
                 selectLab = c('CXCL13', 'OLFM4', 'CHI3L1', 'LOC100061699'),#LOC100061699=S100A8
                 xlab = bquote(~Log[2]~ 'fold change'),
                 pCutoff = 0.05,#I changed 10e-2 to 0.05 (criterium we used for DEG)
                 FCcutoff = 1.0,
                 pointSize = 3.0,
                 labSize = 4.0,
                 colAlpha = 1,
                 legendPosition = 'none',
                 boxedLabels = TRUE,
                 drawConnectors = TRUE,
                 widthConnectors = 0.5,
                 colConnectors = 'black',
                 gridlines.major = TRUE,
                 gridlines.minor = FALSE,
                 border = 'partial',
                 borderWidth = 1.5,
                 borderColour = 'black')
) +
  theme(axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.line = element_line(linewidth = 0.5))
dev.off()


#************************************FIGURE 6**************************************************

# FIGURE 6A

Idents(sc22_DC) <- "integrated_snn_res.0.1"
cellnumber <- dim(sc22_DC)[2]
group.ident <-"integrated_snn_res.0.1"

# Define cluster colors 
DC_colors<-c("#ED7014","#FAC80A","#80400B","#EC9706")
pt.size=0.01
label=FALSE

# Create plot
p_DC <- 
  DimPlot(object = sc22_DC, 
          pt.size = pt.size,
          label = F,
          group.by = group.ident) +
  scale_color_manual(values = DC_colors)  + 
  scale_fill_manual(values = DC_colors)  + 
  theme.c +
  xlab(expression('UMAP'[1])) +
  ylab(expression('UMAP'[2])) +
  theme(legend.position = "noDC",axis.title = element_text(size = 10),axis.text = element_text(size=7)) + labs(title = "")

# Apply custom labels
p_DC <- custom.LabelClusters(plot = p_DC, 
                             id = group.ident, 
                             clusters = levels(sc22_DC@meta.data[,group.ident]),
                             circle.size = 6, 
                             text.size = 4, 
                             colors = DC_colors,
                             repel = T)

# Print plot
tiff(file='DC_UMAP.tiff', width = 2500, height = 2000, units = "px",res = 300)
print(p_DC)
dev.off()


# FIGURE 6B

Idents(sc22_DC) <- "CellSubtype"
pt_DC <- table(Idents(sc22_DC), sc22_DC$orig.ident) 
pt_DC <- as.data.frame(pt_DC)
pt_DC$Var1 <- as.character(pt_DC$Var1)
colnames(pt_DC)<-c("cell_type","horse","Freq")

# Create a new column 'status' based on 'horse_id'
pt_DC$status <- ifelse(pt_DC$horse %in% c("30", "33", "44","52","56","68"), "SEA", "Control")
# Create a new column 'horse_id' based on 'status'
pt_DC$horse_id <- ifelse(pt_DC$status %in% "SEA", paste0("A", pt_DC$horse), paste0("C", pt_DC$horse))
# Create a new column "prop" with cell type proportion for each horse
pt_DC <- pt_DC %>%
  group_by(horse) %>%
  mutate(
    total_count = sum(Freq),
    prop = Freq / total_count * 100
  ) %>%
  ungroup()

# Create plot
DC_boxplot <- ggplot(pt_DC, aes(x = cell_type, y = prop, fill = status)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Percentage", fill = "Group") +
  scale_fill_manual(values = c("#99CCFF", "#CC3333")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 0, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='DC_boxplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
DC_boxplot
dev.off()


# FIGURE 6C

# Reorder levels of 'horse_id' to have control horses on the left
pt_DC$horse_id <- factor(pt_DC$horse_id, levels = c(paste0("C", unique(pt_DC$horse)), paste0("A", unique(pt_DC$horse))))

# Create plot
DC_plot <- ggplot(pt_DC, aes(x = horse_id, y = Freq, fill = cell_type)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x = "", y = "Percentage", fill = "Cluster") +
  scale_fill_manual(values = c("#ED7014","#FAC80A","#80400B","#EC9706")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12))

# Print plot
tiff(file='DC_barplot.tiff', width = 2000, height = 1500, units = "px",res = 300)
DC_plot
dev.off()


#*****************************SUPPLEMENTARY FIGURE 1**************************************************
# In supplementary material & methods

# VIDHYA PLEASE ADD
Power analysis script

#*****************************SUPPLEMENTARY FIGURE 2**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 2A

tiff(file='sc22_featureplots.tiff', width = 4500, height = 2500, units = "px",res = 300)
(((FeaturePlot(sc22, features = "feat_T1")+ labs(title="T cell FES<br>(*CD2, CD3D, CD3E, CD3G*)"))+
    (FeaturePlot(sc22, features = "feat_mac1")+ labs(title="Mo/Ma FES<br>(*CD68, CD163*)"))+
    (FeaturePlot(sc22, features = "feat_neut1")+ labs(title="Neutrophil FES<br>(*CSF3R, LILRA5, RGS2*)"))+
    (FeaturePlot(sc22, features = "feat_mast1")+ labs(title="Mast cell FES<br>(*GCSAML, HPGDS, LTC4S,MS4A2*)"))+
    (FeaturePlot(sc22, features = "feat_B1")+ labs(title="B cell FES<br>(*CD79A, CD79B, MS4A1*)"))+
    (FeaturePlot(sc22, features = "feat_DC1")+ labs(title="DC FES<br>(*CCR7, CD83*)")))+
    plot_layout(ncol = 3)) &
  theme(plot.title=element_markdown(size=16, face = "bold"), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 2B

markers_all_major <- readRDS("markers_all_major.rds")
major_top5 <- markers_all_major %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_major), value = TRUE)), n = 5)

major_top5_cod <- markers_all_major %>% 
  group_by(cluster) %>%
  filter(!gene=="LOC102148763")%>% #removing non-coding gene
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_major), value = TRUE)), n = 5)
# rename genes with their NCBI annotation
major_top5_cod$gene <- recode(major_top5_cod$gene,"LOC102147726"="*IGLL1", "LOC100146200"="*OR7A189P",
                              "LOC100069985"="*CD177")
markers <- c("*IGLL1", "*OR7A189P","*CD177")
markers.loc <- c("LOC102147726","LOC100146200","LOC100069985")
Idents(sc22) <- "MajCellType"
sc22@meta.data[,markers]<-GetAssayData(sc22)[markers.loc,] %>% as.matrix %>% t

tiff(file='sc22_major_dotplot.tiff', width = 4000, height = 1000, units = "px",res = 300)
DotPlot(sc22, features = unique(major_top5_cod$gene)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face="italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
dev.off()


#*****************************SUPPLEMENTARY FIGURE 3**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 3A

tiff(file='neu_featureplots.tiff', width = 3000, height = 1250, units = "px",res = 300)
(((FeaturePlot(sc22_ne, features = "feat_infl1")+ labs(title="Pro-inflammatory genes FES<br>(*IL1A, IL1B, CXCL8*)"))+
    (FeaturePlot(sc22_ne, features = "ISG15")+ labs(title="*ISG15*"))) + 
    plot_layout(ncol = 2)) &
  theme(plot.title=element_markdown(size=16), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 3B

sc22_ne <- SetIdent(sc22_ne, value = sc22_ne$CellSubtype) 

markers_all_ne <- readRDS("markers_all_ne.rds")
ne_top5 <- markers_all_ne %>% group_by(cluster) %>%
  slice_max(get(grep("^avg_log", colnames(markers_all_ne), value = TRUE)), n = 5)

ne_top5_cod <- markers_all_ne %>% 
  group_by(cluster) %>%
  filter(!gene=="ND2") %>% #removing Mt gene
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_ne), value = TRUE)), n = 5)

ne_top5_cod$gene <- recode(ne_top5_cod$gene,"LOC100054211"="*H2BC21","LOC100050797"="*IFITM1")

ne_markers <- c("*H2BC21","*IFITM1")
ne_markers.loc <- c("LOC100054211","LOC100050797")
sc22_ne@meta.data[,ne_markers]<-GetAssayData(sc22_ne)[ne_markers.loc,] %>% as.matrix %>% t

tiff(file='neu_dotplot.tiff', width = 4000, height = 1000, units = "px",res = 300)
DotPlot(sc22_ne, features = unique(ne_top5_cod$gene), 
        idents = c("Neu 0","Neu 1","Neu 2")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face="italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) 
dev.off()


#*****************************SUPPLEMENTARY FIGURE 4**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 4A

tiff(file='B_featureplots.tiff', width = 3000, height = 2500, units = "px",res = 300)
((FeaturePlot(sc22_B, features = "feat_HLA1")+ labs(title="MHC II FES<br>(*DQB, DQB.1, DQA, DQA.1, DRA, DRB*)"))+
    (FeaturePlot(sc22_B, features = "feat_ab1")+ labs(title="Antibody production FES<br>(*TXNDC5, HSP90B1, TENT5C*)"))+
    (FeaturePlot(sc22_B, features = "feat_igm1")+ labs(title="IgM production FES<br>(*JCHAIN, MZB1*))"))+
    (FeaturePlot(sc22_B, features = "LGALS1")+ labs(title="*LGALS1*"))+
    plot_layout(ncol = 2)) &
  theme(plot.title=element_markdown(size=16), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 4B

sc22_B <- SetIdent(sc22_B, value = sc22_B$CellSubtype) 
markers_all_B <- readRDS("markers_all_B.rds")

B_top5 <- markers_all_B %>% group_by(cluster) %>%
  slice_max(get(grep("^avg_log", colnames(markers_all_B), value = TRUE)), n = 5)
B_top5$gene <- recode(B_top5$gene,"LOC102147726"="*IGLL1(1)","LOC100060608"="*IGLL1(2)","LOC111774805"="*IGLL1(3)",
                      "LOC100061331"="*MS4A4A")
B_markers <- c("*IGLL1(1)","*IGLL1(2)","*IGLL1(3)","*MS4A4A")
B_markers.loc <- c("LOC102147726","LOC100060608","LOC111774805","LOC100061331")
sc22_B@meta.data[,B_markers]<-GetAssayData(sc22_B)[B_markers.loc,] %>% as.matrix %>% t

tiff(file='B_dotplot.tiff', width = 3000, height = 1000, units = "px",res = 300)
DotPlot(sc22_B, features = unique(B_top5$gene), 
        idents = c("B0","B1","B2")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
dev.off()


#*****************************SUPPLEMENTARY FIGURE 5**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 5A

tiff(file='T_featureplots_col3.tiff', width = 4500, height = 3750, units = "px",res = 300)
(((FeaturePlot(sc22_T, features = "GZMA")+ labs(title="*GZMA*")) +
    (FeaturePlot(sc22_T, features = "GZMK")+ labs(title="*GZMK*"))+ 
    (FeaturePlot(sc22_T, features = "KLRD1")+ labs(title="*KLRD1*")))+
    (FeaturePlot(sc22_T, features = "KLRB1")+ labs(title="*KLRB1*"))+
    (FeaturePlot(sc22_T, features = "GNLY")+ labs(title="*GNLY*"))+
    (FeaturePlot(sc22_T, features = "ISG15")+ labs(title="*ISG15*"))+
    (FeaturePlot(sc22_T, features = "IL17A")+ labs(title="IL17A"))+
    (FeaturePlot(sc22_T, features = "feat_eq.g2m.genes1")+ labs(title="*G2M genes*")) + 
    plot_layout(ncol = 3)) &
  theme(plot.title=element_markdown(size=16), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 5B

markers_all_T <- readRDS("markers_all_T.rds")

T_top5_cod <- markers_all_T %>% 
  group_by(cluster) %>%
  filter(!gene=="LOC100052742" & !gene=="LOC100071401" & #removing non-coding genes
           !gene=="LOC111775812" & #removing snRNA gene
           !gene=="RPS12" & !gene=="RPS19")%>% #removing RP genes
  slice_max(get(grep("^avg_log", colnames(markers_all_T), value = TRUE)), n = 5)

T_top5_cod$gene <- recode(T_top5_cod$gene,"LOC100051986"="*GZMH", "LOC100065392"="*MTAP", 
                          "LOC100062823"="*KLRC1", "LOC100053968"="*H2BC20", "LOC100059091"="*TUBA1B" )

T_markers <- c("*GZMH", "*MTAP", "*KLRC1","*H2BC20", "*TUBA1B")
T_markers.loc <- c("LOC100051986", "LOC100065392","LOC100062823","LOC100053968","LOC100059091")
sc22_T <- SetIdent(sc22_T, value = sc22_T$CellSubtype) 
sc22_T_fig <- sc22_T
sc22_T_fig@meta.data[,T_markers]<-GetAssayData(sc22_T_fig)[T_markers.loc,] %>% as.matrix %>% t

tiff(file='T_dotplot.tiff', width = 4000, height = 1000, units = "px",res = 300)
DotPlot(sc22_T_fig, features = unique(T_top5_cod$gene), 
        idents = c("T0","T1","T2","T3","T4","T5","T6")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) 
dev.off()


#SUPPLEMENTARY FIGURE 5C

tiff(file='T_violinplot_col3.tiff', width = 3000, height = 1000, units = "px",res = 300)
VlnPlot(sc22_T, features = c("CD4","CD8A","CD8B"), pt.size = 0.1) +
  plot_layout(ncol = 3) &
  theme(axis.title=element_blank(), plot.title=element_text(size=14, face="italic"), legend.position = 'none') &
  scale_fill_brewer(palette="Blues", direction = "-1")
dev.off()


#SUPPLEMENTARY FIGURE 5D

tiff(file='sc22T_RP genes.tiff', width = dpi*12, height = dpi*5, units = "px",res = dpi)
VlnPlot(sc22_T, features = "percent.ribo") & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none") &
  labs(title="") 
dev.off()


#*****************************SUPPLEMENTARY FIGURE 6**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 6A

tiff(file='moma_featureplots_col3.tiff', width = 4500, height = 2500, units = "px",res = 300)
(((FeaturePlot(sc22_moma, features = "feat_AM1")+ labs(title="AM FES<br>(*APOE, CD163, MARCO, MSR1*)"))+
    (FeaturePlot(sc22_moma, features = "feat_eq.g2m.genes1")+ labs(title="G2M genes")))+
    (FeaturePlot(sc22_moma, features = "feat_HLA1")+ labs(title="MHC II FES<br>(*DQB, DQB.1, DQA, DQA.1, DRA, DRB*)"))+
    (FeaturePlot(sc22_moma, features = "feat_T1")+ labs(title="T cell FES<br>(*CD2, CD3D, CD3E, CD3G*)"))+
    (FeaturePlot(sc22_moma, features = "LOC100630729")+ labs(title="**Immunoglobulin Kappa<br>Variable 2-30-like*")) + 
    plot_layout(ncol = 3)) &
  theme(plot.title=element_markdown(size=16), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 6B

markers_all_moma <- readRDS("markers_all_moma.rds")

moma_top5_cod <- markers_all_moma %>% 
  group_by(cluster) %>%
  filter(!gene=="LOC111774319" & !gene=="LOC106781891" & !gene=="LOC111775035" & #removing non-coding genes
           !gene=="ND2" & #removing Mt gene
           !gene=="RPS27")%>% #removing RP genes
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_moma), value = TRUE)), n = 5)

moma_top5_cod$gene <- recode(moma_top5_cod$gene,"LOC100066849"="*CD33","LOC100054684"="*LILRB4","LOC100061154"="*MS4A6A",
                             "LOC100058587"="*H2AC20")
moma_markers <- c("*CD33","*LILRB4","*MS4A6A","*H2AC20")
moma_markers.loc <- c("LOC100066849","LOC100054684","LOC100061154","LOC100058587")

Idents(sc22_moma) <- "CellSubtype" 
sc22_moma@meta.data[,moma_markers]<-GetAssayData(sc22_moma)[moma_markers.loc,] %>% as.matrix %>% t

tiff(file='moma_dotplot.tiff', width = 4000, height = 1500, units = "px",res = 300)
DotPlot(sc22_moma, features = unique(moma_top5_cod$gene)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) 
dev.off()


#SUPPLEMENTARY FIGURE 6C

tiff(file='moma_vlnplot.tiff', width = dpi*12, height = dpi*5, units = "px",res = dpi)
(VlnPlot(sc22_moma, features = "nFeature_RNA") + labs(title="RNA feature count")) +
  (VlnPlot(sc22_moma, features = "percent.mito") + labs(title="Mitochondrial read %")) +
  plot_layout(ncol = 2) & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none") 
dev.off()


#SUPPLEMENTARY FIGURE 6D

tiff(file='moma_RP genes.tiff', width = dpi*12, height = dpi*5, units = "px",res = dpi)
VlnPlot(sc22_moma, features = "percent.ribo") & 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none") &
  labs(title="") 
dev.off()


#*****************************SUPPLEMENTARY FIGURE 7**************************************************
# In supplementary results

#SUPPLEMENTARY FIGURE 7A

tiff(file='DC_featureplots.tiff', width = 4500, height = 2500, units = "px",res = 300)
((FeaturePlot(sc22_DC, features = "feat_cDC21")+ labs(title="cDC2 FES<br>(**CLEC10A, CD1C, FCER1A*))"))+
    (FeaturePlot(sc22_DC, features = "feat_HLA1")+ labs(title="MHC II FES<br>(*DQB, DQB.1, DQA, DQA.1, DRA, DRB*)"))+
    (FeaturePlot(sc22_DC, features = "feat_act1")+ labs(title="Activated DC FES<br>(*CCR7, LAMP3, IDO1, CD83*)"))+
    (FeaturePlot(sc22_DC, features = "feat_AM1")+ labs(title="AM FES<br>(*APOE, CD163, MARCO, MSR1*)")) + 
    (FeaturePlot(sc22_DC, features = "feat_cDC11")+ labs(title="cDC1 FES<br>(*XCR1, CLEC9A, CADM1*)"))+
    plot_layout(ncol = 3)) &
  theme(plot.title=element_markdown(size=16), axis.title = element_text(size=14))
dev.off()


#SUPPLEMENTARY FIGURE 7B

sc22_DC<- SetIdent(sc22_DC, value = sc22_DC$CellSubtype) 
markers_all_DC <- readRDS("markers_all_DC.rds")
DC_top5 <- markers_all_DC %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_DC), value = TRUE)), n = 5)

DC_top5_cod <- markers_all_DC %>% 
  group_by(cluster) %>%
  filter(!gene=="LOC102149228" & !gene=="LOC102148672"& !gene=="LOC102148763")%>% #removing non-coding genes
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_DC), value = TRUE)), n = 5)

tiff(file='DC_dotplot.tiff', width = 4000, height = 1500, units = "px",res = 300)
DotPlot(sc22_DC, features = unique(DC_top5_cod$gene)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 14), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) 
dev.off()


#SUPPLEMENTARY FIGURE 7C

Idents(sc22_DC) <- "CellSubtype" 

tiff(file='DC_RNA.tiff', width = 1500, height = 1500, units = "px",res = 300)
VlnPlot(sc22_DC, features ="nFeature_RNA", pt.size = 0.1) + 
  labs(title = "RNA feature count")& 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position="none") &
  labs(title="")&
  scale_fill_brewer(palette="Oranges", direction = "-1")
dev.off()


#*****************************SUPPLEMENTARY FIGURE 8**************************************************
# In supplementary results

# Read the data file 
dcc_sc2022<-read.table("DCC_sc2022.txt",header=T,sep="\t")

# Create a subset of the data containing only the relevant columns
subset_data <- dcc_sc2022[, c("horse", "status", "type", "lymphocytes", "macrophages", "neutrophils", "mastocytes", "eosinophils")]

# Convert the data from wide to long format using tidyr::pivot_longer
long_data <- pivot_longer(subset_data, cols = c(lymphocytes, macrophages, neutrophils, mastocytes, eosinophils), names_to = "cell_type", values_to = "count")

# Create a new column with the modified horse ID
long_data$horse_id <- paste("Horse", long_data$horse)

# Create separate horizontally stacked bar plots for CTL horses
ctl_plot <- ggplot(long_data[long_data$status == "control", ], aes(x = type, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Fraction of cells (%)", fill = "Cell type") +
  scale_x_discrete(breaks=c("cytopost", "cytopre", "sc"),
                   labels=c("BALF cytology", "Cell suspension cytology", "scRNA-seq")) +
  scale_fill_manual(values = c("#FF0099", "#56B4E9","#F0E442", "#CC0000", "#CC99FF"), guide = "none" ) +
  theme_minimal() +
  facet_wrap(~ horse_id, scales = "free_y", ncol = 1, strip.position = "top") +
  coord_flip()

# Create separate horizontally stacked bar plots for SEA horses
sea_plot <- (ggplot(long_data[long_data$status == "sEA", ], aes(x = type, y = count, fill = cell_type)) +
               geom_bar(stat = "identity", position = "stack") +
               labs(x = "", y = "Fraction of cells (%)", fill = "Cell type") +
               scale_x_discrete(breaks=c("cytopost", "cytopre", "sc"),
                                labels=c("BALF cytology", "Cell suspension cytology", "scRNA-seq")) +
               scale_fill_manual(values = c("#FF0099", "#56B4E9", "#F0E442", "#CC0000", "#CC99FF"), guide = guide_legend(title = "Cell type")) +
               theme_minimal() +
               facet_wrap(~ horse_id, scales = "free_y", ncol = 1, strip.position = "top") +
               coord_flip()) 


# Arrange the plots into two columns using grid.arrange
dpi=300
tiff(file='DCC technique barplot.tiff', width = dpi*10, height = dpi*5, units = "px",res = dpi)
grid.arrange(ctl_plot + labs(title = "Control") + 
               theme(plot.title = element_text(hjust = 0.5, size = 14), 
                     axis.title.x = element_text(size=10), axis.text.y = element_text(size=9)),
             sea_plot + labs(title = "SEA") + 
               theme(plot.title = element_text(hjust = 0.5, size = 14), 
                     axis.title.x = element_text(size=10), axis.text.y = element_text(size=9)),
             ncol = 2, heights = c(0.95, 0.05))
dev.off()