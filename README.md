# SCENICmsthesisscript
SCENIC code for transcriptional network analysis

# ============================================================
# 1. Regulon network inference
# ============================================================
# loading package
library(Seurat)
library(SeuratDisk)
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SingleCellExperiment)
# load main dataset
seu <- LoadH5Seurat(
  "/Users/sriharini/Downloads/Manz_cellbender_annotated_rpca_2000_30.h5seurat",
  assays = "RNA",                # load only active RNA data as Raw counts matrix; Normalized data matrix; Gene names; Cell barcodes
  reductions = c("pca", "umap"), # PCA, UMAP transforms a data set from a high-dimensional space into a low-dimensional space
  graphs = FALSE,
  neighbors = FALSE,
  images = FALSE,                # scenic needs expression matrix and metadata variables
  misc = FALSE,
  verbose = FALSE)
# subsetting on genotype_sex # female Cre/flox control
seu_FCre <- subset(seu, subset = genotype_sex == "F_Cre")
cat("Cells in F_Cre:", ncol(seu_FCre), "\n")
# find highly variable genes n=5000 genes; to reduce memory
seu_FCre <- FindVariableFeatures(seu_FCre, selection.method = "vst", nfeatures = 5000)
hvg <- VariableFeatures(seu_FCre)
# building expression matrix on subsetted FCre dataset
exprMat <- as.matrix(seu_FCre[["RNA"]]@counts[hvg, ]) # rows=genes, columns=cells
cellInfo1 <- seu_FCre@meta.data ## 
rm(seu) # 
# load files (cistarget database paths)
dbDir <- "~/Desktop/cisTarget_databases"
dbFiles <- c("mm9-500bp-upstream-7species.mc9nr.feather",   #TF binding sites near gene promoters
             "mm9-tss-centered-10kb-7species.mc9nr.feather") # 10kb checks broader regions around genes, rest excluded
motifAnnotFile <- file.path(dbDir, "motifs-v9-nr.mgi-m0.001-o0.0.tbl") # mgi "mus musculus"
motifAnnotations_mgi <- importAnnotations(motifAnnotFile) # motif IDs → TF gene names
# initialising SCENIC function; FCre
dir.create("int", showWarnings = FALSE) #output directory settings
scenicOptions <- initializeScenic(
  org = "mgi",                    # Mouse Gene Informatics
  dbDir = dbDir,                  # feather files
  dbs = dbFiles,                  # two database files
  datasetTitle = "F_Cre_HVG",       # analysis name
  nCores = 4)                    # 

saveRDS(scenicOptions, file = "int/scenicOptionsF_Cre.Rds")
### Co-expression network
# Load SCENIC options
scenicOptions <- readRDS("int/scenicOptionsF_Cre.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
# correlation analysis
runCorrelation(exprMat_filtered, scenicOptions)     #correlation method, select TF list
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)      #list of “regulatory links”:TF → target_gene
runGenie3(exprMat_filtered_log,                     # resume function for Genie3
          scenicOptions,
          resumePreviousRun = TRUE,
          nParts = 10)  
# Build and score the GRN
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)   
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)  # final RDS output file saved as 3.4_regulonAUC.Rds



# similar script has been ran for female B cell-specific IL-10 KO mice
# Data visualisation # libraries loaded for visualisation
library(ComplexHeatmap)  # heatmap
library(Seurat)          #PCA, UMAP handling
library(ggplot2)         # violin plot
library(patchwork)       # combine multiple plots

# Load SCENIC AUC output file from Cre/flox control and IL-10 KO mice for comparison analysis
auc_FCre <- getAUC(readRDS("~/Desktop/SCENIC_backup_FCre/3.4_regulonAUC.Rds")) 
auc_FKO  <- getAUC(readRDS("~/Desktop/SCENIC_backup_FKO/3.4_regulonAUC.Rds"))

info_FCre <- seu_FCre@meta.data  # grouping cells by celltype
info_FKO  <- seu_FKO@meta.data
# differential heatmap analysis
# Align regulons by TF name 
cleanTF <- function(x) trimws(gsub("\\s*\\(.*\\)|_extended","", x))   # trim extended name version and keeps only actual TF name 
tf_FCre <- cleanTF(rownames(auc_FCre))   # extract TFs from FCre
tf_FKO  <- cleanTF(rownames(auc_FKO))    # extract TFs from FKO
common  <- intersect(tf_FCre, tf_FKO)    # find common TFs

auc_FCre <- auc_FCre[match(common, tf_FCre), ]  #same TFs across rows in both matrices
auc_FKO  <- auc_FKO[match(common, tf_FKO), ]
rownames(auc_FCre) <- rownames(auc_FKO) <- common

# Cell types of interest
cts <- c(
  "Cycling HSC", 
  "HSC", 
  "Lymphoid T progenitor",
  "CD4 T-cell progenitor", 
  "CD8 T-cell NK cell progenitor",
  "B-cell progenitor", 
  "Early B-cell progenitor"
)
# calculate average auc values for each TF 
avgM <- function(auc, info)
  sapply(cts, \(ct) rowMeans(auc[, info$predicted.id==ct, drop=FALSE], na.rm=TRUE))

avg_FCre <- avgM(auc_FCre, info_FCre)
avg_FKO  <- avgM(auc_FKO,  info_FKO)

# Difference btw Cre and KO
diff <- avg_FKO - avg_FCre   # calculate difference in auc values for each tf between cre and KO
meanDiff <- rowMeans(diff)  #MEAN difference across all cell types

topTFs <- names(sort(abs(meanDiff), decreasing=TRUE))[1:40]  # gives top 40 differential TF

# ============================================================
# 2. Heatmap
# ============================================================
pdf("~/Downloads/SCENIC_FKO_vs_FCre_heatmap.pdf", width=8, height=10)

Heatmap(diff[topTFs, ], name="FKO - FCre",
        col=colorRamp2(c(-0.3,0,0.3), c("blue","white","red")),   # blue-lower tf value in KO; white-no change btw cre and ko; red-higher tf                                                                         value in ko compared to flox control
        cluster_rows=TRUE, cluster_columns=FALSE)                #plot top TF differences

dev.off()

# auc_FCre (aligned)
# auc_FKO (aligned)
# Both have common TF names as rownames
# Combine into ONE matrix
regulonAUC_combined <- cbind(auc_FCre, auc_FKO)   # rows-TFs; column-cell types
cat("Combined matrix dimensions:", dim(regulonAUC_combined), "\n")
# Create Seurat object with dummy counts
num_cells <- ncol(regulonAUC_combined)
fake_counts <- matrix(
  data = rpois(10 * num_cells, lambda = 5),
  nrow = 10,
  ncol = num_cells
)
rownames(fake_counts) <- paste0("gene", 1:10)
colnames(fake_counts) <- colnames(regulonAUC_combined)
seu <- CreateSeuratObject(counts = fake_counts)         # create empty seurat object
# Add SCENIC AUC as assay
seu[["RegulonAUC"]] <- CreateAssayObject(data = regulonAUC_combined)
DefaultAssay(seu) <- "RegulonAUC"    # store regulon matrix
#  Add metadata (condition and cell type)
# Add condition labels FCre or FKO
seu$Condition <- ifelse(
  colnames(seu) %in% colnames(auc_FCre),
  "F-Cre",
  "F-KO"
)
# Add cell types 
fcre_cells <- colnames(auc_FCre)
fko_cells  <- colnames(auc_FKO)
celltype_FCre <- info_FCre[fcre_cells, "predicted.id"]
celltype_FKO  <- info_FKO[fko_cells, "predicted.id"]

seu$CellType <- c(celltype_FCre, celltype_FKO) # each cell assigned to its cell type
# Scale
seu <- ScaleData(seu, assay = "RegulonAUC", verbose = FALSE)   # standardizes regulon activity
# Determine number of PCs
num_features <- nrow(seu[["RegulonAUC"]])
npcs_use <- min(15, num_features - 1)      # 40 TFs to 15 PCs
# Run PCA
seu <- RunPCA(
  seu,
  features = rownames(seu[["RegulonAUC"]]),
  assay = "RegulonAUC",
  npcs = npcs_use,
  verbose = FALSE
)
# Run UMAP
seu <- RunUMAP(seu, reduction = "pca", dims = 1:npcs_use, verbose = FALSE)  # cells with similar regulon activity cluster together
# ============================================================
# 3. Plot UMAPs
# ============================================================
# UMAP by condition
pdf("UMAP_FCre_vs_FKO_condition.pdf", width=6, height=5)
print(
  DimPlot(seu, reduction = "umap", group.by = "Condition", 
          cols = c("F-Cre" = "#0173B2", "F-KO" = "#DE8F05")) +
    ggtitle("UMAP: SCENIC Regulon Activity (F-Cre vs F-KO)")
)
dev.off()

# UMAP by cell type
pdf("UMAP_FCre_vs_FKO_celltype.pdf", width=12, height=5)
print(
  DimPlot(seu, reduction = "umap", group.by = "CellType", 
          label = FALSE, repel = TRUE) +
    ggtitle("UMAP:F-Cre vs F-KO ")
)
dev.off()

# Split by condition
pdf("UMAP_FCre_vs_FKO_split.pdf", width=14, height=5)
print(
  DimPlot(seu, reduction = "umap", group.by = "CellType", 
          split.by = "Condition", label = FALSE) +
    ggtitle("UMAP Split: F-Cre vs F-KO")
)
dev.off()

cat("✓ UMAP plots saved\n")
# Using top TFs from heatmap analysis
cat("\nTop 20 TFs from heatmap:\n")
print(topTFs[1:20])
# Add your TFs of interest
tfs_to_plot <- c("E2f1", "Myc", "E2f3", "Rad21", "Brca1", "Klf3")
# Plot each TF individually, split by condition
for (tf in tfs_to_plot) {
  p <- FeaturePlot(seu, features = tf, reduction = "umap", 
                   split.by = "Condition",
                   min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(paste0("Regulon Activity: ", tf, " (F-Cre vs F-KO)"))
  
  ggsave(paste0("UMAP_", tf, "_FCre_vs_FKO.pdf"), 
         plot = p, width = 10, height = 4.5, dpi = 300)
}

cat("✓ Individual TF plots saved\n")
# Plot multiple TFs in a grid
pdf("UMAPFemale_multiple_TFs_grid.pdf", width=15, height=10)
print(
  FeaturePlot(seu, features = tfs_to_plot[1:min(6, length(tfs_to_plot))], 
              reduction = "umap", ncol = 3,
              min.cutoff = "q05", max.cutoff = "q95")
)
dev.off()

cat("✓ Grid plot saved\n")
============================================================
# 4. Violin plots for top TFs
# ============================================================

# Create violin plots comparing F-Cre vs F-KO
violin_plots <- list()

for (i in 1:min(6, length(tfs_to_plot))) {
  tf <- tfs_to_plot[i]                                         #select TF
  seu[[paste0("TF_", tf)]] <- seu[["RegulonAUC"]]@data[tf, ]   # Add TF activity to metadata
  
  p <- ggplot(seu@meta.data, aes_string(
    x = "Condition",
    y = paste0("TF_", tf),
    fill = "Condition")) +                  # shows didtribution of regulon activity
    
    geom_violin(trim = FALSE, scale = "width") +       # density
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +      # show median and interquatiles
    
    scale_y_continuous(limits = c(0, 1)) +    # AUC values displayed from 0 to 1
    
    scale_fill_manual(values = c("F-Cre" = "#0173B2", "F-KO" = "#DE8F05")) +
    theme_classic() +
    ggtitle(tf) +
    ylab("AUC Score") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
    violin_plots[[i]] <- p
}

# Combine violin plots
pdf("Violin_top_TFs_FCre_vs_FKO.pdf", width=12, height=8)
print(wrap_plots(violin_plots, ncol = 3))
dev.off()

cat("✓ Violin plots saved\n")
