
library(AnnotationHub)
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)



# read
covid_counts <- readRDS("~/Documents/COVID_data/data/GSM4557327_555_1_cell.counts.matrices.rds")
covid <- CreateSeuratObject(counts = covid_counts, project = "COVID", min.cells = 3, min.features = 200)
covid$group <- "COVID"
ctrl_counts <- readRDS("~/Documents/COVID_data/data/GSM4557334_HIP002_cell.counts.matrices.rds")
control <- CreateSeuratObject(counts = ctrl_counts, project = "Control", min.cells = 3, min.features = 200)
control$group <- "Control"
combined <- merge(covid, y = control, add.cell.ids = c("COVID", "Control"))

#preprocessing data
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:10)

# Cluster labels
DimPlot(combined, label = TRUE)

# COVID/Control
DimPlot(combined, group.by = "group")



#################################################
################################################
##ANNOTAION FOR CLUSTER###########################

mat_covid <- GetAssayData(combined, assay = "RNA", layer = "data.exon.COVID")
mat_ctrl  <- GetAssayData(combined, assay = "RNA", layer = "data.exon.Control")
common_genes <- intersect(rownames(mat_covid), rownames(mat_ctrl))
mat_covid <- mat_covid[common_genes, ]
mat_ctrl  <- mat_ctrl[common_genes, ]
data_mat <- cbind(mat_covid, mat_ctrl)
cell_names <- colnames(data_mat)
filtered_metadata <- combined@meta.data[cell_names, ]


combined_sce <- SingleCellExperiment(
  assays = list(logcounts = data_mat),
  colData = filtered_metadata
)


ref <- celldex::HumanPrimaryCellAtlasData()

cluster_preds <- SingleR(
  test = combined_sce,
  ref = ref,
  labels = ref$label.fine,
  clusters = combined$seurat_clusters[cell_names]
)


combined$ClusterAnnotation <- cluster_preds$labels[combined$seurat_clusters]


p <- DimPlot(combined, group.by = "ClusterAnnotation", label = TRUE, repel = TRUE)
#ggsave("UMAP_ClusterAnnotation.png", plot = p, width = 10, height = 8, dpi = 300)
print(p)


cluster_table <- data.frame(Cluster = rownames(cluster_preds), CellType = cluster_preds$labels)
################################################################################################
##################################COMPARE ANNOTATION WITH MANUALLY CURATED MARKERS##############
################################################################################################
marker_genes <- c(
  "CD8A"   = "CD8+ T cells",
  "CD8B"   = "CD8+ T cells",
  "CD14"   = "CD14+ Monocytes",
  "LYZ"    = "CD14+ Monocytes",
  "FCGR3A" = "CD16+ Monocytes",
  "CD34"   = "HSC",
  "GATA2"  = "HSC",
  "MS4A1"  = "Pre-B",
  "CD79A"  = "Pre-B",
  "GNLY"   = "NK cells",
  "NKG7"   = "NK cells"
)


plots <- lapply(names(marker_genes), function(gene) {
  celltype <- marker_genes[gene]
  FeaturePlot(combined, features = gene) + ggtitle(paste0(gene, " – ", celltype))
})

cluster_plot <- DimPlot(combined, group.by = "ClusterAnnotation") + ggtitle("Seurat Clusters")


final_plot <- wrap_plots(c(list(cluster_plot), plots), ncol = 4)


print(final_plot)
##################################################
############################## GO ################
##################################################


combined <- JoinLayers(combined)
markers <- FindMarkers(combined, ident.1 = 0, min.pct = 0.25)


sig_genes <- rownames(markers[markers$p_val_adj < 0.05, ])

entrez_ids <- bitr(sig_genes, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)


ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "MF",        # Biological Process
                pAdjustMethod = "BH")


dotplot(ego, showCategory = 15) + ggtitle("GO Enrichment – Cluster 0")


##################################################
############################## GSEA ################
##################################################



library(msigdbr)

markers <- FindMarkers(combined, ident.1 = 0, logfc.threshold = 0, min.pct = 0.1)
gene_list <- sort(markers$avg_log2FC, decreasing = TRUE)
names(gene_list) <- rownames(markers)


gsea_result <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "SYMBOL",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)



gseaplot2(gsea_result, geneSetID = gsea_result@result$ID[2],title = gsea_result@result$ID[2])


##################################################
############################## GSEA  H C2 C8 ################
##################################################



# C2: curated gene sets (e.g. KEGG, Reactome)
c2_genesets <- msigdbr(species = "Homo sapiens", category = "C2")

# C8: cell type gene sets
#c8_genesets <- msigdbr(species = "Homo sapiens", category = "C8")

# Hallmark: well-defined biological states
h_genesets <- msigdbr(species = "Homo sapiens", category = "H")



library(dplyr)

# Convert to list format: gene sets by term
c2_list <- c2_genesets[, c("gs_name", "gene_symbol")]
h_list <- h_genesets[, c("gs_name", "gene_symbol")]




gsea_hallmark <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = h_list,
  pvalueCutoff = 0.05
)


gsea_c2 <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = c2_list,
  pvalueCutoff = 0.05
)

dotplot(gsea_hallmark, showCategory = 15) + ggtitle("GSEA – Hallmark (Cluster 0)")
dotplot(gsea_c2, showCategory = 15) + ggtitle("GSEA – C2 (Cluster 0)")






res_df <- gsea_c2@result


sig_res <- res_df[res_df$qvalue < 0.05, ]


top_pathway <- sig_res[which.max(sig_res$NES), "ID"]
es_value <- sig_res[which.max(sig_res$NES), "enrichmentScore"]


getwd()

# 

