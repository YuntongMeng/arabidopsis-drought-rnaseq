# ============================================================
# WT drought vs watered RNA-seq analysis
# Arabidopsis thaliana
#
# Input:
#   counts/WT_6samples_raw_counts.tsv
#   scripts/sample_metadata.tsv
#
# Output:
#   PCA, sample-distance heatmap, MA plot, volcano plot,
#   top DEG heatmap, GO enrichment plots and result tables
# ============================================================


# ------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(scales)
library(dplyr)
library(AnnotationDbi)
library(org.At.tair.db)
library(clusterProfiler)


# ------------------------------------------------------------
# 2. Read count matrix and sample metadata
# ------------------------------------------------------------

counts <- read.delim(
  "counts/WT_6samples_raw_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

metadata_all <- read.delim(
  "scripts/sample_metadata.tsv"
)


# ------------------------------------------------------------
# 3. Prepare WT sample metadata
# ------------------------------------------------------------

# Keep only samples present in the WT count matrix
metadata_wt <- metadata_all[
  metadata_all$sample_id %in% colnames(counts),
]

# Match metadata order exactly to count matrix columns
metadata_wt <- metadata_wt[
  match(colnames(counts), metadata_wt$sample_id),
]

rownames(metadata_wt) <- metadata_wt$sample_id

# Confirm sample order is correct
stopifnot(
  all(colnames(counts) == rownames(metadata_wt))
)

# Set watered as the reference condition
metadata_wt$condition <- factor(
  metadata_wt$condition,
  levels = c("watered", "drought")
)


# ------------------------------------------------------------
# 4. Construct DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata_wt,
  design = ~ condition
)


# ------------------------------------------------------------
# 5. Prefilter low-count genes
# ------------------------------------------------------------

# Keep genes with at least 10 total reads across all six samples
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]


# ------------------------------------------------------------
# 6. Run DESeq2 differential expression analysis
# ------------------------------------------------------------

dds <- DESeq(dds)


# ------------------------------------------------------------
# 7. Variance-stabilizing transformation
# ------------------------------------------------------------

# VST is used for exploratory analyses such as PCA and clustering.
# Differential expression testing itself is still performed on raw counts.
vsd <- vst(
  dds,
  blind = FALSE
)


# ------------------------------------------------------------
# 8. PCA
# ------------------------------------------------------------

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
) +
  ggtitle("PCA of WT drought and watered samples")

ggsave(
  "results/WT_PCA.png",
  plot = pca_plot,
  width = 7,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 9. Sample-distance heatmap
# ------------------------------------------------------------

sampleDists <- dist(
  t(assay(vsd))
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

pheatmap(
  sampleDistMatrix,
  filename = "results/WT_sample_distance_heatmap.png",
  width = 8,
  height = 7
)


# ------------------------------------------------------------
# 10. Differential expression results
# ------------------------------------------------------------

# Contrast: drought vs watered
# Positive log2FC = higher expression under drought
# Negative log2FC = lower expression under drought
res <- results(
  dds,
  contrast = c(
    "condition",
    "drought",
    "watered"
  )
)

res_ordered <- res[
  order(res$padj),
]


# ------------------------------------------------------------
# 11. Convert results to data frame
# ------------------------------------------------------------

res_df <- as.data.frame(res)

res_df$TAIR <- rownames(
  res_df
)


# ------------------------------------------------------------
# 12. Annotate TAIR gene IDs
# ------------------------------------------------------------

anno <- AnnotationDbi::select(
  org.At.tair.db,
  keys = res_df$TAIR,
  keytype = "TAIR",
  columns = c(
    "SYMBOL",
    "GENENAME"
  )
)

# Collapse one-to-many annotations while retaining all aliases
anno_collapsed <- anno %>%
  group_by(TAIR) %>%
  summarise(
    SYMBOL = paste(
      unique(na.omit(SYMBOL)),
      collapse = "; "
    ),
    GENENAME = paste(
      unique(na.omit(GENENAME)),
      collapse = "; "
    ),
    .groups = "drop"
  )

# Replace empty strings with NA
anno_collapsed$SYMBOL[
  anno_collapsed$SYMBOL == ""
] <- NA

anno_collapsed$GENENAME[
  anno_collapsed$GENENAME == ""
] <- NA

# Join annotations to DESeq2 results
res_annotated <- res_df %>%
  left_join(
    anno_collapsed,
    by = "TAIR"
  )


# ------------------------------------------------------------
# 13. Define differentially expressed genes
# ------------------------------------------------------------

# Main DEG criterion:
# adjusted p-value < 0.05
# absolute log2 fold change >= 1
deg_annotated <- res_annotated %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  ) %>%
  arrange(padj)


# ------------------------------------------------------------
# 14. MA plot
# ------------------------------------------------------------

ma_df <- as.data.frame(res)

ma_df$gene_id <- rownames(
  ma_df
)

ma_df$status <- "Not significant"

ma_df$status[
  !is.na(ma_df$padj) &
    ma_df$padj < 0.05 &
    ma_df$log2FoldChange >= 1
] <- "Up"

ma_df$status[
  !is.na(ma_df$padj) &
    ma_df$padj < 0.05 &
    ma_df$log2FoldChange <= -1
] <- "Down"

ma_df$status <- factor(
  ma_df$status,
  levels = c(
    "Down",
    "Not significant",
    "Up"
  )
)

n_up <- sum(
  ma_df$status == "Up",
  na.rm = TRUE
)

n_down <- sum(
  ma_df$status == "Down",
  na.rm = TRUE
)

p_ma <- ggplot(
  ma_df,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.2
  ) +
  geom_hline(
    yintercept = 0
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  scale_x_log10(
    labels = comma
  ) +
  scale_color_manual(
    values = c(
      "Down" = "#3B82F6",
      "Not significant" = "grey75",
      "Up" = "#EF4444"
    ),
    labels = c(
      paste0(
        "Down (",
        n_down,
        ")"
      ),
      "Not significant",
      paste0(
        "Up (",
        n_up,
        ")"
      )
    )
  ) +
  coord_cartesian(
    ylim = c(-8, 8)
  ) +
  labs(
    title = "WT drought response",
    subtitle = "DESeq2 differential expression: drought vs watered",
    x = "Mean normalized count",
    y = "log2 fold change",
    color = NULL
  ) +
  theme_classic()

ggsave(
  "results/WT_MA_plot.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 15. Volcano plot
# ------------------------------------------------------------

volcano_df <- as.data.frame(res)

volcano_df$gene_id <- rownames(
  volcano_df
)

volcano_df$status <- "Not significant"

volcano_df$status[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange >= 1
] <- "Up"

volcano_df$status[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange <= -1
] <- "Down"

volcano_df$status <- factor(
  volcano_df$status,
  levels = c(
    "Down",
    "Not significant",
    "Up"
  )
)

volcano_df$minus_log10_padj <- -log10(
  volcano_df$padj
)

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Down" = "#3B82F6",
      "Not significant" = "grey75",
      "Up" = "#EF4444"
    ),
    labels = c(
      paste0(
        "Down (",
        n_down,
        ")"
      ),
      "Not significant",
      paste0(
        "Up (",
        n_up,
        ")"
      )
    )
  ) +
  labs(
    title = "WT drought response",
    subtitle = "DESeq2 differential expression: drought vs watered",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_classic()

ggsave(
  "results/WT_volcano_plot.png",
  plot = p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 16. Top-30 DEG heatmap
# ------------------------------------------------------------

# Select the 30 DEGs with the smallest adjusted p-values
top30 <- head(
  deg_annotated$TAIR,
  30
)

# Extract VST expression values
mat_top30 <- assay(vsd)[
  top30,
]

# Standardize each gene across samples
# Each row is converted to a z-score
mat_top30_z <- t(
  scale(
    t(mat_top30)
  )
)

# Prepare column annotation
annotation_col <- data.frame(
  condition = colData(vsd)$condition
)

rownames(annotation_col) <- colnames(
  vsd
)

# Create readable row labels
# Use the first available gene symbol; fall back to TAIR ID
label_df <- deg_annotated %>%
  filter(
    TAIR %in% top30
  ) %>%
  dplyr::select(
    TAIR,
    SYMBOL
  )

label_df$display_name <- ifelse(
  is.na(label_df$SYMBOL),
  label_df$TAIR,
  sub(
    ";.*$",
    "",
    label_df$SYMBOL
  )
)

row_labels <- label_df$display_name[
  match(
    rownames(mat_top30_z),
    label_df$TAIR
  )
]

pheatmap(
  mat_top30_z,
  annotation_col = annotation_col,
  labels_row = row_labels,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  filename = "results/WT_top30_DEG_heatmap.png",
  width = 8,
  height = 9
)


# ------------------------------------------------------------
# 17. Prepare upregulated and downregulated genes
# ------------------------------------------------------------

up_genes <- deg_annotated %>%
  filter(
    log2FoldChange > 0
  ) %>%
  pull(TAIR)

down_genes <- deg_annotated %>%
  filter(
    log2FoldChange < 0
  ) %>%
  pull(TAIR)


# ------------------------------------------------------------
# 18. Define GO enrichment background
# ------------------------------------------------------------

# Use genes that were actually tested by DESeq2 and have valid p-values.
# This is a more appropriate enrichment background than the entire
# Arabidopsis annotation database.
go_universe <- res_annotated %>%
  filter(
    !is.na(pvalue)
  ) %>%
  pull(TAIR)


# ------------------------------------------------------------
# 19. GO Biological Process enrichment
# ------------------------------------------------------------

# Upregulated genes
ego_up <- enrichGO(
  gene = up_genes,
  universe = go_universe,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Downregulated genes
ego_down <- enrichGO(
  gene = down_genes,
  universe = go_universe,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# ------------------------------------------------------------
# 20. GO enrichment dotplots
# ------------------------------------------------------------

p_up <- dotplot(
  ego_up,
  showCategory = 10,
  title = "GO enrichment of upregulated genes"
) +
  labs(
    color = "Adjusted p-value"
  )

ggsave(
  "results/WT_GO_enrichment_upregulated.png",
  plot = p_up,
  width = 8,
  height = 6,
  dpi = 300
)

p_down <- dotplot(
  ego_down,
  showCategory = 10,
  title = "GO enrichment of downregulated genes"
) +
  labs(
    color = "Adjusted p-value"
  )

ggsave(
  "results/WT_GO_enrichment_downregulated.png",
  plot = p_down,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 21. Export result tables
# ------------------------------------------------------------

# Complete DESeq2 results with annotations
write.table(
  res_annotated,
  file = "results/WT_DESeq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Significant DEGs
write.table(
  deg_annotated,
  file = "results/WT_DEGs_padj0.05_log2FC1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# GO enrichment results
write.table(
  as.data.frame(ego_up),
  file = "results/WT_GO_enrichment_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(ego_down),
  file = "results/WT_GO_enrichment_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Final summary
# ------------------------------------------------------------

cat(
  "\nWT drought vs watered analysis complete\n"
)

cat(
  "Genes retained after prefilter:",
  nrow(dds),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(deg_annotated),
  "\n"
)

cat(
  "Upregulated:",
  sum(
    deg_annotated$log2FoldChange > 0
  ),
  "\n"
)

cat(
  "Downregulated:",
  sum(
    deg_annotated$log2FoldChange < 0
  ),
  "\n"
)