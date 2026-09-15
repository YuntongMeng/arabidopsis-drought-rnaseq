# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Full factorial DESeq2 analysis: genotype × condition
#
# Samples:
#   WT watered      n = 3
#   WT drought      n = 3
#   BRL3 watered    n = 2
#   BRL3 drought    n = 3
#
# Design:
#   ~ genotype * condition
#
# Main question:
#   Does the transcriptional response to drought differ
#   between BRL3 and WT?
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(ggplot2)
library(scales)
library(pheatmap)
library(org.At.tair.db)
library(AnnotationDbi)
library(clusterProfiler)


# ------------------------------------------------------------
# 2. Load count matrix and sample metadata
# ------------------------------------------------------------

counts_11 <- read.delim(
  "counts/all_11samples_raw_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

metadata_all <- read.delim(
  "scripts/sample_metadata.tsv"
)


# Reorder metadata to match count-matrix columns
metadata_11 <- metadata_all[
  match(colnames(counts_11), metadata_all$sample_id),
]

rownames(metadata_11) <- metadata_11$sample_id


# Confirm that sample order is identical
stopifnot(
  all(colnames(counts_11) == rownames(metadata_11))
)


# ------------------------------------------------------------
# 3. Define experimental factors
# ------------------------------------------------------------

# WT and watered are set as reference levels.
#
# Therefore:
#   genotype_BRL3_vs_WT
#     = BRL3 vs WT under watered conditions
#
#   condition_drought_vs_watered
#     = drought vs watered in WT
#
#   genotypeBRL3.conditiondrought
#     = genotype × drought interaction

metadata_11$genotype <- factor(
  metadata_11$genotype,
  levels = c("WT", "BRL3")
)

metadata_11$condition <- factor(
  metadata_11$condition,
  levels = c("watered", "drought")
)


# ------------------------------------------------------------
# 4. Create DESeq2 dataset
# ------------------------------------------------------------

dds_11 <- DESeqDataSetFromMatrix(
  countData = counts_11,
  colData = metadata_11,
  design = ~ genotype * condition
)


# ------------------------------------------------------------
# 5. Prefilter low-count genes
# ------------------------------------------------------------

# Keep genes with at least 10 total reads across all 11 samples.
# This removes genes with extremely little information before
# model fitting.

keep_11 <- rowSums(counts(dds_11)) >= 10

dds_11 <- dds_11[keep_11,]


# ------------------------------------------------------------
# 6. Fit DESeq2 model
# ------------------------------------------------------------

dds_11 <- DESeq(dds_11)


# Inspect model coefficients
resultsNames(dds_11)


# ------------------------------------------------------------
# 7. Variance-stabilizing transformation
# ------------------------------------------------------------

# VST is used for visualization and sample-level exploration.
# Differential-expression testing is still performed using
# the negative-binomial model on count data.

vsd_11 <- vst(
  dds_11,
  blind = FALSE
)


# ------------------------------------------------------------
# 8. PCA of all 11 samples
# ------------------------------------------------------------

pca_11 <- plotPCA(
  vsd_11,
  intgroup = c("genotype", "condition")
) +
  ggtitle("PCA of WT and BRL3 samples")


ggsave(
  "results/FULL_PCA.png",
  plot = pca_11,
  width = 7,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 9. Sample-distance heatmap
# ------------------------------------------------------------

sampleDists_11 <- dist(
  t(assay(vsd_11))
)

sampleDistMatrix_11 <- as.matrix(
  sampleDists_11
)

rownames(sampleDistMatrix_11) <- colnames(vsd_11)
colnames(sampleDistMatrix_11) <- colnames(vsd_11)


annotation_col_11 <- data.frame(
  genotype = colData(vsd_11)$genotype,
  condition = colData(vsd_11)$condition
)

rownames(annotation_col_11) <- colnames(vsd_11)


pheatmap(
  sampleDistMatrix_11,
  annotation_col = annotation_col_11,
  filename = "results/FULL_sample_distance_heatmap.png",
  width = 8,
  height = 7
)


# ============================================================
# 10. Genotype × drought interaction analysis
# ============================================================

# Interaction coefficient:
#
# (BRL3 drought - BRL3 watered)
# -
# (WT drought - WT watered)
#
# Positive interaction:
#   drought response is more positive in BRL3 than in WT
#
# Negative interaction:
#   drought response is more negative / less positive
#   in BRL3 than in WT


res_interaction <- results(
  dds_11,
  name = "genotypeBRL3.conditiondrought"
)

summary(res_interaction)


res_interaction_df <- as.data.frame(
  res_interaction
)


# ------------------------------------------------------------
# 11. Define significant interaction genes
# ------------------------------------------------------------

# Same final DEG thresholds as the Stage 1 WT analysis:
#   adjusted p-value < 0.05
#   absolute log2 fold change >= 1

interaction_deg <- res_interaction_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )


# ------------------------------------------------------------
# 12. Add Arabidopsis gene annotation
# ------------------------------------------------------------

res_interaction_df$TAIR <- rownames(
  res_interaction_df
)


anno_interaction <- AnnotationDbi::select(
  org.At.tair.db,
  keys = res_interaction_df$TAIR,
  keytype = "TAIR",
  columns = c(
    "SYMBOL",
    "GENENAME"
  )
)


# Some TAIR IDs map to multiple aliases.
# Collapse multiple mappings into one row per TAIR gene.

anno_interaction_collapsed <- anno_interaction %>%
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


# Replace empty annotation strings with NA
anno_interaction_collapsed$SYMBOL[
  anno_interaction_collapsed$SYMBOL == ""
] <- NA

anno_interaction_collapsed$GENENAME[
  anno_interaction_collapsed$GENENAME == ""
] <- NA


res_interaction_annotated <- res_interaction_df %>%
  left_join(
    anno_interaction_collapsed,
    by = "TAIR"
  )


interaction_deg_annotated <- res_interaction_annotated %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  ) %>%
  arrange(padj)


# Export complete and significant interaction results

write.table(
  res_interaction_annotated,
  file = "results/FULL_interaction_DESeq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  interaction_deg_annotated,
  file = "results/FULL_interaction_DEGs_padj0.05_log2FC1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Interaction MA plot
# ------------------------------------------------------------

ma_interaction_df <- as.data.frame(
  res_interaction
)


ma_interaction_df$status <- "Not significant"

ma_interaction_df$status[
  !is.na(ma_interaction_df$padj) &
    ma_interaction_df$padj < 0.05 &
    ma_interaction_df$log2FoldChange >= 1
] <- "Positive interaction"

ma_interaction_df$status[
  !is.na(ma_interaction_df$padj) &
    ma_interaction_df$padj < 0.05 &
    ma_interaction_df$log2FoldChange <= -1
] <- "Negative interaction"


ma_interaction_df$status <- factor(
  ma_interaction_df$status,
  levels = c(
    "Negative interaction",
    "Not significant",
    "Positive interaction"
  )
)


p_ma_interaction <- ggplot(
  ma_interaction_df,
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
      "Negative interaction" = "#3B82F6",
      "Not significant" = "grey75",
      "Positive interaction" = "#EF4444"
    )
  ) +
  labs(
    title = "Genotype × drought interaction",
    subtitle = "BRL3 drought response relative to WT drought response",
    x = "Mean normalized count",
    y = "Interaction log2 fold change",
    color = NULL
  ) +
  theme_classic()


ggsave(
  "results/FULL_interaction_MA_plot.png",
  plot = p_ma_interaction,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 14. Interaction volcano plot
# ------------------------------------------------------------

volcano_interaction_df <- as.data.frame(
  res_interaction
)


volcano_interaction_df$status <- "Not significant"

volcano_interaction_df$status[
  !is.na(volcano_interaction_df$padj) &
    volcano_interaction_df$padj < 0.05 &
    volcano_interaction_df$log2FoldChange >= 1
] <- "Positive interaction"

volcano_interaction_df$status[
  !is.na(volcano_interaction_df$padj) &
    volcano_interaction_df$padj < 0.05 &
    volcano_interaction_df$log2FoldChange <= -1
] <- "Negative interaction"


volcano_interaction_df$status <- factor(
  volcano_interaction_df$status,
  levels = c(
    "Negative interaction",
    "Not significant",
    "Positive interaction"
  )
)


volcano_interaction_df$minus_log10_padj <- -log10(
  volcano_interaction_df$padj
)


p_volcano_interaction <- ggplot(
  volcano_interaction_df,
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
      "Negative interaction" = "#3B82F6",
      "Not significant" = "grey75",
      "Positive interaction" = "#EF4444"
    )
  ) +
  labs(
    title = "Genotype × drought interaction",
    subtitle = "BRL3 drought response relative to WT drought response",
    x = "Interaction log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_classic()


ggsave(
  "results/FULL_interaction_volcano_plot.png",
  plot = p_volcano_interaction,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 15. Top 30 interaction-gene heatmap
# ------------------------------------------------------------

# Select the 30 most statistically significant interaction genes.

top30_interaction <- head(
  interaction_deg_annotated$TAIR,
  30
)


mat_interaction_top30 <- assay(vsd_11)[
  top30_interaction,
]


# Row-wise z-score:
# standardize each gene across the 11 samples.
#
# Therefore heatmap colors represent expression relative
# to each gene's own mean, not absolute expression between genes.

mat_interaction_top30_z <- t(
  scale(
    t(mat_interaction_top30)
  )
)


# Create readable gene labels.
# Use the first gene symbol when available;
# otherwise retain the TAIR identifier.

label_interaction_df <- interaction_deg_annotated %>%
  filter(
    TAIR %in% top30_interaction
  ) %>%
  dplyr::select(
    TAIR,
    SYMBOL
  )


label_interaction_df$display_name <- ifelse(
  is.na(label_interaction_df$SYMBOL),
  label_interaction_df$TAIR,
  sub(
    ";.*$",
    "",
    label_interaction_df$SYMBOL
  )
)


row_labels_interaction <- label_interaction_df$display_name[
  match(
    rownames(mat_interaction_top30_z),
    label_interaction_df$TAIR
  )
]


pheatmap(
  mat_interaction_top30_z,
  annotation_col = annotation_col_11,
  labels_row = row_labels_interaction,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  filename = "results/FULL_interaction_top30_DEG_heatmap.png",
  width = 8,
  height = 9
)


# ============================================================
# 16. GO enrichment of interaction genes
# ============================================================

# Positive and negative interaction genes are analyzed separately
# because the direction of the interaction has different
# biological interpretations.


positive_interaction_genes <- interaction_deg_annotated %>%
  filter(
    log2FoldChange > 0
  ) %>%
  pull(TAIR)


negative_interaction_genes <- interaction_deg_annotated %>%
  filter(
    log2FoldChange < 0
  ) %>%
  pull(TAIR)


# Use genes that received a valid DESeq2 p-value as the
# enrichment background rather than the entire genome.

interaction_universe <- res_interaction_annotated %>%
  filter(
    !is.na(pvalue)
  ) %>%
  pull(TAIR)


# ------------------------------------------------------------
# 17. Positive-interaction GO enrichment
# ------------------------------------------------------------

ego_positive_interaction <- enrichGO(
  gene = positive_interaction_genes,
  universe = interaction_universe,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


write.table(
  as.data.frame(ego_positive_interaction),
  file = "results/FULL_interaction_GO_enrichment_positive.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Positive interaction genes produced significant GO terms,
# so a dot plot is generated.

p_go_positive_interaction <- enrichplot::dotplot(
  ego_positive_interaction,
  showCategory = 10,
  title = "GO enrichment of positive interaction genes"
)


ggsave(
  "results/FULL_interaction_GO_enrichment_positive.png",
  plot = p_go_positive_interaction,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 18. Negative-interaction GO enrichment
# ------------------------------------------------------------

ego_negative_interaction <- enrichGO(
  gene = negative_interaction_genes,
  universe = interaction_universe,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# No significant GO BP terms were detected under the current
# thresholds. The empty result table is still exported to
# document that the analysis was performed.

write.table(
  as.data.frame(ego_negative_interaction),
  file = "results/FULL_interaction_GO_enrichment_negative.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 19. Additional factorial-model contrasts
# ============================================================

# These contrasts provide context for the interaction analysis.
# No additional figures are generated because the interaction
# remains the primary biological question.


# ------------------------------------------------------------
# 19a. BRL3 drought vs watered
# ------------------------------------------------------------

# BRL3 drought response =
# WT drought main effect + genotype:drought interaction

res_BRL3_drought <- results(
  dds_11,
  contrast = list(
    c(
      "condition_drought_vs_watered",
      "genotypeBRL3.conditiondrought"
    )
  )
)


res_BRL3_drought_df <- as.data.frame(
  res_BRL3_drought
)


BRL3_drought_DEG <- res_BRL3_drought_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )


# ------------------------------------------------------------
# 19b. BRL3 vs WT under watered conditions
# ------------------------------------------------------------

# Because watered is the reference condition,
# this is directly represented by the genotype main effect.

res_BRL3_vs_WT_watered <- results(
  dds_11,
  name = "genotype_BRL3_vs_WT"
)


res_BRL3_vs_WT_watered_df <- as.data.frame(
  res_BRL3_vs_WT_watered
)


BRL3_vs_WT_watered_DEG <- res_BRL3_vs_WT_watered_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )


# ------------------------------------------------------------
# 19c. BRL3 vs WT under drought conditions
# ------------------------------------------------------------

# Genotype difference under drought =
# genotype main effect + genotype:drought interaction

res_BRL3_vs_WT_drought <- results(
  dds_11,
  contrast = list(
    c(
      "genotype_BRL3_vs_WT",
      "genotypeBRL3.conditiondrought"
    )
  )
)


res_BRL3_vs_WT_drought_df <- as.data.frame(
  res_BRL3_vs_WT_drought
)


BRL3_vs_WT_drought_DEG <- res_BRL3_vs_WT_drought_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )


# ============================================================
# 20. Summary of the five core comparisons
# ============================================================

# WT drought vs watered from the full factorial model
res_WT_drought_full <- results(
  dds_11,
  name = "condition_drought_vs_watered"
)

res_WT_drought_full_df <- as.data.frame(res_WT_drought_full)

WT_drought_full_DEG <- res_WT_drought_full_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )
comparison_summary <- data.frame(
  comparison = c(
    "WT drought vs watered",
    "BRL3 drought vs watered",
    "BRL3 vs WT under watered",
    "BRL3 vs WT under drought",
    "Genotype x drought interaction"
  ),
  
  total_DEGs = c(
    nrow(WT_drought_full_DEG),
    nrow(BRL3_drought_DEG),
    nrow(BRL3_vs_WT_watered_DEG),
    nrow(BRL3_vs_WT_drought_DEG),
    nrow(interaction_deg_annotated)
  ),
  
  positive_or_up = c(
    sum(WT_drought_full_DEG$log2FoldChange > 0),
    sum(BRL3_drought_DEG$log2FoldChange > 0),
    sum(BRL3_vs_WT_watered_DEG$log2FoldChange > 0),
    sum(BRL3_vs_WT_drought_DEG$log2FoldChange > 0),
    sum(interaction_deg_annotated$log2FoldChange > 0)
  ),
  
  negative_or_down = c(
    sum(WT_drought_full_DEG$log2FoldChange < 0),
    sum(BRL3_drought_DEG$log2FoldChange < 0),
    sum(BRL3_vs_WT_watered_DEG$log2FoldChange < 0),
    sum(BRL3_vs_WT_drought_DEG$log2FoldChange < 0),
    sum(interaction_deg_annotated$log2FoldChange < 0)
  )
)


comparison_summary


write.table(
  comparison_summary,
  file = "results/FULL_comparison_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# End of full factorial DESeq2 analysis
# ============================================================