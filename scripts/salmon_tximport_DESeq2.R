# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Salmon -> tximport -> DESeq2, and comparison with the
# HISAT2 + featureCounts workflow
#
# Question:
#   Do the conclusions of the factorial analysis depend on the
#   upstream method?
#
#   Workflow A (existing): HISAT2 genome alignment
#                          -> featureCounts gene counts
#   Workflow B (this one): Salmon transcript quantification
#                          -> tximport gene-level summary
#
#   Both are analysed with the same DESeq2 model
#   (~ genotype * condition), the same prefilter and the same
#   DEG thresholds, so differences come from the upstream step.
#
# Input:
#   salmon_quant/<SRR>/quant.sf           (11 samples)
#   salmon_quant/tx2gene_TAIR10.63.tsv    (from the TAIR10.63 GTF)
#   salmon_quant/salmon_mapping_summary.tsv
#   scripts/sample_metadata.tsv
#   counts/all_11samples_featureCounts.txt.summary
#   results/FULL_interaction_DESeq2_all_results.tsv
#   results/FULL_WT_drought_DESeq2_all_results.tsv
#   results/FULL_BRL3_drought_DESeq2_all_results.tsv
#
# Output (all new files):
#   counts/salmon_11samples_gene_counts.tsv
#   results/salmon_gene_TPM.tsv
#   results/salmon_<contrast>_DESeq2_all_results.tsv
#   results/salmon_comparison_summary.tsv
#   results/comparison_log2FC_correlation.tsv
#   results/comparison_log2FC_scatter.png
#   results/comparison_DEG_overlap.tsv
#   results/comparison_DEG_overlap_upset.png
#   results/comparison_mapping_rates.tsv
#   results/comparison_pipeline_specific_DEGs.tsv
#   results/comparison_pipeline_specific_summary.tsv
#
# Run from the project root (the folder with the .Rproj file).
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(org.At.tair.db)
library(AnnotationDbi)
library(UpSetR)


# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

results_dir <- "results"
counts_dir <- "counts"
salmon_dir <- "salmon_quant"

# Same DEG thresholds as the existing analysis.
padj_cutoff <- 0.05
lfc_cutoff <- 1


# ------------------------------------------------------------
# 3. Sample metadata and Salmon files
# ------------------------------------------------------------

metadata_11 <- read.delim(
  "scripts/sample_metadata.tsv"
)

rownames(metadata_11) <- metadata_11$sample_id

# WT and watered are the reference levels, exactly as in
# scripts/FULL_factorial_DESeq2_analysis.R.

metadata_11$genotype <- factor(
  metadata_11$genotype,
  levels = c("WT", "BRL3")
)

metadata_11$condition <- factor(
  metadata_11$condition,
  levels = c("watered", "drought")
)

quant_files <- file.path(
  salmon_dir,
  metadata_11$sra_run,
  "quant.sf"
)

names(quant_files) <- metadata_11$sample_id

stopifnot(
  all(file.exists(quant_files))
)


# ------------------------------------------------------------
# 4. Transcript-to-gene table
# ------------------------------------------------------------

# Written by scripts/upstream_salmon_workflow.sh from
# Arabidopsis_thaliana.TAIR10.63.gtf, the same annotation that
# featureCounts used. That script has already checked that the
# cDNA + ncRNA FASTA and the GTF contain identical transcripts.

tx2gene <- read.delim(
  file.path(salmon_dir, "tx2gene_TAIR10.63.tsv"),
  header = FALSE,
  col.names = c("transcript_id", "gene_id")
)


# Version check on the R side: every transcript quantified by
# Salmon must be in tx2gene, and the resulting gene IDs must be
# the gene IDs of the featureCounts matrix.

quant_transcripts <- read.delim(
  quant_files[1]
)$Name

stopifnot(
  all(quant_transcripts %in% tx2gene$transcript_id),
  all(tx2gene$transcript_id %in% quant_transcripts)
)

featurecounts_genes <- read.delim(
  file.path(counts_dir, "all_11samples_raw_counts.tsv"),
  check.names = FALSE
)$gene_id

stopifnot(
  setequal(
    unique(tx2gene$gene_id),
    featurecounts_genes
  )
)

cat(
  "Annotation check passed:",
  nrow(tx2gene), "transcripts,",
  length(unique(tx2gene$gene_id)), "genes,",
  "identical gene set to the featureCounts matrix.\n"
)


# ------------------------------------------------------------
# 5. Import Salmon results and summarize to gene level
# ------------------------------------------------------------

# tximport adds up the transcripts of each gene and returns:
#   counts     estimated number of reads per gene
#   abundance  TPM per gene
#   length     average transcript length per gene and sample,
#              weighted by transcript abundance
#
# countsFromAbundance = "no" (the default): counts stay on the
# original scale, and the gene lengths are passed to DESeq2 as
# an offset. This corrects for changes in gene length between
# samples caused by different isoform usage.

txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)

stopifnot(
  all(colnames(txi$counts) == rownames(metadata_11))
)


# ------------------------------------------------------------
# 6. Export gene-level counts and TPM
# ------------------------------------------------------------

# Estimated counts are not integers, because reads compatible
# with several transcripts are divided between them.

salmon_counts_out <- data.frame(
  gene_id = rownames(txi$counts),
  round(txi$counts, 3),
  check.names = FALSE
)

write.table(
  salmon_counts_out,
  file = file.path(counts_dir, "salmon_11samples_gene_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# TPM (transcripts per million), gene level.
#
# IMPORTANT: this table is for display and for comparing the
# expression level of different genes within a sample. It is
# NOT used for differential expression testing. TPM is a
# relative measure that sums to one million in every sample, so
# a strong change in a few highly expressed genes shifts the
# TPM of all other genes. Differential expression is tested by
# DESeq2 on the estimated counts (section 7).

salmon_tpm_out <- data.frame(
  gene_id = rownames(txi$abundance),
  round(txi$abundance, 3),
  check.names = FALSE
)

write.table(
  salmon_tpm_out,
  file = file.path(results_dir, "salmon_gene_TPM.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. DESeq2 with the same factorial model
# ------------------------------------------------------------

dds_salmon <- DESeqDataSetFromTximport(
  txi,
  colData = metadata_11,
  design = ~ genotype * condition
)

# Same prefilter as the featureCounts analysis.
keep_salmon <- rowSums(counts(dds_salmon)) >= 10

dds_salmon <- dds_salmon[keep_salmon, ]

dds_salmon <- DESeq(dds_salmon)

resultsNames(dds_salmon)

stopifnot(
  identical(
    resultsNames(dds_salmon),
    c(
      "Intercept",
      "genotype_BRL3_vs_WT",
      "condition_drought_vs_watered",
      "genotypeBRL3.conditiondrought"
    )
  )
)


# ------------------------------------------------------------
# 8. Extract the five contrasts
# ------------------------------------------------------------

# Defined exactly as in scripts/FULL_factorial_DESeq2_analysis.R.

res_salmon <- list(
  WT_drought = results(
    dds_salmon,
    name = "condition_drought_vs_watered"
  ),
  BRL3_drought = results(
    dds_salmon,
    contrast = list(
      c(
        "condition_drought_vs_watered",
        "genotypeBRL3.conditiondrought"
      )
    )
  ),
  BRL3_vs_WT_watered = results(
    dds_salmon,
    name = "genotype_BRL3_vs_WT"
  ),
  BRL3_vs_WT_drought = results(
    dds_salmon,
    contrast = list(
      c(
        "genotype_BRL3_vs_WT",
        "genotypeBRL3.conditiondrought"
      )
    )
  ),
  interaction = results(
    dds_salmon,
    name = "genotypeBRL3.conditiondrought"
  )
)


# ------------------------------------------------------------
# 9. Annotate and export complete Salmon result tables
# ------------------------------------------------------------

anno <- AnnotationDbi::select(
  org.At.tair.db,
  keys = rownames(dds_salmon),
  keytype = "TAIR",
  columns = c(
    "SYMBOL",
    "GENENAME"
  )
)

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

anno_collapsed$SYMBOL[
  anno_collapsed$SYMBOL == ""
] <- NA

anno_collapsed$GENENAME[
  anno_collapsed$GENENAME == ""
] <- NA


annotate_results <- function(res) {
  res_df <- as.data.frame(res)

  res_df$TAIR <- rownames(res_df)

  res_df %>%
    left_join(
      anno_collapsed,
      by = "TAIR"
    )
}


# The three contrasts that also have a complete featureCounts
# result table are exported and compared gene by gene.

compared_contrasts <- c(
  interaction = "Genotype x drought interaction",
  WT_drought = "WT drought vs watered",
  BRL3_drought = "BRL3 drought vs watered"
)

salmon_tables <- list()

for (contrast_name in names(compared_contrasts)) {
  salmon_tables[[contrast_name]] <- annotate_results(
    res_salmon[[contrast_name]]
  )

  write.table(
    salmon_tables[[contrast_name]],
    file = file.path(
      results_dir,
      paste0("salmon_", contrast_name, "_DESeq2_all_results.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 10. DEG counts for the five comparisons, both workflows
# ------------------------------------------------------------

is_deg <- function(padj, log2FoldChange) {
  !is.na(padj) &
    padj < padj_cutoff &
    abs(log2FoldChange) >= lfc_cutoff
}

count_degs <- function(res) {
  res_df <- as.data.frame(res)

  deg <- res_df[
    is_deg(res_df$padj, res_df$log2FoldChange),
  ]

  c(
    total_DEGs = nrow(deg),
    positive_or_up = sum(deg$log2FoldChange > 0),
    negative_or_down = sum(deg$log2FoldChange < 0)
  )
}

hisat2_summary <- read.delim(
  file.path(results_dir, "FULL_comparison_summary.tsv")
)

# Same row order as results/FULL_comparison_summary.tsv.
stopifnot(
  identical(
    hisat2_summary$comparison,
    c(
      "WT drought vs watered",
      "BRL3 drought vs watered",
      "BRL3 vs WT under watered",
      "BRL3 vs WT under drought",
      "Genotype x drought interaction"
    )
  )
)

salmon_deg_counts <- t(
  sapply(
    res_salmon[
      c(
        "WT_drought",
        "BRL3_drought",
        "BRL3_vs_WT_watered",
        "BRL3_vs_WT_drought",
        "interaction"
      )
    ],
    count_degs
  )
)

salmon_comparison_summary <- data.frame(
  comparison = hisat2_summary$comparison,
  HISAT2_total_DEGs = hisat2_summary$total_DEGs,
  HISAT2_positive_or_up = hisat2_summary$positive_or_up,
  HISAT2_negative_or_down = hisat2_summary$negative_or_down,
  Salmon_total_DEGs = salmon_deg_counts[, "total_DEGs"],
  Salmon_positive_or_up = salmon_deg_counts[, "positive_or_up"],
  Salmon_negative_or_down = salmon_deg_counts[, "negative_or_down"],
  row.names = NULL
)

salmon_comparison_summary

write.table(
  salmon_comparison_summary,
  file = file.path(results_dir, "salmon_comparison_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 11. Gene-by-gene comparison of the two workflows
# ============================================================

# Read the complete featureCounts-based result tables and join
# them with the Salmon tables by TAIR ID. A gene missing from
# one table was removed by the count prefilter in that workflow.

hisat2_files <- c(
  interaction = "FULL_interaction_DESeq2_all_results.tsv",
  WT_drought = "FULL_WT_drought_DESeq2_all_results.tsv",
  BRL3_drought = "FULL_BRL3_drought_DESeq2_all_results.tsv"
)

stats_columns <- c(
  "TAIR",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "padj"
)

# dplyr::count is written out in full: DESeq2 attaches
# matrixStats, which has its own count() and can mask it.
transcripts_per_gene <- tx2gene %>%
  dplyr::count(
    gene_id,
    name = "n_transcripts"
  )

comparison_tables <- list()

for (contrast_name in names(compared_contrasts)) {
  hisat2_table <- read.delim(
    file.path(results_dir, hisat2_files[[contrast_name]]),
    quote = ""
  )

  joined <- full_join(
    hisat2_table[, stats_columns],
    salmon_tables[[contrast_name]][, stats_columns],
    by = "TAIR",
    suffix = c("_HISAT2", "_Salmon")
  )

  joined$contrast <- contrast_name

  joined$DEG_HISAT2 <- is_deg(
    joined$padj_HISAT2,
    joined$log2FoldChange_HISAT2
  )

  joined$DEG_Salmon <- is_deg(
    joined$padj_Salmon,
    joined$log2FoldChange_Salmon
  )

  joined$DEG_status <- case_when(
    joined$DEG_HISAT2 & joined$DEG_Salmon ~ "Both",
    joined$DEG_HISAT2 ~ "HISAT2 only",
    joined$DEG_Salmon ~ "Salmon only",
    TRUE ~ "Neither"
  )

  comparison_tables[[contrast_name]] <- joined
}

comparison_all <- bind_rows(comparison_tables)


# ------------------------------------------------------------
# 12. Correlation of log2 fold changes
# ------------------------------------------------------------

# Two gene sets per contrast:
#   all shared genes: every gene with a log2FC in both workflows
#   DEG in either:    genes significant in at least one workflow
#
# Low-count genes have very noisy fold changes in both
# workflows, so the correlation over all shared genes is a
# conservative figure.

correlation_row <- function(df, contrast_name, gene_set) {
  data.frame(
    contrast = contrast_name,
    gene_set = gene_set,
    n_genes = nrow(df),
    pearson_r = cor(
      df$log2FoldChange_HISAT2,
      df$log2FoldChange_Salmon,
      method = "pearson"
    ),
    spearman_rho = cor(
      df$log2FoldChange_HISAT2,
      df$log2FoldChange_Salmon,
      method = "spearman"
    )
  )
}

log2fc_correlation <- bind_rows(
  lapply(
    names(compared_contrasts),
    function(contrast_name) {
      shared <- comparison_tables[[contrast_name]] %>%
        filter(
          !is.na(log2FoldChange_HISAT2),
          !is.na(log2FoldChange_Salmon)
        )

      bind_rows(
        correlation_row(
          shared,
          contrast_name,
          "all shared genes"
        ),
        correlation_row(
          shared %>% filter(DEG_status != "Neither"),
          contrast_name,
          "DEG in either workflow"
        )
      )
    }
  )
)

log2fc_correlation

write.table(
  log2fc_correlation,
  file = file.path(results_dir, "comparison_log2FC_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. log2 fold change scatter plot
# ------------------------------------------------------------

scatter_df <- comparison_all %>%
  filter(
    !is.na(log2FoldChange_HISAT2),
    !is.na(log2FoldChange_Salmon)
  )

scatter_df$DEG_status <- factor(
  scatter_df$DEG_status,
  levels = c(
    "Neither",
    "Both",
    "HISAT2 only",
    "Salmon only"
  )
)

# Draw "Neither" first so the DEG points stay visible.
scatter_df <- scatter_df %>%
  arrange(DEG_status)

panel_labels <- log2fc_correlation %>%
  filter(
    gene_set == "all shared genes"
  ) %>%
  mutate(
    panel = paste0(
      compared_contrasts[contrast],
      "\nPearson r = ", sprintf("%.3f", pearson_r),
      ", Spearman rho = ", sprintf("%.3f", spearman_rho),
      " (", format(n_genes, big.mark = ","), " genes)"
    )
  )

scatter_df$panel <- factor(
  panel_labels$panel[
    match(scatter_df$contrast, panel_labels$contrast)
  ],
  levels = panel_labels$panel
)

# A few genes have an extreme fold change in one workflow
# (|log2FC| > 10). This happens when a gene has zero counts in
# every sample of one group, so the fold change is not
# estimable and DESeq2 reports a very large value. They are
# kept in all tables and correlations, but the plot is limited
# to +/- 10 so that the other genes remain readable.

scatter_limit <- 10

n_outside_plot <- sum(
  abs(scatter_df$log2FoldChange_HISAT2) > scatter_limit |
    abs(scatter_df$log2FoldChange_Salmon) > scatter_limit
)

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = log2FoldChange_HISAT2,
    y = log2FoldChange_Salmon,
    color = DEG_status
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.5,
    size = 0.8
  ) +
  facet_wrap(
    ~ panel,
    ncol = 1
  ) +
  coord_cartesian(
    xlim = c(-scatter_limit, scatter_limit),
    ylim = c(-scatter_limit, scatter_limit)
  ) +
  scale_color_manual(
    values = c(
      "Neither" = "grey80",
      "Both" = "#EF4444",
      "HISAT2 only" = "#3B82F6",
      "Salmon only" = "#F59E0B"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2.5, alpha = 1)
    )
  ) +
  labs(
    title = "Gene-level log2 fold changes: HISAT2 + featureCounts vs Salmon",
    subtitle = paste0(
      "Same DESeq2 model in both workflows. Dashed line: identical estimates.\n",
      n_outside_plot, " points with |log2FC| > ", scatter_limit,
      " are outside the plotted range"
    ),
    x = "log2 fold change, HISAT2 + featureCounts",
    y = "log2 fold change, Salmon + tximport",
    color = paste0(
      "DEG (padj < ", padj_cutoff,
      ", |log2FC| >= ", lfc_cutoff, ")"
    )
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(hjust = 0)
  )

ggsave(
  file.path(results_dir, "comparison_log2FC_scatter.png"),
  plot = p_scatter,
  width = 8,
  height = 14,
  dpi = 300
)


# ------------------------------------------------------------
# 14. DEG overlap
# ------------------------------------------------------------

deg_overlap <- comparison_all %>%
  group_by(contrast) %>%
  summarise(
    HISAT2_DEGs = sum(DEG_HISAT2),
    Salmon_DEGs = sum(DEG_Salmon),
    shared_DEGs = sum(DEG_HISAT2 & DEG_Salmon),
    HISAT2_only = sum(DEG_HISAT2 & !DEG_Salmon),
    Salmon_only = sum(DEG_Salmon & !DEG_HISAT2),
    # Shared DEGs whose fold change has the same sign.
    shared_same_direction = sum(
      DEG_HISAT2 & DEG_Salmon &
        sign(log2FoldChange_HISAT2) == sign(log2FoldChange_Salmon)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    percent_of_HISAT2_DEGs_shared = round(100 * shared_DEGs / HISAT2_DEGs, 1),
    percent_of_Salmon_DEGs_shared = round(100 * shared_DEGs / Salmon_DEGs, 1),
    # Jaccard index: shared / union.
    jaccard_index = round(
      shared_DEGs / (HISAT2_DEGs + Salmon_DEGs - shared_DEGs),
      3
    )
  )

deg_overlap <- deg_overlap[
  match(names(compared_contrasts), deg_overlap$contrast),
]

as.data.frame(deg_overlap)

write.table(
  deg_overlap,
  file = file.path(results_dir, "comparison_DEG_overlap.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# UpSet plot: six DEG sets (three contrasts x two workflows).
# Each bar is the number of genes in exactly that combination
# of sets.

deg_sets <- list()

for (contrast_name in names(compared_contrasts)) {
  df <- comparison_tables[[contrast_name]]

  deg_sets[[paste0(contrast_name, " | HISAT2")]] <- df$TAIR[df$DEG_HISAT2]
  deg_sets[[paste0(contrast_name, " | Salmon")]] <- df$TAIR[df$DEG_Salmon]
}

png(
  file.path(results_dir, "comparison_DEG_overlap_upset.png"),
  width = 11,
  height = 6.5,
  units = "in",
  res = 300
)

print(
  upset(
    fromList(deg_sets),
    sets = rev(names(deg_sets)),
    keep.order = TRUE,
    nintersects = 25,
    order.by = "freq",
    mainbar.y.label = "Genes in this combination",
    sets.x.label = "DEGs per set",
    text.scale = 1.2
  )
)

dev.off()


# ------------------------------------------------------------
# 15. featureCounts assigned rate vs Salmon mapping rate
# ------------------------------------------------------------

# Both rates use the number of sequenced reads as denominator.
#
#   HISAT2 + featureCounts: reads assigned to a gene.
#     A read is not assigned when it is unmapped, maps to
#     several genomic positions (multi-mapping, discarded by
#     default), overlaps no annotated exon on the expected
#     strand, or overlaps more than one gene (ambiguous).
#
#   Salmon: reads assigned to an annotated transcript.
#     Multi-mapping reads are kept and divided between the
#     compatible transcripts. Reads that match the genome decoy
#     better than any transcript are not counted.

salmon_mapping <- read.delim(
  file.path(salmon_dir, "salmon_mapping_summary.tsv")
)

fc_summary <- read.delim(
  file.path(counts_dir, "all_11samples_featureCounts.txt.summary"),
  check.names = FALSE,
  row.names = 1
)

# Column names look like "alignment/SRR7779219.sorted.bam".
colnames(fc_summary) <- sub(
  "\\.sorted\\.bam$",
  "",
  basename(colnames(fc_summary))
)

stopifnot(
  all(salmon_mapping$sra_run %in% colnames(fc_summary))
)

mapping_rates <- salmon_mapping %>%
  transmute(
    sample_id,
    sra_run,
    total_reads = observed_fragments,
    featureCounts_assigned = unlist(
      fc_summary["Assigned", sra_run]
    ),
    featureCounts_assigned_percent = round(
      100 * featureCounts_assigned / total_reads,
      2
    ),
    Salmon_mapped = mapped_fragments,
    Salmon_mapping_rate_percent = round(
      100 * Salmon_mapped / total_reads,
      2
    ),
    difference_percentage_points = round(
      Salmon_mapping_rate_percent - featureCounts_assigned_percent,
      2
    ),
    Salmon_library_type = library_type
  )

mapping_rates

write.table(
  mapping_rates,
  file = file.path(results_dir, "comparison_mapping_rates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Genes significant in only one workflow
# ------------------------------------------------------------

# For every gene that is a DEG in exactly one workflow, look at
# what the other workflow reports for it and sort it into one
# of four descriptive categories:
#
#   1. not tested in the other workflow
#        The gene has no adjusted p-value there: removed by the
#        count prefilter or by DESeq2 independent filtering
#        (low counts), or flagged as a count outlier.
#
#   2. below the fold-change threshold in the other workflow
#        Significant there too (padj < 0.05), but |log2FC| < 1.
#
#   3. above the padj threshold in the other workflow
#        Same direction and |log2FC| >= 1 there, but padj >= 0.05.
#
#   4. different estimate
#        Everything else: the two workflows disagree more
#        clearly about this gene.
#
# Categories 2 and 3 are threshold effects: the two workflows
# agree, and the gene falls just on the other side of a cut-off.
#
# Two extra columns help to interpret the remaining cases:
#
#   log2_baseMean_ratio  log2(Salmon baseMean / HISAT2 baseMean).
#        Far from 0 means the two workflows count clearly
#        different numbers of reads for the gene. Typical
#        causes: reads shared with a paralog or another gene
#        family member (discarded by featureCounts, divided by
#        Salmon), or reads in a region where two genes overlap.
#
#   n_transcripts        annotated isoforms of the gene. With
#        several isoforms, Salmon's gene-level estimate depends
#        on how reads are divided between them, and tximport
#        corrects for changes in average transcript length.
#
# These are plausible explanations, not proof for any single gene.

pipeline_specific <- comparison_all %>%
  filter(
    DEG_status %in% c("HISAT2 only", "Salmon only")
  ) %>%
  mutate(
    other_padj = ifelse(
      DEG_status == "HISAT2 only",
      padj_Salmon,
      padj_HISAT2
    ),
    other_log2FC = ifelse(
      DEG_status == "HISAT2 only",
      log2FoldChange_Salmon,
      log2FoldChange_HISAT2
    ),
    own_log2FC = ifelse(
      DEG_status == "HISAT2 only",
      log2FoldChange_HISAT2,
      log2FoldChange_Salmon
    ),
    category = case_when(
      is.na(other_padj) ~
        "1. not tested in the other workflow",
      other_padj < padj_cutoff &
        sign(other_log2FC) == sign(own_log2FC) ~
        "2. below the fold-change threshold in the other workflow",
      sign(other_log2FC) == sign(own_log2FC) &
        abs(other_log2FC) >= lfc_cutoff ~
        "3. above the padj threshold in the other workflow",
      TRUE ~
        "4. different estimate"
    ),
    log2_baseMean_ratio = round(
      log2(baseMean_Salmon / baseMean_HISAT2),
      3
    ),
    large_count_difference = !is.na(log2_baseMean_ratio) &
      abs(log2_baseMean_ratio) >= 1
  ) %>%
  left_join(
    transcripts_per_gene,
    by = c("TAIR" = "gene_id")
  ) %>%
  left_join(
    anno_collapsed[, c("TAIR", "SYMBOL")],
    by = "TAIR"
  ) %>%
  dplyr::select(
    contrast,
    TAIR,
    SYMBOL,
    significant_in = DEG_status,
    category,
    baseMean_HISAT2,
    baseMean_Salmon,
    log2_baseMean_ratio,
    large_count_difference,
    n_transcripts,
    log2FoldChange_HISAT2,
    log2FoldChange_Salmon,
    padj_HISAT2,
    padj_Salmon
  ) %>%
  arrange(
    contrast,
    significant_in,
    category,
    TAIR
  )

write.table(
  pipeline_specific,
  file = file.path(results_dir, "comparison_pipeline_specific_DEGs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


pipeline_specific_summary <- pipeline_specific %>%
  group_by(
    contrast,
    significant_in,
    category
  ) %>%
  summarise(
    n_genes = n(),
    n_large_count_difference = sum(large_count_difference),
    n_multi_isoform = sum(n_transcripts > 1),
    .groups = "drop"
  )

as.data.frame(pipeline_specific_summary)

write.table(
  pipeline_specific_summary,
  file = file.path(results_dir, "comparison_pipeline_specific_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

cat("\nSalmon + tximport + DESeq2 analysis complete\n")

cat("\nGenes retained after prefilter (Salmon):", nrow(dds_salmon), "\n")

cat("\nDEG counts, both workflows:\n")
print(salmon_comparison_summary)

cat("\nlog2 fold change correlation:\n")
print(log2fc_correlation)

cat("\nDEG overlap:\n")
print(as.data.frame(deg_overlap))

cat("\nMean featureCounts assigned rate (%):",
    round(mean(mapping_rates$featureCounts_assigned_percent), 2), "\n")

cat("Mean Salmon mapping rate (%):",
    round(mean(mapping_rates$Salmon_mapping_rate_percent), 2), "\n")


# ============================================================
# End of Salmon / tximport / DESeq2 analysis
# ============================================================
