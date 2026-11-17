# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Full factorial model: export complete result tables for the
# two drought-vs-watered contrasts
#
# Why this script exists:
#   scripts/FULL_factorial_DESeq2_analysis.R saves the complete
#   result table only for the interaction term. The complete
#   tables for "WT drought vs watered" and "BRL3 drought vs
#   watered" from the SAME factorial model are needed as input
#   for GSEA and for the HISAT2-vs-Salmon comparison.
#
#   results/WT_DESeq2_all_results.tsv is not a substitute: it
#   comes from the standalone WT-only model (~ condition).
#
# This script refits exactly the same model as the main
# factorial script, checks that it reproduces the published
# numbers, and writes two NEW files. It does not modify any
# existing result file.
#
# Input:
#   counts/all_11samples_raw_counts.tsv
#   scripts/sample_metadata.tsv
#
# Output:
#   results/FULL_WT_drought_DESeq2_all_results.tsv
#   results/FULL_BRL3_drought_DESeq2_all_results.tsv
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(org.At.tair.db)
library(AnnotationDbi)


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

metadata_11 <- metadata_all[
  match(colnames(counts_11), metadata_all$sample_id),
]

rownames(metadata_11) <- metadata_11$sample_id

stopifnot(
  all(colnames(counts_11) == rownames(metadata_11))
)


# ------------------------------------------------------------
# 3. Fit the same factorial model as the main script
# ------------------------------------------------------------

# Reference levels, design, and prefilter are identical to
# scripts/FULL_factorial_DESeq2_analysis.R.

metadata_11$genotype <- factor(
  metadata_11$genotype,
  levels = c("WT", "BRL3")
)

metadata_11$condition <- factor(
  metadata_11$condition,
  levels = c("watered", "drought")
)

dds_11 <- DESeqDataSetFromMatrix(
  countData = counts_11,
  colData = metadata_11,
  design = ~ genotype * condition
)

keep_11 <- rowSums(counts(dds_11)) >= 10

dds_11 <- dds_11[keep_11, ]

dds_11 <- DESeq(dds_11)

resultsNames(dds_11)


# ------------------------------------------------------------
# 4. Extract the three contrasts
# ------------------------------------------------------------

# Interaction: only used here to verify the refit.
res_interaction <- results(
  dds_11,
  name = "genotypeBRL3.conditiondrought"
)

# WT drought vs watered = condition main effect
# (WT is the reference genotype).
res_WT_drought <- results(
  dds_11,
  name = "condition_drought_vs_watered"
)

# BRL3 drought vs watered =
# condition main effect + genotype:drought interaction.
res_BRL3_drought <- results(
  dds_11,
  contrast = list(
    c(
      "condition_drought_vs_watered",
      "genotypeBRL3.conditiondrought"
    )
  )
)


# ------------------------------------------------------------
# 5. Verify that the refit reproduces the existing results
# ------------------------------------------------------------

# 5a. Interaction statistics must match the saved table.

existing_interaction <- read.delim(
  "results/FULL_interaction_DESeq2_all_results.tsv",
  quote = ""
)

stopifnot(
  nrow(existing_interaction) == nrow(res_interaction),
  all(existing_interaction$TAIR == rownames(res_interaction))
)

for (column in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
  stopifnot(
    isTRUE(
      all.equal(
        existing_interaction[[column]],
        res_interaction[[column]],
        tolerance = 1e-8
      )
    )
  )
}


# 5b. DEG counts must match results/FULL_comparison_summary.tsv.

count_degs <- function(res) {
  res_df <- as.data.frame(res)

  deg <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) >= 1
    )

  c(
    total_DEGs = nrow(deg),
    positive_or_up = sum(deg$log2FoldChange > 0),
    negative_or_down = sum(deg$log2FoldChange < 0)
  )
}

existing_summary <- read.delim(
  "results/FULL_comparison_summary.tsv"
)

check_against_summary <- function(res, comparison_name) {
  expected <- existing_summary[
    existing_summary$comparison == comparison_name,
    c("total_DEGs", "positive_or_up", "negative_or_down")
  ]

  stopifnot(
    nrow(expected) == 1,
    all(count_degs(res) == unlist(expected))
  )
}

check_against_summary(res_WT_drought, "WT drought vs watered")
check_against_summary(res_BRL3_drought, "BRL3 drought vs watered")
check_against_summary(res_interaction, "Genotype x drought interaction")

cat("\nRefit reproduces the existing interaction table and DEG counts.\n")


# ------------------------------------------------------------
# 6. Add Arabidopsis gene annotation
# ------------------------------------------------------------

# Same annotation and one-row-per-gene collapsing as the main
# factorial script, so the new tables have identical columns to
# results/FULL_interaction_DESeq2_all_results.tsv.

annotate_results <- function(res) {
  res_df <- as.data.frame(res)

  res_df$TAIR <- rownames(res_df)

  anno <- AnnotationDbi::select(
    org.At.tair.db,
    keys = res_df$TAIR,
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

  res_df %>%
    left_join(
      anno_collapsed,
      by = "TAIR"
    )
}

res_WT_drought_annotated <- annotate_results(res_WT_drought)
res_BRL3_drought_annotated <- annotate_results(res_BRL3_drought)


# ------------------------------------------------------------
# 7. Export the two new complete result tables
# ------------------------------------------------------------

write.table(
  res_WT_drought_annotated,
  file = "results/FULL_WT_drought_DESeq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_BRL3_drought_annotated,
  file = "results/FULL_BRL3_drought_DESeq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Final summary
# ------------------------------------------------------------

cat("\nGenes in the fitted model:", nrow(dds_11), "\n")

cat("\nWT drought vs watered:\n")
print(count_degs(res_WT_drought))

cat("\nBRL3 drought vs watered:\n")
print(count_degs(res_BRL3_drought))

cat("\nGenotype x drought interaction:\n")
print(count_degs(res_interaction))


# ============================================================
# End of factorial contrast export
# ============================================================
