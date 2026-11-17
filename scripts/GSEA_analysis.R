# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Gene set enrichment analysis (GSEA) of GO Biological Process
#
# Why GSEA in addition to ORA:
#   ORA (enrichGO) only sees genes that pass the DEG thresholds.
#   The interaction contrast has 137 DEGs, which is a small
#   input. GSEA is cut-off free: it ranks ALL tested genes and
#   asks whether the genes of a GO term sit together near the
#   top or bottom of the ranking. It can therefore detect
#   pathways whose genes shift modestly but consistently.
#
# Input (complete DESeq2 result tables; DESeq2 is not rerun):
#   results/FULL_interaction_DESeq2_all_results.tsv
#   results/FULL_WT_drought_DESeq2_all_results.tsv
#   results/FULL_BRL3_drought_DESeq2_all_results.tsv
#
#   All three come from the same factorial model
#   (~ genotype * condition). The last two are written by
#   scripts/FULL_factorial_export_contrasts.R.
#
# Output (all new files, prefix GSEA_):
#   results/GSEA_<contrast>_results.tsv
#   results/GSEA_<contrast>_dotplot.png
#   results/GSEA_<contrast>_enrichment_plot_positive.png
#   results/GSEA_<contrast>_enrichment_plot_negative.png
#   results/GSEA_summary.tsv
#   results/GSEA_vs_ORA_comparison.tsv
#   results/GSEA_vs_ORA_comparison.png
#
# Run from the project root (the folder with the .Rproj file).
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(org.At.tair.db)
library(clusterProfiler)
library(enrichplot)


# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

out_dir <- "results"

# Random seed. GSEA p-values are estimated by random sampling,
# so the seed is fixed before every gseGO() call.
gsea_seed <- 20261116

# Gene-set size limits, counted AFTER restricting each GO term
# to the genes present in the ranked list.
#   minGSSize = 15: very small sets give unstable scores.
#   maxGSSize = 500: very large sets are too general to interpret.
gsea_min_size <- 15
gsea_max_size <- 500

# Permutation scheme: "multilevel".
#   The ranking is fixed and gene labels are permuted (gene-set
#   permutation, the "preranked" form of GSEA). Sample labels
#   are NOT permuted: with 2-3 replicates per group that is not
#   possible. "multilevel" is the adaptive multilevel splitting
#   Monte Carlo scheme of fgsea, which estimates very small
#   p-values without a fixed number of permutations.
#
#   Consequence for interpretation: gene-label permutation
#   assumes genes are independent. Genes in one pathway are
#   usually co-expressed, so these p-values are optimistic and
#   the number of significant terms should not be over-read.
gsea_method <- "multilevel"

# Significance threshold for reporting.
gsea_padj_cutoff <- 0.05

# Number of terms per direction in the dot plot, and number of
# terms per enrichment plot.
n_dotplot_terms <- 10
n_enrichment_plot_terms <- 3


# ------------------------------------------------------------
# 3. Define the contrasts
# ------------------------------------------------------------

# For the two drought contrasts, NES > 0 means the gene set is
# shifted towards higher expression under drought.
#
# For the interaction, the sign has a different meaning:
#
#   (BRL3 drought - BRL3 watered) - (WT drought - WT watered)
#
#   NES > 0: the drought response of the gene set is MORE
#            POSITIVE in BRL3 than in WT.
#   NES < 0: the drought response is MORE NEGATIVE, or LESS
#            POSITIVE, in BRL3 than in WT.
#
# This is not the same as up- and downregulated. A set that is
# induced by drought in both genotypes, but less strongly in
# BRL3, has a negative interaction NES.

contrasts <- list(
  interaction = list(
    file = "results/FULL_interaction_DESeq2_all_results.tsv",
    title = "Genotype x drought interaction",
    positive_label = "NES > 0: drought response\nmore positive in BRL3 than WT",
    negative_label = "NES < 0: drought response\nmore negative in BRL3 than WT"
  ),
  WT_drought = list(
    file = "results/FULL_WT_drought_DESeq2_all_results.tsv",
    title = "WT drought vs watered",
    positive_label = "NES > 0: higher under drought",
    negative_label = "NES < 0: lower under drought"
  ),
  BRL3_drought = list(
    file = "results/FULL_BRL3_drought_DESeq2_all_results.tsv",
    title = "BRL3 drought vs watered",
    positive_label = "NES > 0: higher under drought",
    negative_label = "NES < 0: lower under drought"
  )
)

stopifnot(
  all(file.exists(sapply(contrasts, function(x) x$file)))
)


# ------------------------------------------------------------
# 4. Build the ranked gene list
# ------------------------------------------------------------

# Ranking metric: the DESeq2 Wald statistic (stat column),
#   stat = log2FoldChange / lfcSE.
# It carries the direction of the change and is scaled by its
# uncertainty, so noisy low-count genes with large fold changes
# do not dominate the ends of the ranking.
#
# Steps:
#   - drop genes with stat = NA
#   - keep one row per TAIR ID
#   - sort in decreasing order

make_ranked_list <- function(file) {
  res <- read.delim(
    file,
    quote = ""
  )

  stopifnot(
    all(c("TAIR", "stat") %in% colnames(res))
  )

  n_total <- nrow(res)

  res <- res %>%
    filter(
      !is.na(stat),
      !is.na(TAIR)
    )

  n_after_na <- nrow(res)

  # If a TAIR ID occurred more than once, keep the row with the
  # largest absolute statistic. (The DESeq2 tables have one row
  # per gene, so this normally removes nothing.)
  res <- res %>%
    arrange(
      dplyr::desc(abs(stat))
    ) %>%
    distinct(
      TAIR,
      .keep_all = TRUE
    )

  n_after_dedup <- nrow(res)

  # Sort by stat; TAIR ID as a second key makes the order of
  # tied genes deterministic.
  res <- res %>%
    arrange(
      dplyr::desc(stat),
      TAIR
    )

  ranked <- res$stat
  names(ranked) <- res$TAIR

  cat(
    "  genes in table:", n_total,
    "| after removing NA stat:", n_after_na,
    "| after de-duplication:", n_after_dedup,
    "| tied stat values:", sum(duplicated(ranked)),
    "\n"
  )

  ranked
}


# ------------------------------------------------------------
# 5. Run gseGO for one ranked list
# ------------------------------------------------------------

# pvalueCutoff = 1 keeps every tested GO term in the object.
# BH adjustment is done over all tested terms either way, so
# p.adjust values are unaffected. Keeping all terms lets the
# ORA comparison look up the GSEA result of a term even when it
# is not significant. Significant terms are filtered afterwards.

run_gsea <- function(ranked) {
  set.seed(gsea_seed)

  gseGO(
    geneList = ranked,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    minGSSize = gsea_min_size,
    maxGSSize = gsea_max_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    method = gsea_method,
    verbose = FALSE
  )
}


# ------------------------------------------------------------
# 6. Run all contrasts
# ------------------------------------------------------------

# About 30 seconds per contrast.
#
# Two warnings are expected and harmless:
#   - "There are ties in the preranked stats": a handful of
#     genes share the same stat value.
#   - "P-values are less than 1e-10": the most significant
#     terms hit the lower limit of the p-value estimate, so
#     their p-values are reported as upper bounds.

ranked_lists <- list()
gsea_objects <- list()
gsea_all_terms <- list()
gsea_significant <- list()

for (contrast_name in names(contrasts)) {
  cat("\n", contrast_name, "\n", sep = "")

  ranked_lists[[contrast_name]] <- make_ranked_list(
    contrasts[[contrast_name]]$file
  )

  gsea_objects[[contrast_name]] <- run_gsea(
    ranked_lists[[contrast_name]]
  )

  gsea_all_terms[[contrast_name]] <- as.data.frame(
    gsea_objects[[contrast_name]]
  ) %>%
    arrange(
      pvalue,
      dplyr::desc(abs(NES))
    )

  gsea_significant[[contrast_name]] <- gsea_all_terms[[contrast_name]] %>%
    filter(
      p.adjust < gsea_padj_cutoff
    )

  cat(
    "  GO terms tested:", nrow(gsea_all_terms[[contrast_name]]),
    "| significant:", nrow(gsea_significant[[contrast_name]]),
    "| NES > 0:", sum(gsea_significant[[contrast_name]]$NES > 0),
    "| NES < 0:", sum(gsea_significant[[contrast_name]]$NES < 0),
    "\n"
  )
}


# ------------------------------------------------------------
# 7. Export result tables
# ------------------------------------------------------------

# One table per contrast with the significant terms.
# core_enrichment lists the leading-edge genes (TAIR IDs): the
# genes that drive the enrichment score of the term.

for (contrast_name in names(contrasts)) {
  write.table(
    gsea_significant[[contrast_name]],
    file = file.path(
      out_dir,
      paste0("GSEA_", contrast_name, "_results.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


gsea_summary <- data.frame(
  contrast = names(contrasts),
  genes_ranked = sapply(ranked_lists, length),
  GO_terms_tested = sapply(gsea_all_terms, nrow),
  significant_terms = sapply(gsea_significant, nrow),
  NES_positive = sapply(gsea_significant, function(x) sum(x$NES > 0)),
  NES_negative = sapply(gsea_significant, function(x) sum(x$NES < 0)),
  row.names = NULL
)

gsea_summary

write.table(
  gsea_summary,
  file = file.path(out_dir, "GSEA_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Dot plots, split by the sign of NES
# ------------------------------------------------------------

# The top terms of each direction are shown in separate panels.
# x = NES, point size = number of genes of the term in the
# ranked list, colour = adjusted p-value.

plot_gsea_dotplot <- function(significant, contrast) {
  plot_df <- significant %>%
    mutate(
      direction = ifelse(
        NES > 0,
        contrast$positive_label,
        contrast$negative_label
      )
    ) %>%
    group_by(direction) %>%
    slice_min(
      order_by = pvalue,
      n = n_dotplot_terms,
      with_ties = FALSE
    ) %>%
    ungroup()

  plot_df$direction <- factor(
    plot_df$direction,
    levels = c(
      contrast$negative_label,
      contrast$positive_label
    )
  )

  plot_df$Description <- factor(
    plot_df$Description,
    levels = unique(
      plot_df$Description[order(plot_df$NES)]
    )
  )

  ggplot(
    plot_df,
    aes(
      x = NES,
      y = Description,
      size = setSize,
      color = p.adjust
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ direction,
      scales = "free",
      ncol = 1
    ) +
    scale_color_gradient(
      low = "#EF4444",
      high = "#3B82F6",
      trans = "log10"
    ) +
    scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = 45)
    ) +
    labs(
      title = paste0("GSEA of GO Biological Process: ", contrast$title),
      subtitle = paste0(
        "Top ", n_dotplot_terms,
        " terms per direction, adjusted p-value < ", gsea_padj_cutoff
      ),
      x = "Normalized enrichment score (NES)",
      y = NULL,
      size = "Gene set size",
      color = "Adjusted p-value"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(hjust = 0)
    )
}


for (contrast_name in names(contrasts)) {
  if (nrow(gsea_significant[[contrast_name]]) == 0) {
    cat("No significant terms for", contrast_name, "- dot plot skipped.\n")
    next
  }

  ggsave(
    file.path(
      out_dir,
      paste0("GSEA_", contrast_name, "_dotplot.png")
    ),
    plot = plot_gsea_dotplot(
      gsea_significant[[contrast_name]],
      contrasts[[contrast_name]]
    ),
    width = 9,
    height = 9,
    dpi = 300
  )
}


# ------------------------------------------------------------
# 9. Enrichment plots (running score + barcode)
# ------------------------------------------------------------

# For the most significant terms of each direction:
#   top panel:    running enrichment score along the ranking
#   middle panel: position of the term's genes (barcode)
#   bottom panel: the ranking metric (Wald statistic)
#
# A positive-NES term peaks on the left (genes concentrated at
# the top of the ranking); a negative-NES term dips on the right.

for (contrast_name in names(contrasts)) {
  significant <- gsea_significant[[contrast_name]]

  for (direction in c("positive", "negative")) {
    if (direction == "positive") {
      top_terms <- significant %>% filter(NES > 0)
    } else {
      top_terms <- significant %>% filter(NES < 0)
    }

    top_terms <- head(
      top_terms,
      n_enrichment_plot_terms
    )

    if (nrow(top_terms) == 0) {
      cat(
        "No significant", direction, "NES terms for",
        contrast_name, "- enrichment plot skipped.\n"
      )
      next
    }

    p_enrichment <- gseaplot2(
      gsea_objects[[contrast_name]],
      geneSetID = top_terms$ID,
      title = paste0(
        contrasts[[contrast_name]]$title,
        ": top ", direction, " NES terms"
      ),
      pvalue_table = FALSE
    )

    ggsave(
      file.path(
        out_dir,
        paste0(
          "GSEA_", contrast_name,
          "_enrichment_plot_", direction, ".png"
        )
      ),
      plot = p_enrichment,
      width = 9,
      height = 7,
      dpi = 300
    )
  }
}


# ============================================================
# 10. Compare GSEA with the existing ORA results
# ============================================================

# ORA tables already in results/ (not modified here):
#
#   Interaction, factorial model:
#     FULL_interaction_GO_enrichment_positive.tsv
#     FULL_interaction_GO_enrichment_negative.tsv
#
#   WT drought vs watered, standalone WT-only model:
#     WT_GO_enrichment_upregulated.tsv
#     WT_GO_enrichment_downregulated.tsv
#
# Note on the WT comparison: the ORA used DEGs from the WT-only
# model (1,480 DEGs), while this GSEA ranks genes from the
# factorial model. The biological contrast is the same, but the
# two inputs are not from the same fit.
#
# There is no ORA for BRL3 drought vs watered, so that contrast
# is not part of this comparison.

ora_files <- list(
  interaction = c(
    positive = "results/FULL_interaction_GO_enrichment_positive.tsv",
    negative = "results/FULL_interaction_GO_enrichment_negative.tsv"
  ),
  WT_drought = c(
    positive = "results/WT_GO_enrichment_upregulated.tsv",
    negative = "results/WT_GO_enrichment_downregulated.tsv"
  )
)


read_ora <- function(files) {
  ora_list <- lapply(
    names(files),
    function(direction) {
      ora <- read.delim(
        files[[direction]],
        quote = ""
      )

      if (nrow(ora) == 0) {
        return(NULL)
      }

      data.frame(
        ID = ora$ID,
        Description = ora$Description,
        ORA_direction = direction,
        ORA_p.adjust = ora$p.adjust,
        ORA_gene_count = ora$Count
      )
    }
  )

  bind_rows(ora_list)
}


compare_gsea_ora <- function(contrast_name) {
  ora <- read_ora(
    ora_files[[contrast_name]]
  )

  gsea_all <- gsea_all_terms[[contrast_name]] %>%
    transmute(
      ID,
      GSEA_Description = Description,
      GSEA_setSize = setSize,
      GSEA_NES = NES,
      GSEA_p.adjust = p.adjust
    )

  gsea_sig_ids <- gsea_significant[[contrast_name]]$ID

  # Every term that is significant in at least one method.
  # A term can appear twice in ORA (once per direction); those
  # rows are kept separately.
  comparison <- full_join(
    ora,
    gsea_all %>%
      filter(
        ID %in% union(gsea_sig_ids, ora$ID)
      ),
    by = "ID"
  )

  comparison %>%
    mutate(
      contrast = contrast_name,
      Description = ifelse(
        is.na(Description),
        GSEA_Description,
        Description
      ),
      significant_in_ORA = !is.na(ORA_p.adjust),
      # NA here means the term was outside the GSEA size limits
      # and was therefore not tested.
      tested_in_GSEA = !is.na(GSEA_NES),
      significant_in_GSEA = ID %in% gsea_sig_ids,
      GSEA_direction = case_when(
        is.na(GSEA_NES) ~ NA_character_,
        GSEA_NES > 0 ~ "positive",
        TRUE ~ "negative"
      ),
      same_direction = ifelse(
        significant_in_ORA & tested_in_GSEA,
        ORA_direction == GSEA_direction,
        NA
      ),
      category = case_when(
        significant_in_ORA & significant_in_GSEA & same_direction ~
          "both",
        significant_in_ORA & significant_in_GSEA ~
          "both, opposite direction",
        significant_in_ORA & !tested_in_GSEA ~
          "ORA only (not tested in GSEA: set size outside limits)",
        significant_in_ORA ~
          "ORA only",
        TRUE ~
          "GSEA only"
      )
    ) %>%
    dplyr::select(
      contrast,
      ID,
      Description,
      category,
      ORA_direction,
      ORA_p.adjust,
      ORA_gene_count,
      GSEA_direction,
      GSEA_NES,
      GSEA_p.adjust,
      GSEA_setSize
    ) %>%
    arrange(
      category,
      GSEA_p.adjust
    )
}


gsea_vs_ora <- bind_rows(
  lapply(
    names(ora_files),
    compare_gsea_ora
  )
)

write.table(
  gsea_vs_ora,
  file = file.path(out_dir, "GSEA_vs_ORA_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Counts per category
gsea_vs_ora_counts <- gsea_vs_ora %>%
  dplyr::count(
    contrast,
    category
  )

gsea_vs_ora_counts


# The ORA terms and what GSEA says about each of them
gsea_vs_ora %>%
  filter(
    !is.na(ORA_p.adjust),
    contrast == "interaction"
  )


# ------------------------------------------------------------
# 11. GSEA vs ORA comparison figure
# ------------------------------------------------------------

# One point per GO term tested by GSEA.
#   x = NES, y = -log10 adjusted GSEA p-value.
#   Points above the dashed line are significant in GSEA.
#   Colour shows what ORA found for the same term.
#
# Reading the figure:
#   red above the line   = found by both methods
#   orange above the line = found by GSEA only
#   blue below the line  = found by ORA only
#
# ORA terms that GSEA did not test (gene set outside the size
# limits) have no NES and cannot be drawn; their number is
# given in the panel label.

comparison_plot_df <- bind_rows(
  lapply(
    names(ora_files),
    function(contrast_name) {
      ora_ids <- unique(
        read_ora(ora_files[[contrast_name]])$ID
      )

      n_ora_not_tested <- sum(
        !ora_ids %in% gsea_all_terms[[contrast_name]]$ID
      )

      gsea_all_terms[[contrast_name]] %>%
        mutate(
          panel = paste0(
            contrasts[[contrast_name]]$title,
            "\n(", n_ora_not_tested,
            " ORA terms not tested by GSEA are not shown)"
          ),
          in_ORA = ID %in% ora_ids,
          in_GSEA = p.adjust < gsea_padj_cutoff,
          category = case_when(
            in_ORA & in_GSEA ~ "Both",
            in_GSEA ~ "GSEA only",
            in_ORA ~ "ORA only",
            TRUE ~ "Neither"
          )
        )
    }
  )
)

comparison_plot_df$category <- factor(
  comparison_plot_df$category,
  levels = c(
    "Neither",
    "GSEA only",
    "ORA only",
    "Both"
  )
)

# Draw the highlighted categories on top of the others.
comparison_plot_df <- comparison_plot_df %>%
  arrange(category)


# Label the ORA terms: up to 8 per panel and category, chosen
# by the strongest GSEA result.
comparison_label_df <- comparison_plot_df %>%
  filter(
    category %in% c("Both", "ORA only")
  ) %>%
  group_by(
    panel,
    category
  ) %>%
  slice_min(
    order_by = pvalue,
    n = 8,
    with_ties = FALSE
  ) %>%
  ungroup()


p_gsea_vs_ora <- ggplot(
  comparison_plot_df,
  aes(
    x = NES,
    y = -log10(p.adjust),
    color = category
  )
) +
  geom_point(
    aes(size = category),
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = -log10(gsea_padj_cutoff),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0
  ) +
  ggrepel::geom_text_repel(
    data = comparison_label_df,
    aes(label = Description),
    size = 2.8,
    bg.color = "white",
    bg.r = 0.15,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ panel,
    ncol = 1
  ) +
  scale_color_manual(
    values = c(
      "Neither" = "grey85",
      "GSEA only" = "#F59E0B",
      "ORA only" = "#3B82F6",
      "Both" = "#EF4444"
    )
  ) +
  scale_size_manual(
    values = c(
      "Neither" = 0.8,
      "GSEA only" = 1.2,
      "ORA only" = 2.2,
      "Both" = 2.2
    )
  ) +
  labs(
    title = "GO Biological Process terms: GSEA compared with ORA",
    subtitle = paste0(
      "Each point is one GO term. Dashed line: GSEA adjusted p-value = ",
      gsea_padj_cutoff
    ),
    x = "GSEA normalized enrichment score (NES)",
    y = "-log10 GSEA adjusted p-value",
    color = "Significant in",
    size = "Significant in"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(hjust = 0)
  )


ggsave(
  file.path(out_dir, "GSEA_vs_ORA_comparison.png"),
  plot = p_gsea_vs_ora,
  width = 9,
  height = 11,
  dpi = 300
)


# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\nGSEA complete\n\n")

print(gsea_summary)

cat("\nGSEA vs ORA:\n")

print(as.data.frame(gsea_vs_ora_counts))


# ============================================================
# End of GSEA analysis
# ============================================================
