#!/usr/bin/env bash

# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Upstream workflow: QC, alignment, BAM processing, counting
#
# Dataset:
#   GSE119382 Arabidopsis root RNA-seq
#
# Reference genome:
#   TAIR10
#
# Gene annotation:
#   Ensembl Plants release 63 / Araport11
#
# Notes:
# - Raw FASTQ files are not stored in this repository.
# - Reference files and alignment files are excluded from Git.
# - This script documents the upstream commands used to generate
#   the gene-level count matrix used for downstream DESeq2 analysis.
# ============================================================

set -euo pipefail


# ------------------------------------------------------------
# 1. Create output directories
# ------------------------------------------------------------

mkdir -p qc/full_fastqc
mkdir -p qc/full_multiqc
mkdir -p alignment
mkdir -p counts


# ------------------------------------------------------------
# 2. FastQC
# ------------------------------------------------------------

# Run FastQC on all single-end FASTQ files.

fastqc \
  raw_data/*.fastq.gz \
  --outdir qc/full_fastqc


# ------------------------------------------------------------
# 3. MultiQC
# ------------------------------------------------------------

# Summarize FastQC reports.

multiqc \
  qc/full_fastqc \
  --outdir qc/full_multiqc \
  --filename full_multiqc_report.html


# ------------------------------------------------------------
# 4. HISAT2 alignment
# ------------------------------------------------------------

# Align each single-end sample to TAIR10 and pipe directly into
# samtools sort to avoid storing large intermediate SAM files.
#
# HISAT2 index prefix:
#   reference/hisat2_index/TAIR10

for fq in raw_data/*.fastq.gz
do
    sample=$(basename "$fq" .fastq.gz)

    hisat2 \
      -x reference/hisat2_index/TAIR10 \
      -U "$fq" \
      --summary-file "alignment/${sample}_hisat2_summary.txt" \
    | samtools sort \
      -o "alignment/${sample}.sorted.bam"

    samtools index \
      "alignment/${sample}.sorted.bam"

    samtools flagstat \
      "alignment/${sample}.sorted.bam" \
      > "alignment/${sample}.flagstat.txt"
done


# ------------------------------------------------------------
# 5. featureCounts
# ------------------------------------------------------------

# Empirical strandedness testing supported reverse-stranded
# counting, therefore -s 2 is used.
#
# Reads are summarized at gene level using exon features and
# gene_id attributes from the GTF.

featureCounts \
  -T 6 \
  -s 2 \
  -t exon \
  -g gene_id \
  -a reference/Arabidopsis_thaliana.TAIR10.63.gtf \
  -o counts/all_11samples_featureCounts.txt \
  alignment/SRR7779219.sorted.bam \
  alignment/SRR7779220.sorted.bam \
  alignment/SRR7779221.sorted.bam \
  alignment/SRR7779222.sorted.bam \
  alignment/SRR7779223.sorted.bam \
  alignment/SRR7779224.sorted.bam \
  alignment/SRR7779225.sorted.bam \
  alignment/SRR7779226.sorted.bam \
  alignment/SRR7779227.sorted.bam \
  alignment/SRR7779228.sorted.bam \
  alignment/SRR7779229.sorted.bam


# ------------------------------------------------------------
# 6. Extract gene-level count matrix
# ------------------------------------------------------------

# featureCounts output contains annotation columns before the
# sample count columns. Extract gene_id plus the 11 sample counts.

awk 'BEGIN{FS=OFS="\t"}
     !/^#/ {
       print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17
     }' \
  counts/all_11samples_featureCounts.txt \
  > counts/all_11samples_raw_counts_unlabeled.tsv


echo "Upstream RNA-seq workflow completed."
