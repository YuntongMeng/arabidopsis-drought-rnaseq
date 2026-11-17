#!/usr/bin/env bash

# ============================================================
# Arabidopsis drought RNA-seq reanalysis
# Alternative upstream workflow: Salmon transcript quantification
#
# Purpose:
#   Quantify the same 11 single-end samples with a
#   transcriptome-based method (Salmon selective alignment) so
#   that the results can be compared with the genome-based
#   workflow (HISAT2 + featureCounts) in
#   scripts/upstream_rnaseq_workflow.sh.
#
# Dataset:
#   GSE119382 Arabidopsis root RNA-seq
#
# Reference (same annotation version as the featureCounts GTF):
#   Ensembl Plants release 63, TAIR10
#     cDNA FASTA   protein-coding transcripts
#     ncRNA FASTA  non-coding transcripts
#     genome FASTA used only as decoy sequence
#     GTF          used only to check IDs and build tx2gene
#
#   cDNA and ncRNA are merged because featureCounts counted
#   every gene in the GTF, coding and non-coding. With both
#   files the Salmon reference covers the same gene set.
#
# Notes:
# - Raw FASTQ files, reference files and the Salmon index are
#   not stored in this repository.
# - Run from the project root with the conda environment that
#   provides salmon activated:
#
#     conda activate rnaseq
#     bash scripts/upstream_salmon_workflow.sh
#
# - Paths and thread count can be overridden with environment
#   variables, for example:
#
#     THREADS=4 bash scripts/upstream_salmon_workflow.sh
# ============================================================

set -euo pipefail


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

THREADS="${THREADS:-8}"

REF_DIR="${REF_DIR:-reference}"
FASTQ_DIR="${FASTQ_DIR:-raw_data}"
INDEX_DIR="${INDEX_DIR:-${REF_DIR}/salmon_index}"
QUANT_DIR="${QUANT_DIR:-salmon_quant}"
METADATA="${METADATA:-scripts/sample_metadata.tsv}"

CDNA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
NCRNA="${REF_DIR}/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz"
GENOME="${REF_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GTF="${REF_DIR}/Arabidopsis_thaliana.TAIR10.63.gtf"

# The reference FASTA files were downloaded from:
#   https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/cdna/
#   https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/ncrna/

# Salmon k-mer length for the index.
KMER="${KMER:-31}"

# featureCounts was run with -s 2 (reverse-stranded). For
# single-end reads the matching Salmon library type is SR.
EXPECTED_LIBTYPE="SR"


# ------------------------------------------------------------
# 2. Check inputs and create output directories
# ------------------------------------------------------------

for required_file in "$CDNA" "$NCRNA" "$GENOME" "$GTF" "$METADATA"
do
    if [ ! -s "$required_file" ]
    then
        echo "ERROR: required input not found: $required_file" >&2
        exit 1
    fi
done

mkdir -p "$QUANT_DIR"


# ------------------------------------------------------------
# 3. Record the Salmon version
# ------------------------------------------------------------

salmon --version | tee "${QUANT_DIR}/salmon_version.txt"


# ------------------------------------------------------------
# 4. Check that the FASTA and GTF describe the same transcripts
# ------------------------------------------------------------

# A cDNA/ncRNA FASTA from a different Ensembl release than the
# GTF would give transcript or gene IDs that do not match the
# featureCounts gene IDs. The check compares the complete
# transcript-to-gene mapping from both sources and stops if
# they differ in either direction.
#
# FASTA header example:
#   >AT1G01010.1 cdna chromosome:TAIR10:1:3631:5899:1 gene:AT1G01010 ...
#
# The transcript-to-gene table taken from the GTF is saved as
# tx2gene for tximport.

TX2GENE="${QUANT_DIR}/tx2gene_TAIR10.63.tsv"
FASTA_TX2GENE="${QUANT_DIR}/.fasta_tx2gene.tmp"

gzip -dc "$CDNA" "$NCRNA" \
| awk '/^>/ {
         transcript = substr($1, 2)
         gene = ""
         for (i = 2; i <= NF; i++) {
           if ($i ~ /^gene:/) gene = substr($i, 6)
         }
         print transcript "\t" gene
       }' \
| LC_ALL=C sort -u \
> "$FASTA_TX2GENE"

awk -F'\t' '$3 == "transcript" {
       match($9, /gene_id "[^"]+"/)
       gene = substr($9, RSTART + 9, RLENGTH - 10)
       match($9, /transcript_id "[^"]+"/)
       transcript = substr($9, RSTART + 15, RLENGTH - 16)
       print transcript "\t" gene
     }' "$GTF" \
| LC_ALL=C sort -u \
> "$TX2GENE"

n_fasta=$(wc -l < "$FASTA_TX2GENE" | tr -d ' ')
n_gtf=$(wc -l < "$TX2GENE" | tr -d ' ')
n_only_fasta=$(LC_ALL=C comm -23 "$FASTA_TX2GENE" "$TX2GENE" | wc -l | tr -d ' ')
n_only_gtf=$(LC_ALL=C comm -13 "$FASTA_TX2GENE" "$TX2GENE" | wc -l | tr -d ' ')
n_genes=$(cut -f2 "$TX2GENE" | LC_ALL=C sort -u | wc -l | tr -d ' ')

echo "Transcripts in cDNA + ncRNA FASTA: ${n_fasta}"
echo "Transcripts in GTF:                ${n_gtf}"
echo "Genes in GTF:                      ${n_genes}"
echo "Transcript-gene pairs only in FASTA: ${n_only_fasta}"
echo "Transcript-gene pairs only in GTF:   ${n_only_gtf}"

if [ "$n_fasta" -eq 0 ] || [ "$n_only_fasta" -ne 0 ] || [ "$n_only_gtf" -ne 0 ]
then
    echo "ERROR: the FASTA files and the GTF do not describe the same" >&2
    echo "       transcripts. Check that both come from Ensembl Plants" >&2
    echo "       release 63." >&2
    exit 1
fi

rm -f "$FASTA_TX2GENE"

echo "Reference check passed: FASTA and GTF are consistent."


# ------------------------------------------------------------
# 5. Build the decoy-aware Salmon index
# ------------------------------------------------------------

# The index contains the transcripts (cDNA + ncRNA) followed by
# the genome as decoy. A read that matches an unannotated
# genomic region better than any transcript is then assigned to
# the decoy and discarded, instead of being forced onto a
# similar-looking transcript. This is the closest equivalent to
# HISAT2 aligning against the whole genome.
#
# Salmon requires the decoy sequences to come last in the FASTA
# and their names to be listed in a separate file.
#
# --keepDuplicates: by default Salmon keeps only one copy of
# transcripts with identical sequences. About 550 TAIR10
# transcripts are exact duplicates of another one, some of them
# from different genes. Dropping them would remove transcripts
# (and possibly whole genes) from the output and give all reads
# to the retained copy. Keeping them leaves every transcript of
# the GTF in quant.sf, and reads from identical sequences are
# shared between them.

if [ -d "$INDEX_DIR" ] && [ -n "$(ls -A "$INDEX_DIR" 2>/dev/null)" ]
then
    echo "Salmon index already exists, skipping build: $INDEX_DIR"
else
    DECOYS="${REF_DIR}/salmon_decoys.txt"
    GENTROME="${REF_DIR}/salmon_gentrome.fa.gz"

    gzip -dc "$GENOME" \
    | awk '/^>/ { print substr($1, 2) }' \
    > "$DECOYS"

    # Concatenated gzip files form a valid gzip file.
    cat "$CDNA" "$NCRNA" "$GENOME" > "$GENTROME"

    salmon index \
      --transcripts "$GENTROME" \
      --decoys "$DECOYS" \
      --index "$INDEX_DIR" \
      --kmerLen "$KMER" \
      --keepDuplicates \
      --threads "$THREADS"

    rm -f "$GENTROME"
fi


# ------------------------------------------------------------
# 6. Quantify each sample
# ------------------------------------------------------------

# -l A lets Salmon detect the library type from the first
# mappings. The detected type is checked in the next step.
#
# One output directory per SRA run:
#   salmon_quant/<SRR>/quant.sf
#
# The metadata file is read on file descriptor 3 so that salmon
# cannot consume it through standard input.

while IFS=$'\t' read -r sample_id sra_run genotype condition replicate <&3
do
    if [ "$sample_id" = "sample_id" ]
    then
        continue
    fi

    fastq="${FASTQ_DIR}/${sra_run}.fastq.gz"

    if [ ! -s "$fastq" ]
    then
        echo "ERROR: FASTQ not found for ${sample_id}: $fastq" >&2
        exit 1
    fi

    echo "Quantifying ${sample_id} (${sra_run})"

    salmon quant \
      --index "$INDEX_DIR" \
      --libType A \
      --unmatedReads "$fastq" \
      --output "${QUANT_DIR}/${sra_run}" \
      --threads "$THREADS"
done 3< "$METADATA"


# ------------------------------------------------------------
# 7. Check library type and summarize mapping rates
# ------------------------------------------------------------

# Salmon writes, for every sample:
#   lib_format_counts.json   detected library type
#   logs/salmon_quant.log    observed and mapped fragments
#
# The detected library type must be SR (single-end, reverse
# strand), which corresponds to featureCounts -s 2. The script
# stops if any sample was detected differently.
#
# The mapping rate is the percentage of reads assigned to an
# annotated transcript. Reads that match the genome decoy
# better than any transcript are not counted as mapped.

SUMMARY="${QUANT_DIR}/salmon_mapping_summary.tsv"

printf 'sample_id\tsra_run\tlibrary_type\tobserved_fragments\tmapped_fragments\tmapping_rate_percent\n' \
> "$SUMMARY"

libtype_ok=1

while IFS=$'\t' read -r sample_id sra_run genotype condition replicate <&3
do
    if [ "$sample_id" = "sample_id" ]
    then
        continue
    fi

    sample_dir="${QUANT_DIR}/${sra_run}"

    library_type=$(
      sed -n 's/.*"expected_format": *"\([A-Z]*\)".*/\1/p' \
        "${sample_dir}/lib_format_counts.json"
    )

    observed=$(
      sed -n 's/^observed fragments: *//p' \
        "${sample_dir}/logs/salmon_quant.log"
    )

    mapped=$(
      sed -n 's/^mapped fragments: *//p' \
        "${sample_dir}/logs/salmon_quant.log"
    )

    mapping_rate=$(
      sed -n 's/^mapping rate: *\([0-9.]*\)%.*/\1/p' \
        "${sample_dir}/logs/salmon_quant.log"
    )

    if [ -z "$library_type" ] || [ -z "$observed" ] || [ -z "$mapped" ] || [ -z "$mapping_rate" ]
    then
        echo "ERROR: could not read Salmon logs for ${sample_id} (${sra_run})" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample_id" "$sra_run" "$library_type" "$observed" "$mapped" "$mapping_rate" \
    >> "$SUMMARY"

    if [ "$library_type" != "$EXPECTED_LIBTYPE" ]
    then
        echo "ERROR: ${sample_id} detected as ${library_type}, expected ${EXPECTED_LIBTYPE}" >&2
        libtype_ok=0
    fi
done 3< "$METADATA"

cat "$SUMMARY"

if [ "$libtype_ok" -ne 1 ]
then
    echo "ERROR: library type check failed." >&2
    exit 1
fi

echo "Library type check passed: all samples detected as ${EXPECTED_LIBTYPE}."


echo "Salmon upstream workflow completed."
