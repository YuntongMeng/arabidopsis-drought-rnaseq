# Arabidopsis Drought RNA-seq Reanalysis

An independent reanalysis of publicly available *Arabidopsis thaliana* RNA-seq data to investigate transcriptional responses to drought and the effects of BRL3 overexpression.

The project aims to build a reproducible workflow from raw sequencing reads to differential expression analysis and biological interpretation.

## Dataset

RNA-seq data were obtained from **NCBI Gene Expression Omnibus (GEO), accession GSE119382**.

The original experiment used a factorial design comparing:

- wild type (Col-0) and 35S:BRL3 plants
- watered and drought conditions
- biological replicates from *Arabidopsis thaliana* roots

The public GEO record currently contains 11 RNA-seq samples.

## Analysis

The analysis is organized in two stages.

### Stage 1 — Wild-type drought response

Six wild-type samples are used to establish the RNA-seq workflow:

- 3 drought biological replicates
- 3 watered biological replicates

### Stage 2 — BRL3 overexpression and drought response

The analysis will be extended to the full public dataset to investigate:

- transcriptional responses to drought
- expression changes associated with BRL3 overexpression
- genotype × condition interaction effects
- genes whose drought response differs between wild-type and BRL3-overexpressing plants

The workflow is:

```text
FASTQ
  ↓
FastQC / MultiQC
  ↓
HISAT2 alignment to TAIR10
  ↓
SAMtools
  ↓
featureCounts
  ↓
Gene-level raw count matrix
  ↓
DESeq2
  ↓
PCA / differential expression
  ↓
Functional analysis and biological interpretation
```

Reference genome: **TAIR10**  
Gene annotation: **Ensembl Plants release 63 / Araport11**

## Repository Structure

```text
counts/       Gene-level count matrices and summaries
qc/           Quality-control reports
results/      Downstream analysis results
scripts/      Sample metadata and analysis scripts
```

Large raw sequencing files, reference files, and alignment files are excluded from version control.

## Tools

FastQC · MultiQC · HISAT2 · SAMtools · featureCounts · R · DESeq2 · Git

## Data and Attribution Statement

This repository contains an independent computational reanalysis of publicly available sequencing data.

The original experimental design, biological hypotheses, plant material, sample preparation, sequencing, and primary data generation were performed by the authors of the original study. No authorship or ownership of those experimental contributions is claimed here.

Unless otherwise stated, scripts, data processing, statistical analyses, visualizations, and interpretations presented in this repository refer to the independent reanalysis performed in this project.

## References

Fàbregas, N., Lozano-Elena, F., Blasco-Escámez, D. *et al.*  
Overexpression of the vascular brassinosteroid receptor BRL3 confers drought resistance without penalizing plant growth.  
*Nature Communications* **9**, 4680 (2018).  
https://doi.org/10.1038/s41467-018-06861-3

NCBI Gene Expression Omnibus.  
**GSE119382: Transcriptomic study of Arabidopsis roots overexpressing the brassinosteroid receptor BRL3, in control conditions and under severe drought.**  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119382
