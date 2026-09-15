# Arabidopsis Drought RNA-seq Reanalysis

An independent reanalysis of publicly available *Arabidopsis thaliana* RNA-seq data to investigate transcriptional responses to drought and the effects of BRL3 overexpression.

The main goal of this project was to build and document a reproducible RNA-seq workflow from raw sequencing reads to differential expression analysis and preliminary biological interpretation.

## Dataset

RNA-seq data were obtained from **NCBI Gene Expression Omnibus (GEO), accession GSE119382**.

The original experiment used a factorial design comparing:

- wild type (Col-0) and 35S:BRL3 plants
- watered and drought conditions
- biological replicates from *Arabidopsis thaliana* roots

The public GEO record currently contains 11 RNA-seq samples:

- WT watered: 3 replicates
- WT drought: 3 replicates
- BRL3 watered: 2 replicates
- BRL3 drought: 3 replicates

The unequal number of replicates was retained in the analysis and handled directly using the DESeq2 generalized linear model rather than through sample removal or imputation.

## Analysis

The analysis was performed in two stages.

### Stage 1 — Wild-type drought response

Six wild-type samples were first used to establish the RNA-seq workflow:

- 3 drought biological replicates
- 3 watered biological replicates

Differential expression was evaluated using DESeq2 with watered samples as the reference condition.

Using:

- adjusted p-value < 0.05
- |log2 fold change| ≥ 1

the standalone WT-only analysis identified **1,480 drought-responsive genes**, including:

- 904 upregulated genes
- 576 downregulated genes

PCA and sample-distance analysis showed clear separation between watered and drought-treated wild-type samples.

GO Biological Process enrichment highlighted processes associated with water deprivation, abiotic stress responses, and related metabolic processes.

### Stage 2 — Full genotype × condition analysis

The full set of 11 public samples was analyzed using the DESeq2 factorial model:

```r
~ genotype * condition
```

This design allows genotype, drought treatment, and genotype × condition interaction effects to be evaluated within the same statistical model.

The following comparisons were examined using the same DEG threshold:

| Comparison                          | Total DEGs | Positive / Up | Negative / Down |
| ----------------------------------- | ---------: | ------------: | --------------: |
| WT drought vs watered               |       1621 |           986 |             635 |
| BRL3 drought vs watered             |       2257 |          1337 |             920 |
| BRL3 vs WT under watered conditions |        539 |           105 |             434 |
| BRL3 vs WT under drought conditions |        591 |           196 |             395 |
| Genotype × drought interaction      |        137 |            57 |              80 |

The WT drought-vs-watered contrast yields **1,621 DEGs in the full factorial model**, compared with **1,480 DEGs in the standalone WT-only analysis**. The biological contrast is the same, but the two analyses use different DESeq2 model structures and dispersion-estimation frameworks, which can change which genes pass the statistical thresholds.

The first four comparisons describe direct expression differences between experimental groups.

The interaction term addresses a different question:

> Does the transcriptional response to drought differ between BRL3-overexpressing and wild-type plants?

The interaction corresponds to:

```text
(BRL3 drought − BRL3 watered)
−
(WT drought − WT watered)
```

Using adjusted p-value < 0.05 and |log2 fold change| ≥ 1, **137 genes showed significant genotype × drought interaction effects**.

Among these:

- 57 showed positive interaction effects
- 80 showed negative interaction effects

A positive interaction indicates that the drought-associated expression change was more positive in BRL3 than in WT, whereas a negative interaction indicates a more negative or less positive drought response in BRL3.

These terms should not be interpreted simply as gene upregulation or downregulation.

## Biological Interpretation

The reanalysis showed broad transcriptional responses to drought in both wild-type and BRL3-overexpressing plants.

Under the thresholds used here, more drought-responsive genes were detected in BRL3 than in WT. However, DEG counts alone should not be interpreted as evidence that one genotype has a globally stronger drought response, because statistical detection also depends on expression level, effect size, dispersion, and replicate structure.

The genotype × drought interaction analysis was therefore used to identify genes whose drought response differed between genotypes.

GO Biological Process enrichment of the **positive interaction genes** identified several stress-related processes, including:

- jasmonic acid metabolic process
- long-chain fatty acid metabolic process
- response to cold
- cold acclimation
- response to water deprivation
- response to water

Genes associated with water-deprivation responses included *KIN1*, *ANAC072*, *AtGolS2*, *LTI30*, *ANAC019*, and *STMP6*.

Jasmonic-acid-related enrichment was also observed among the positive interaction genes. This provides a possible link between BRL3-dependent transcriptional differences and stress-associated signaling or metabolism, but the enrichment result alone is not sufficient to establish a regulatory mechanism.

No significantly enriched GO Biological Process terms were detected for the negative interaction gene set under the enrichment thresholds used in this analysis.

Overall, these results are consistent with genotype-dependent differences in the drought transcriptional response. The biological interpretation presented here is intended as an exploratory interpretation of the reanalyzed RNA-seq data rather than as mechanistic validation of BRL3 function.

## Workflow

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
PCA / clustering / differential expression
  ↓
Genotype × condition interaction analysis
  ↓
GO enrichment and biological interpretation
```

Reference genome: **TAIR10**  
Gene annotation: **Ensembl Plants release 63 / Araport11**

## Reproducing the Analysis

The repository includes executable scripts documenting both upstream RNA-seq processing and downstream statistical analysis.

### Upstream processing

Given the raw FASTQ files and prepared TAIR10 reference files, the upstream workflow can be run with:

```bash
conda activate rnaseq
bash scripts/upstream_rnaseq_workflow.sh
```

The upstream workflow performs:

- FastQC quality assessment
- MultiQC report generation
- HISAT2 alignment
- SAMtools sorting and BAM processing
- reverse-stranded gene-level counting with featureCounts
- generation of the 11-sample raw count matrix

The workflow was checked using shell syntax and dependency validation and was additionally smoke-tested by rerunning the FASTQ → HISAT2 → SAMtools path on one sample.

Large raw FASTQ files, genome reference files, HISAT2 indices, and BAM files are excluded from version control.

### Downstream analysis

The standalone WT analysis is implemented in:

```text
scripts/WT_DESeq2_analysis.R
```

The full genotype × condition analysis is implemented in:

```text
scripts/FULL_factorial_DESeq2_analysis.R
```

The full analysis can be rerun from the project root with:

```bash
Rscript scripts/FULL_factorial_DESeq2_analysis.R
```

The full factorial script calculates all five reported comparisons directly from the fitted DESeq2 model, including the WT drought-vs-watered contrast. Summary values are generated programmatically rather than manually entered.

### Software environment

Command-line software versions used for the analysis are recorded in:

```text
scripts/software_versions.txt
```

The R environment, including versions of DESeq2, ggplot2, clusterProfiler, AnnotationDbi, and other dependencies, is recorded in:

```text
scripts/R_sessionInfo.txt
```

Key software versions include:

- FastQC 0.12.1
- MultiQC 1.35
- HISAT2 2.2.3
- SAMtools 1.24
- featureCounts 2.1.1
- R 4.6.1
- DESeq2 1.52.0
- clusterProfiler 4.20.0
- org.At.tair.db 3.22.0

## Repository Structure

```text
counts/       Gene-level count matrices and summaries
qc/           Quality-control reports
results/      Differential expression, figures, and enrichment results
scripts/      Workflow scripts, analysis scripts, and environment records
```

Large raw sequencing files, reference files, and alignment files are excluded from version control.

The downstream analysis is organized into separate scripts for the initial wild-type analysis and the full factorial analysis.

## Tools

FastQC · MultiQC · HISAT2 · SAMtools · featureCounts · R · DESeq2 · ggplot2 · clusterProfiler · Git

## Data and Attribution Statement

This repository contains an independent computational reanalysis of publicly available sequencing data.

The original experimental design, biological hypotheses, plant material, sample preparation, sequencing, and primary data generation were performed by the authors of the original study. No authorship or ownership of those experimental contributions is claimed here.

The biological conclusions reported in the original publication should likewise be attributed to the original authors.

Unless otherwise stated, scripts, data processing, statistical analyses, visualizations, and interpretations presented in this repository refer specifically to the independent computational reanalysis performed in this project.

This project should therefore be interpreted as a bioinformatics training and reproducibility exercise using an existing public dataset rather than as an independent experimental replication of the original study.

## References

Fàbregas, N., Lozano-Elena, F., Blasco-Escámez, D. *et al.*  
Overexpression of the vascular brassinosteroid receptor BRL3 confers drought resistance without penalizing plant growth.  
*Nature Communications* **9**, 4680 (2018).  
[https://doi.org/10.1038/s41467-018-06861-3](https://doi.org/10.1038/s41467-018-06861-3)

NCBI Gene Expression Omnibus.  
**GSE119382: Transcriptomic study of Arabidopsis roots overexpressing the brassinosteroid receptor BRL3, in control conditions and under severe drought.**  
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119382](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119382)