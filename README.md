# Arabidopsis drought RNA-seq reanalysis

An independent reanalysis of public *Arabidopsis thaliana* root RNA-seq data, looking at the transcriptional response to drought and how BRL3 overexpression changes it.

My main goal was to build and document a reproducible RNA-seq workflow, from raw reads through differential expression to a first biological interpretation.

Related repositories: [Project 2](https://github.com/YuntongMeng/arabidopsis-root-singlecell) · [Project 3](https://github.com/YuntongMeng/arabidopsis-root-multiome). Part of [Arabidopsis root omics](https://github.com/YuntongMeng/arabidopsis-root-omics).

## Key findings

- The full factorial model identified **1,621 drought-responsive genes in WT**
  and **2,257 in BRL3-overexpressing roots** at adjusted p-value < 0.05 and
  |log2 fold change| ≥ 1.
- **137 genes showed a significant genotype × drought interaction**, indicating
  that their drought-associated expression change differed between BRL3 and WT.
- Positive interaction genes were enriched for seven stress-related GO terms,
  including water-deprivation, cold-response, jasmonic-acid, and fatty-acid
  processes. These results provide candidate pathways rather than evidence of a
  BRL3 regulatory mechanism.
- Cut-off-free gene set enrichment analysis (GSEA) recovered all seven of these
  terms and added a broader pattern: cell-cycle, chromatin and
  ribosome-biogenesis gene sets that drought shifted down in WT did not reach
  significance in BRL3.
- Requantifying the reads with Salmon gave essentially the same answer: 94–97%
  of the featureCounts-based DEGs were found again, all in the same direction,
  and gene-level fold changes correlated at Pearson r = 0.93–0.98.

## Analysis overview

The repository contains seven executable scripts. The first three are the
original analysis:

| Script | Purpose | Main validated result |
|---|---|---|
| [`scripts/upstream_rnaseq_workflow.sh`](scripts/upstream_rnaseq_workflow.sh) | Runs FastQC/MultiQC, aligns single-end reads to TAIR10, processes BAM files, and generates reverse-stranded gene counts | Documents the complete 11-sample FASTQ-to-count workflow; the alignment path was smoke-tested on one sample |
| [`scripts/WT_DESeq2_analysis.R`](scripts/WT_DESeq2_analysis.R) | Establishes the downstream workflow using the six WT samples | 1,480 WT drought-responsive genes in the standalone WT model |
| [`scripts/FULL_factorial_DESeq2_analysis.R`](scripts/FULL_factorial_DESeq2_analysis.R) | Fits the 11-sample genotype × condition model and performs all five contrasts, interaction analysis, and GO enrichment | 137 genotype × drought interaction genes; seven enriched positive-interaction GO terms |

The scripts represent two related downstream analyses rather than sequential
steps: the WT-only script establishes the initial workflow, while the full
factorial script answers the main genotype-dependent question.

The other four were added later to check those results with a cut-off-free
enrichment method and with a second upstream workflow:

| Script | Purpose | Main validated result |
|---|---|---|
| [`scripts/FULL_factorial_export_contrasts.R`](scripts/FULL_factorial_export_contrasts.R) | Refits the factorial model and exports the complete WT and BRL3 drought-vs-watered result tables, which the original script did not save | Reproduces the saved interaction table and all three DEG counts exactly |
| [`scripts/GSEA_analysis.R`](scripts/GSEA_analysis.R) | GSEA of GO Biological Process for three contrasts, and comparison with the existing ORA results | 366, 339 and 248 significant terms; all seven interaction ORA terms recovered |
| [`scripts/upstream_salmon_workflow.sh`](scripts/upstream_salmon_workflow.sh) | Quantifies the same FASTQ files with Salmon against the Ensembl Plants release 63 transcriptome, with the genome as decoy | Library type detected as `SR` in all 11 samples; 95.2–96.1% mapping rate |
| [`scripts/salmon_tximport_DESeq2.R`](scripts/salmon_tximport_DESeq2.R) | Summarizes Salmon to gene level with tximport, fits the same DESeq2 model, and compares both workflows | 94–97% of DEGs shared; log2 fold change r = 0.93–0.98 |

The full workflow:

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

The two later checks branch off this workflow:

```text
DESeq2 results for all genes  →  GSEA  →  compared with GO ORA

FASTQ  →  Salmon  →  tximport  →  same DESeq2 model
                                    →  compared with HISAT2 + featureCounts
```

Reference genome: **TAIR10**  
Gene annotation: **Ensembl Plants release 63 / Araport11**  
Transcriptome for Salmon: **Ensembl Plants release 63 cDNA + ncRNA**, the same annotation release

## Inputs and provenance

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

I kept the unequal replicate numbers and modelled them directly in the DESeq2 generalized linear model, rather than dropping or imputing samples.

The exact mapping between sample names, SRA run accessions, genotype, condition,
and replicate is recorded in
[`scripts/sample_metadata.tsv`](scripts/sample_metadata.tsv). The six runs used
for the initial WT-only workflow are also listed in
[`scripts/wt_runs.txt`](scripts/wt_runs.txt).

The tracked downstream input is the generated
[`counts/all_11samples_raw_counts.tsv`](counts/all_11samples_raw_counts.tsv).
Raw FASTQ files, the TAIR10 genome and HISAT2 index, the Ensembl Plants release
63/Araport11 GTF, and BAM files are excluded from Git because of their size.
The same applies to the cDNA/ncRNA FASTA files, the Salmon index and Salmon's
per-transcript `quant.sf` files; the Salmon run logs and summary tables in
[`salmon_quant`](salmon_quant) are tracked.
To rerun from FASTQ, these inputs need to be restored in the paths used by the
upstream scripts.

## Running the analysis

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

I checked the workflow for shell syntax and dependencies, and smoke-tested it by rerunning the FASTQ → HISAT2 → SAMtools path on one sample.

The Salmon workflow uses the same FASTQ files. It also needs the Ensembl Plants
release 63 cDNA and ncRNA FASTA files in `reference/` (the download URLs are in
the script):

```bash
conda activate rnaseq
bash scripts/upstream_salmon_workflow.sh
```

It checks that the transcript FASTA files and the GTF describe exactly the
same transcripts, builds a decoy-aware index, quantifies each sample with
automatic library-type detection, and stops unless every sample is detected as
`SR` (single-end, reverse strand), the Salmon equivalent of featureCounts
`-s 2`. This complete workflow was run on all 11 samples.

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

The later analyses run in this order, after the full factorial script:

```bash
Rscript scripts/FULL_factorial_export_contrasts.R
Rscript scripts/GSEA_analysis.R
Rscript scripts/salmon_tximport_DESeq2.R
```

The first one checks with `stopifnot()` that its refit reproduces the saved
interaction table and DEG counts before it writes anything.
`salmon_tximport_DESeq2.R` needs the Salmon `quant.sf` files, which are not
tracked in Git; its gene-level output is tracked as
[`counts/salmon_11samples_gene_counts.tsv`](counts/salmon_11samples_gene_counts.tsv).

### Software environment

Command-line software versions used for the analysis are recorded in
[`scripts/software_versions.txt`](scripts/software_versions.txt).

The R environment, including versions of DESeq2, ggplot2, clusterProfiler, AnnotationDbi, and other dependencies, is recorded in:
[`scripts/R_sessionInfo.txt`](scripts/R_sessionInfo.txt).

Key software versions include:

- FastQC 0.12.1
- MultiQC 1.35
- HISAT2 2.2.3
- SAMtools 1.24
- featureCounts 2.1.1
- Salmon 2.11.0
- R 4.6.1
- DESeq2 1.52.0
- tximport 1.40.0
- clusterProfiler 4.20.0
- enrichplot 1.32.0
- org.At.tair.db 3.22.0
- UpSetR 1.4.1

## Main analysis

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

A positive interaction means the drought-associated expression change was more positive in BRL3 than in WT. A negative interaction means it was more negative, or less positive, in BRL3. So "positive" and "negative" here are not the same as up- and downregulated.

### Stage 3 — Cut-off-free enrichment (GSEA)

The GO enrichment above is an over-representation analysis (ORA): it only sees
genes that pass the DEG thresholds, and the interaction has just 137 of them.
Gene set enrichment analysis (GSEA) uses every tested gene instead. It ranks
all genes and asks whether the members of a GO term sit together near the top
or the bottom of the ranking, so it can detect pathways whose genes change
modestly but consistently.

Settings:

- ranking metric: the DESeq2 Wald statistic (`stat`, log2 fold change divided
  by its standard error), which keeps the direction and down-weights noisy
  low-count genes; genes with `stat` = NA removed, one value per TAIR ID
- `clusterProfiler::gseGO()` with `org.At.tair.db`, keyType `TAIR`, GO
  Biological Process
- gene-set size 15–500 genes, Benjamini–Hochberg adjustment, adjusted p-value
  < 0.05
- `method = "multilevel"`: the ranking is fixed and gene labels are permuted
  (preranked GSEA), with fgsea's adaptive multilevel Monte Carlo p-values;
  `set.seed(20261116)` before each run
- all three contrasts come from the same factorial model, so the WT and BRL3
  drought tables were exported first with
  `FULL_factorial_export_contrasts.R`

Each contrast tested 1,662 GO terms:

| Contrast | Significant terms | NES > 0 | NES < 0 |
| --- | ---: | ---: | ---: |
| Genotype × drought interaction | 366 | 318 | 48 |
| WT drought vs watered | 339 | 163 | 176 |
| BRL3 drought vs watered | 248 | 195 | 53 |

For the drought contrasts, a positive normalized enrichment score (NES) means
the gene set shifted towards higher expression under drought. For the
interaction, the sign follows the same logic as the interaction log2 fold
change: a positive NES means the drought response of the gene set was more
positive in BRL3 than in WT. It does not mean the gene set was upregulated.

**Comparison with ORA.** ORA results exist for the interaction and for WT
drought vs watered:

| Contrast | Found by both | GSEA only | ORA only |
| --- | ---: | ---: | ---: |
| Genotype × drought interaction | 7 | 359 | 0 |
| WT drought vs watered | 62 | 277 | 19 |

No term was significant in opposite directions in the two methods.

- All seven positive-interaction ORA terms are also significant in GSEA, with
  positive NES.
- Most GSEA-only terms are large gene sets in which many genes shift a little
  in the same direction. Ribosome biogenesis in WT drought is a clear example:
  NES −2.45, although only 1 of its 310 tested genes passes the DEG
  thresholds.
- Of the 19 WT terms found only by ORA, 4 were outside the GSEA size limits
  (for example, response to salt stress has more than 500 genes). The other 15
  are mostly driven by a few strongly changed genes while the rest of the set
  does not move, or moves the other way; a few just miss the GSEA cut-off. "Chemical homeostasis" has 25 downregulated
  DEGs, but its full set of 385 genes has a slightly positive NES.
- The WT ORA used DEGs from the standalone WT-only model, while GSEA ranks
  genes from the factorial model, so the two WT inputs are not identical.

[`GSEA_vs_ORA_comparison.png`](results/GSEA_vs_ORA_comparison.png) shows every
tested term, coloured by which method found it.

### Stage 4 — Second upstream workflow: HISAT2 + featureCounts vs Salmon

To check whether the conclusions depend on the upstream method, I quantified
the same FASTQ files with Salmon and analysed the result with the same DESeq2
model, prefilter and DEG thresholds. Salmon maps reads to transcript sequences
instead of the genome and produces no BAM files. It also keeps multi-mapping
reads and divides them between compatible transcripts, whereas featureCounts
discards them. Both workflows use the same Ensembl Plants release 63
annotation, checked transcript by transcript.

| Contrast | HISAT2 DEGs | Salmon DEGs | Shared | % of HISAT2 DEGs | log2 FC Pearson r |
| --- | ---: | ---: | ---: | ---: | ---: |
| Genotype × drought interaction | 137 | 151 | 129 | 94.2 | 0.932 |
| WT drought vs watered | 1621 | 1614 | 1560 | 96.2 | 0.970 |
| BRL3 drought vs watered | 2257 | 2288 | 2195 | 97.3 | 0.977 |

Every shared DEG has the same direction in both workflows. On average Salmon
assigned 95.7% of reads and featureCounts 94.2%. Most genes that are
significant in only one workflow are threshold effects, where the other
workflow gives a similar estimate just on the other side of a cut-off. A few
genes with many multi-mapping reads, such as *AT2G05250*, get unreliable Salmon
estimates.

Overall, the main results do not depend on the upstream method.

The method comparison, mapping rates, all five DEG counts, the breakdown of
workflow-specific genes and a note on the Salmon TPM table are in
[`docs/upstream_comparison.md`](docs/upstream_comparison.md).

## Biological interpretation

Both wild-type and BRL3-overexpressing plants showed broad transcriptional responses to drought.

More drought-responsive genes passed the thresholds in BRL3 than in WT. That alone does not mean BRL3 plants have a stronger drought response overall, because whether a gene passes also depends on its expression level, effect size, dispersion and the replicate structure.

This is why I used the genotype × drought interaction to find genes whose drought response differed between genotypes.

GO Biological Process enrichment of the **positive interaction genes** identified several stress-related processes, including:

- jasmonic acid metabolic process
- long-chain fatty acid metabolic process
- response to cold
- cold acclimation
- response to water deprivation
- response to water

Genes associated with water-deprivation responses included *KIN1*, *ANAC072*, *AtGolS2*, *LTI30*, *ANAC019*, and *STMP6*.

Jasmonic-acid-related terms were also enriched among the positive interaction genes, which points to a possible link between BRL3-dependent expression differences and stress-related signaling or metabolism.

The negative interaction genes had no significantly enriched GO Biological Process terms at the thresholds used here.

GSEA added a broader pattern that ORA could not see. Gene sets for the mitotic
cell cycle, chromosome and chromatin organization, DNA repair and ribosome
biogenesis had clearly negative NES in WT drought (−1.7 to −2.5), did not reach
significance in BRL3 drought, and had the strongest positive interaction NES
(1.9 to 2.5). Read together: drought shifted these growth-related gene sets
down in WT roots, and this shift was weaker or absent in BRL3. Here the
positive interaction means "less repressed in BRL3", not "induced in BRL3".

On the negative side, the brassinosteroid homeostasis and brassinosteroid
metabolic process gene sets shifted down in BRL3 drought (NES about −2) but not
in WT, and *BRI1* is among the strongest negative interaction genes. Because
BRL3 is itself a brassinosteroid receptor, this is a plausible feedback signal
worth noting, but it is a candidate observation, not evidence of a mechanism.

Overall, the results are consistent with genotype-dependent differences in the drought transcriptional response.

## Interpretation boundary

This is a computational reanalysis of an existing bulk root RNA-seq experiment.
It can identify genotype-associated expression differences and candidate
pathways. It cannot assign those differences to particular root cell types, and
GO enrichment cannot show how BRL3 produces them. Project 2 uses an independent
[normal-root single-cell atlas](https://github.com/YuntongMeng/arabidopsis-root-singlecell)
to add cellular context, but not direct evidence of cell-type-specific drought
effects.

## Limitations

- **BRL3 watered has only two replicates.** The third replicate is missing from
  the public record. Every comparison that uses this group has lower power and
  less stable estimates: BRL3 drought vs watered, BRL3 vs WT under watered
  conditions, and the genotype × drought interaction. DESeq2 shrinks each
  gene's dispersion towards a trend fitted across all genes, which partly
  compensates by borrowing information between genes, but the mean of this
  group still rests on two samples for every gene. One unusual sample can shift
  a whole gene set, so conclusions that involve this group, including the
  broad GSEA interaction signal, are less certain than the WT drought results.
- **GSEA p-values are optimistic.** With two or three replicates per group,
  sample labels cannot be permuted, so gene labels are permuted instead. That
  assumes genes are independent, while genes in a pathway are often
  co-expressed. The hundreds of significant GSEA terms are best read as a
  ranked list of candidate pathways, many of them overlapping GO terms, not as
  hundreds of independent findings.
- **The two WT enrichment analyses use different models.** WT ORA uses DEGs
  from the standalone WT-only model, WT GSEA uses the factorial model.
- **Salmon estimates are less stable for genes with many shared reads.** A few
  genes, such as *AT2G05250*, get extreme fold changes because multi-mapping
  reads are assigned differently between samples.

## Outputs and repository structure

- [`counts`](counts): gene-level count matrices and featureCounts summaries.
  `counts/SRR7779222_s2.txt` is from an early strandedness check: I counted one
  sample with all three featureCounts strandedness settings and chose
  reverse-stranded counting (`-s 2`) from the results. Only that output is kept.
  `counts/salmon_11samples_gene_counts.tsv` holds the Salmon gene-level
  estimated counts.
- [`qc`](qc): generated MultiQC quality-control reports for WT, BRL3, and the
  complete dataset.
- [`salmon_quant`](salmon_quant): Salmon run logs per sample, the mapping-rate
  and library-type summary, the transcript-to-gene table, and the Salmon
  version.
- [`results`](results): differential-expression tables, comparison summaries,
  PCA and distance plots, DEG visualizations, and GO enrichment outputs.
  Files starting with `GSEA_` are the GSEA results, `salmon_` the Salmon-based
  DESeq2 results, and `comparison_` the HISAT2-vs-Salmon comparison.
- [`scripts`](scripts): upstream and downstream workflows, sample metadata, and
  software-environment records.
- [`docs`](docs): the detailed HISAT2-vs-Salmon comparison.

Key outputs include the
[`full comparison summary`](results/FULL_comparison_summary.tsv),
[`137 interaction genes`](results/FULL_interaction_DEGs_padj0.05_log2FC1.tsv),
[`positive-interaction GO results`](results/FULL_interaction_GO_enrichment_positive.tsv),
the [`GSEA summary`](results/GSEA_summary.tsv) and
[`GSEA vs ORA comparison`](results/GSEA_vs_ORA_comparison.tsv), and the
[`HISAT2 vs Salmon DEG overlap`](results/comparison_DEG_overlap.tsv) and
[`fold-change scatter plot`](results/comparison_log2FC_scatter.png).

## Data and attribution

This is an independent computational reanalysis of public sequencing data. The experimental design, hypotheses, plant material, sample preparation, sequencing and the conclusions of the original publication belong to its authors (Fàbregas et al., 2018). The scripts, statistics, figures and interpretations here come from this reanalysis.

I did this as a bioinformatics training and reproducibility exercise, not as an experimental replication of the original study.

## References

Fàbregas, N., Lozano-Elena, F., Blasco-Escámez, D. *et al.*  
Overexpression of the vascular brassinosteroid receptor BRL3 confers drought resistance without penalizing plant growth.  
*Nature Communications* **9**, 4680 (2018).  
[https://doi.org/10.1038/s41467-018-06861-3](https://doi.org/10.1038/s41467-018-06861-3)

NCBI Gene Expression Omnibus.  
**GSE119382: Transcriptomic study of Arabidopsis roots overexpressing the brassinosteroid receptor BRL3, in control conditions and under severe drought.**  
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119382](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119382)
