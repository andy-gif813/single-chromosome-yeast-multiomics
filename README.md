# Appliedbio_ABCC_group11

Reproducing the genomic and transcriptomic analyses of the single-chromosome yeast strain **SY14** versus its parental wild-type strain **BY4742**, based on:

> Shao, Y. et al. (2018) *Creating a functional single-chromosome yeast.* **Nature** 560, 331–335.

This repository contains the scripts, intermediate tables and figures produced by Group 11 for the Applied Bioinformatics (ABCC) group project.

---

## Project overview

SY14 was engineered by fusing the 16 native chromosomes of *Saccharomyces cerevisiae* BY4742 into a single chromosome. Our aims were to:

1. **PacBio long reads** – assemble the SY14 genome de novo and confirm the chromosome fusion by whole-genome alignment to the reference.
2. **Illumina DNA-seq** – map short reads to the SY14 assembly / WT genome and identify small insertions and deletions (indels).
3. **RNA-seq** – quantify gene expression in SY14 vs BY4742 (3 biological replicates each) and identify differentially expressed genes (DEGs).

```
PacBio subreads ──► Flye assembly ──► LAST / minimap2 vs reference ──► fusion & large-deletion analysis
Illumina DNA    ──► FastQC ──► BWA-MEM ──► bcftools mpileup/call ──► indel filtering
RNA-seq         ──► FastQC ──► Trimmomatic ──► Bowtie2 / RSEM ──► DESeq2 ──► DEGs, heatmap, correlation
```

---

## Data

| Data type | Accession | Samples |
|---|---|---|
| PacBio subreads (SY14) | SRR6823435 | 1 |
| Illumina WGS (SY14) | SRR6825081 | 1 (paired-end) |
| RNA-seq | SRR7059704–SRR7059709 | BY4742 ×3, SY14 ×3 |

RNA-seq sample mapping:

| Sample | SRR | Condition |
|---|---|---|
| BY4742_R1 | SRR7059707 | wild_type |
| BY4742_R2 | SRR7059706 | wild_type |
| BY4742_R3 | SRR7059709 | wild_type |
| SY14_R1 | SRR7059708 | engineered |
| SY14_R2 | SRR7059705 | engineered |
| SY14_R3 | SRR7059704 | engineered |

Reference genome: *S. cerevisiae* S288C R64 (used for RNA-seq quantification and genome comparison).
Raw data are **not** included in this repository; download them from NCBI SRA (e.g. `prefetch` + `fasterq-dump`).

---

## Repository structure

```
.
├── pacbio/
│   ├── flye_full_dataset_fixed.slurm   # de novo assembly of SY14 with Flye
│   ├── pacbio_compare_SY14_WT.sh       # minimap2 (map-pb) alignment of subreads to reference
│   ├── last_alignment.slurm            # LAST whole-genome alignment, SY14 assembly vs BY4742/S288C
│   ├── analyze_indel.py                # detect large candidate deletions from .coords alignments
│   ├── visualize_fusion.py             # plot chromosome-fusion evidence from alignment links
│   └── chromosome_fusion_visualization.png
├── DNA_seq/
│   ├── fastqc/                         # FastQC reports for SRR6825081
│   ├── SY14Illunina_mapping_SY14assemble.sh  # BWA-MEM mapping, flagstat, depth
│   ├── detect_indel.sh                 # bcftools mpileup + call, keep indels
│   └── filter_indel.sh                 # filter and count indels
└── RNA_seq/
    ├── data_prepare/
    │   ├── fastqc.sh                   # read QC
    │   ├── trim.sh                     # Trimmomatic adapter/quality trimming
    │   ├── submit_bowtie2.sh           # wrapper to submit the Bowtie2 job
    │   ├── bowtie_process_all_samples.sh  # Bowtie2 alignment for all 6 samples
    │   ├── rsem_run_all_trimmed.sh     # RSEM quantification (trimmed reads)
    │   └── rsem_run_all_untrimmed.sh   # RSEM quantification (untrimmed reads, for comparison)
    ├── DESeq/
    │   ├── DESeq2.Rmd / DESeq2.html    # differential expression analysis (code + rendered report)
    │   ├── sample_info.txt             # sample metadata
    │   ├── count_matrix_raw.csv        # RSEM expected counts
    │   ├── count_matrix_tpm.csv        # TPM matrix
    │   ├── DESeq2_results.csv          # full DESeq2 output
    │   ├── Verified_ORF_results.csv    # significant DEGs restricted to SGD Verified ORFs
    │   ├── ORF.csv                     # SGD Verified ORF list
    │   ├── pearson.csv                 # sample–sample Pearson correlation
    │   └── result.png                  # heatmap of DEGs
    ├── Gene_count_ggplot.R             # stacked bar chart of DEG categories
    └── Gene_count_ggplot.png
```

---

## Software and environment

All shell scripts were written for the King's College London **CREATE HPC** (SLURM scheduler, `module load` and conda environments).

| Step | Tool (version) |
|---|---|
| Read QC | FastQC 0.12.1 |
| Trimming | Trimmomatic 0.39 |
| Long-read assembly | Flye |
| Long-read alignment | minimap2, LAST |
| Short-read DNA mapping | BWA-MEM, SAMtools 1.17 |
| Variant calling | bcftools 1.12, htslib 1.12 |
| RNA-seq alignment / quantification | Bowtie2 2.5.1, RSEM |
| Differential expression | R, DESeq2, pheatmap, RColorBrewer, ggplot2 |
| Plotting (Python) | matplotlib, seaborn, pandas, numpy |

---

## How to run

> ⚠️ Scripts contain absolute paths on the HPC (e.g. `/scratch/grp/msc_appbio/Group11_ABCC/...`). Edit the path variables at the top of each script before running.

### 1. PacBio – assembly and fusion analysis
```bash
cd pacbio
sbatch flye_full_dataset_fixed.slurm        # assembly → 04_results/flye_assembly_full/assembly.fasta
sbatch last_alignment.slurm                 # SY14 assembly vs reference
sbatch pacbio_compare_SY14_WT.sh            # read-level alignment with minimap2
python3 analyze_indel.py <alignment.coords> # candidate large deletions (500 bp – 200 kb)
python3 visualize_fusion.py                 # needs circos_links.txt
```

### 2. Illumina DNA-seq – indel detection
```bash
cd DNA_seq
sbatch SY14Illunina_mapping_SY14assemble.sh
sbatch detect_indel.sh
sbatch filter_indel.sh
```
Indel filter: `QUAL ≥ 30`, `10 ≤ FORMAT/DP ≤ 500`, `IMF ≥ 0.3`.

### 3. RNA-seq – differential expression
```bash
cd RNA_seq/data_prepare
sbatch fastqc.sh
sbatch trim.sh            # ILLUMINACLIP:TruSeq3-PE:2:30:10, SLIDINGWINDOW:4:20, MINLEN:36
bash   submit_bowtie2.sh
bash   rsem_run_all_trimmed.sh
```
Then open `RNA_seq/DESeq/DESeq2.Rmd` in RStudio and knit.
DEG thresholds: `padj < 0.001` and `|log2FoldChange| ≥ 1`, with BY4742 (wild_type) as the reference level; genes with zero total counts are removed.

---

## Results

<!-- Fill in your own numbers / one-sentence conclusions -->

**Chromosome fusion (PacBio)**

![Chromosome fusion](pacbio/chromosome_fusion_visualization.png)

- Number of contigs in the SY14 assembly: _TODO_
- Assembly size / N50: _TODO_
- Evidence that all 16 chromosomes are joined into one: _TODO_

**Indels (Illumina)**

- Raw indels called: _TODO_; after filtering: _TODO_

**Differential expression (RNA-seq)**

![DEG heatmap](RNA_seq/DESeq/result.png)
![DEG categories](RNA_seq/Gene_count_ggplot.png)

- Significant DEGs among Verified ORFs: 10 up-regulated, 43 down-regulated
  - TPE (telomere-proximal) genes: 1 up / 21 down
  - Stress-response genes: 2 up / 1 down
  - Others: 7 up / 21 down
- Replicates are highly correlated (see `pearson.csv`).
- Comparison with the original paper: _TODO_

---

## Team members

- Andy Chou
- Mengyue Liu
- Yifei Zhang

All analyses (PacBio, Illumina DNA-seq and RNA-seq) were carried out jointly by all three members.

---

## Reference

Shao, Y., Lu, N., Wu, Z. *et al.* Creating a functional single-chromosome yeast. *Nature* **560**, 331–335 (2018). https://doi.org/10.1038/s41586-018-0382-x
