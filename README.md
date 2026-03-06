# Detecting IBD Segments: GERMLINE vs Beagle

**CSE 284 – Personal Genomics / Bioinformatics, UCSD WI26**
**Authors:** Abishai Ebenezer (A69045190), Inna Amogolonova (A16725376)
**Project Option:** Option 2 – Apply two or more methods and compare on real data

## Overview

This project benchmarks two methods for detecting Identity-by-Descent (IBD) segments in genomic data:

| Tool | Input | Algorithm |
|---|---|---|
| **GERMLINE** | Phased haplotype PED/MAP | Hash-based haplotype matching |
| **Beagle 4.1** | Phased VCF | HMM-based Refined IBD |

We apply both tools to the **1000 Genomes Phase 3 ASW** (African Ancestry in Southwest USA) population on chromosomes 13 and 22, then compare their outputs. We validate against a known mother-child pair (NA20317 → NA20318) confirmed from the 1000G pedigree.

## Repository Structure

```
.
├── data/
│   ├── raw/          # 1000G VCFs, panel, pedigree (large files, not tracked in git)
│   └── processed/    # Per-chromosome ASW VCFs and GERMLINE haploid PEDs
├── results/
│   ├── chr13/        # GERMLINE .match and Beagle .ibd.gz for chr13
│   ├── chr22/        # GERMLINE .match and Beagle .ibd.gz for chr22
│   ├── summary/      # Aggregate TSV metrics
│   └── figures/      # All plots (PNG)
├── scripts/
│   ├── vcf_to_germline_hap.py        # Convert phased VCF → GERMLINE haploid PED
│   ├── summarize_ibd_results.py      # Aggregate and compare IBD outputs
│   ├── plot_analysis.py              # Generate all figures
│   └── extract_parent_child_pairs.py # Extract trios from 1000G pedigree
├── tools/
│   └── germline-master/   # GERMLINE source (compiled by setup.sh)
├── environment.yml    # Conda environment specification
├── setup.sh           # One-time setup: installs tools and conda env
└── run_analysis.sh    # End-to-end pipeline
```

## Quick Start

```bash
# 1. One-time setup (creates conda env, downloads Beagle jar, downloads + compiles GERMLINE)
bash setup.sh

# 2. Activate the conda environment
conda activate cse284-ibd

# 3. Place the two required panel/pedigree files in data/raw/:
#      integrated_call_samples_v3.20130502.ALL.panel
#      integrated_call_samples_v3.20250704.ALL.ped

# 4. Run the full pipeline (downloads chromosome VCFs automatically)
bash run_analysis.sh
```

`setup.sh` handles everything: conda environment creation, Beagle jar download from the Browning lab, and GERMLINE download + compilation from GitHub. The pipeline downloads 1000G chromosome VCFs from the EBI FTP automatically if they are not present.

To run specific chromosomes only:

```bash
CHROMOSOMES="13 22" bash run_analysis.sh
```

## Dataset

- **Source:** 1000 Genomes Project Phase 3 ([internationalgenome.org](https://www.internationalgenome.org/data/))
- **FTP:** `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`
- **Population:** ASW (African Ancestry in Southwest USA) — 61 individuals
- **Chromosomes:** 13 and 22 (biallelic SNPs only)
- **Phasing:** Pre-phased by SHAPEIT2 in the original 1000G release
- **Known parent-child pair:** NA20317 (mother) → NA20318 (child), from the 1000G pedigree

> **Why chromosomes 13 and 22?**
> Both are acrocentric — their short p-arms are not sequenced in 1000G, so the VCF begins after the centromere (~19 Mb for chr13, ~16 Mb for chr22). This avoids a GERMLINE artifact where centromeric high-LD regions on metacentric chromosomes are incorrectly called as IBD.

## Dependencies

| Tool | Version | How obtained |
|---|---|---|
| Beagle | 4.1 (21Jan17.6cc) | Downloaded automatically by `setup.sh` |
| GERMLINE | master | Downloaded and compiled by `setup.sh` |
| bcftools | 1.19 | conda env (`environment.yml`) or `module load bcftools/1.19` |
| Python | 3.8 | conda env |
| matplotlib | 3.5 | conda env |
| Java | 8+ | System (`/usr/bin/java`) or module |

## Pipeline Steps

`run_analysis.sh` runs the following for each chromosome (default: `13 22`):

1. Download 1000G Phase 3 VCF if not present
2. Subset to ASW samples, biallelic SNPs (`bcftools view -m2 -M2 -v snps`)
3. Convert phased VCF to GERMLINE haploid PED (`vcf_to_germline_hap.py`)
4. Run GERMLINE (`-min_m 1`)
5. Run Beagle IBD (`ibdcm=0.5`, `ibdlod=1e-6`)
6. Summarize results across chromosomes (`summarize_ibd_results.py`)
7. Generate figures (`plot_analysis.py`)

All steps are idempotent — re-running skips already-completed steps.

## Key Results

### Chromosome 13 (115 Mb)

| Metric | GERMLINE | Beagle |
|---|---|---|
| Total segments (all pairs) | 517 | 171 |
| Total IBD detected | 773 Mb | 107 Mb |
| Parent-child segments | 26 | 13 |
| Parent-child IBD | 54.3 Mb (47% of chr) | 8.2 Mb |
| Jaccard overlap | 9.4% | — |

### Chromosome 22 (51 Mb)

| Metric | GERMLINE | Beagle |
|---|---|---|
| Total segments (all pairs) | 151 | 27 |
| Total IBD detected | 213 Mb | 17 Mb |
| Parent-child segments | 4 | 2 |
| Parent-child IBD | 6.6 Mb | 1.1 Mb |
| Jaccard overlap | 4.7% | — |

GERMLINE consistently detects more and longer IBD segments than Beagle. The low Jaccard overlap (~5–9%) reflects algorithmic differences: GERMLINE uses exact haplotype matching while Beagle uses probabilistic HMM inference. Running GERMLINE on unphased diploid input produces zero parent-child segments — phased haploid input is required.

## Figures

All figures are written to `results/figures/`:

| Figure | Description |
|---|---|
| `genome_karyogram_NA20317_NA20318.png` | IBD segments for the parent-child pair across chr13 and chr22 |
| `seg_length_hist.png` | Distribution of IBD segment lengths (GERMLINE vs Beagle, both chromosomes) |
| `method_overlap.png` | Total IBD detected per chromosome per method |

## LLM Usage

Claude (Anthropic) was used to assist with Python and shell script development and documentation. The analysis design, method selection, parameter choices, and scientific interpretation were developed by the project authors.
