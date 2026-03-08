# Detecting IBD Segments: germline2 vs Beagle

**CSE 284 – Personal Genomics / Bioinformatics, UCSD WI26**
**Authors:** Abishai Ebenezer (A69045190), Inna Amogolonova (A16725376)
**Project Option:** Option 2 – Apply two or more methods and compare on real data

## Overview

This project benchmarks two methods for detecting Identity-by-Descent (IBD) segments in genomic data:

| Tool | Input | Algorithm |
|---|---|---|
| **germline2** | Phased haplotype HAPS/SAMPLE + genetic map | Hash-based haplotype matching |
| **Beagle 4.1** | Phased VCF | HMM-based Refined IBD |

We apply both tools to the **1000 Genomes Phase 3 ASW** (African Ancestry in Southwest USA) population on chromosomes 21 and 22, then compare their outputs. We validate against a known parent-child pair (NA20317 → NA20318) confirmed from the 1000G pedigree.

## Dependencies 
- Only available for Linux 
- Java 8+ 
- wget 
- conda 
- gcc 
- make 

Please ensure you met the above requirements before running the tool. 

## Instructions to Run 

```bash
# 1. One-time setup (creates conda env, downloads Beagle jar, installs germline2)
bash setup.sh

# 2. Activate the conda environment
conda activate cse284-ibd

# 3. Run the pipeline on chr22 by default 
bash run_full.sh

# 4. Clean up before running the pipeline again (optional)
bash cleanup.sh
```

`setup.sh` handles conda environment creation, Beagle download, metadata download, and fresh germline2 install/compilation. Run the pipeline within the created conda environment.

To run specific chromosomes:

```bash
bash run_full.sh --chromosomes '21 22'
```

To run on all autosomes:

```bash
bash run_full.sh --chromosomes all
```

Helper scripts can also be run separately (after setup and conda environment activation). By default, the commands below run on chr22. Note, the commands must be run in the following order to avoid errors: 
```bash
# 1. Processing data 
bash scripts/prep_data.sh

# 2. Run Germline2
bash scripts/run_germline.sh

# 3. Run Beagle
bash scripts/run_beagle.sh

# 4. Plot fugures and summarize analysis
bash scripts/analyze.sh
```


## Dataset

- **Source:** 1000 Genomes Project Phase 3 ([internationalgenome.org](https://www.internationalgenome.org/data/))
- **FTP:** `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`
- **Population:** ASW (African Ancestry in Southwest USA) — 61 individuals
- **Chromosomes:** 21 and 22 (biallelic SNPs only)
- **Phasing:** Pre-phased by SHAPEIT2 in the original 1000G release
- **Known parent-child pair:** NA20317 (parent) → NA20318 (child), from the 1000G pedigree


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

### Figures

All figures are written to `results/figures/`:

| Figure | Description |
|---|---|
| `genome_karyogram_NA20317_NA20318.png` | IBD segments for the parent-child pair across specified chromosomes |
| `seg_length_hist.png` | Distribution of IBD segment lengths (GERMLINE vs Beagle, specified chromosomes) |
| `method_overlap.png` | Total IBD detected per chromosome per method |

## LLM Usage

Claude (Anthropic) was used to assist with minor Python and shell script development and documentation. The analysis design, method selection, parameter choices, and scientific interpretation were developed by the project authors.

## TODO: 
- [ ] run on different beagle tolerances and pick one 
- [ ] document choosing beagle tolerance 
- [ ] test on full genome (autosomes)
- [ ] edit key findings section 
- [ ] edit README to reflect full genome 