# Detecting IBD Segments: GERMLINE vs Beagle

**CSE 284 – Personal Genomics / Bioinformatics, UCSD WI26**
**Authors:** Abishai Ebenezer (A69045190), Inna Amogolonova (A16725376)
**Project Option:** Option 2 – Apply two or more methods for a task discussed in class and compare results on real data

## Overview

This project benchmarks two widely-used methods for detecting Identity-by-Descent (IBD) segments in genomic data:

| Tool | Input | Algorithm |
|---|---|---|
| **GERMLINE** | Phased haplotype PED/MAP | Hash-based haplotype matching |
| **Beagle 4.1** | Phased VCF | HMM-based Refined IBD |

We apply both tools to the **1000 Genomes Phase 3 ASW** (African Ancestry in Southwest USA) population on chromosome 22, then compare their outputs against each other and against **PLINK `--genome`** as a relatedness baseline.

## Repository Structure

```
.
├── data/
│   ├── raw/                     # 1000G VCF, panel, and pedigree files
│   └── processed/               # Filtered VCF, PLINK PED/MAP, haploid HAP/MAP
├── results/
│   ├── asw_germline_phased.match  # GERMLINE IBD segments (phased)
│   ├── asw_beagle.ibd.gz          # Beagle IBD segments
│   ├── asw_plink.genome           # PLINK pairwise relatedness
│   ├── sanity/                    # Parent-child validation results
│   ├── summary/                   # Aggregated metrics (TSV files)
│   └── figures/                   # All plots (PNG)
├── scripts/
│   ├── extract_parent_child_pairs.py  # Extract trios from 1000G pedigree
│   ├── vcf_to_germline_hap.py         # Convert phased VCF → GERMLINE haploid PED
│   ├── summarize_ibd_results.py       # Aggregate and compare IBD outputs
│   ├── plot_analysis.py               # Generate all figures
│   └── plot_pair_karyogram.py         # Per-pair chromosome karyogram
├── tools/
│   └── germline-master/           # GERMLINE source and compiled binary
├── run_analysis.sh                # End-to-end pipeline script
└── README.md
```

## Dependencies

| Tool | Version | Notes |
|---|---|---|
| PLINK | 1.9 | `plink` in PATH |
| Beagle | 4.1 (21Jan17.6cc) | `tools/beagle.21Jan17.6cc.jar` |
| GERMLINE | 1.5.1 (2019-04-04) | Compiled binary at `tools/germline-master/germline` |
| bcftools | 1.19 | via `module load bcftools/1.19` |
| Python | 3.8+ | `matplotlib` required |
| Java | 8+ | For Beagle |

## Dataset

- **Source:** 1000 Genomes Project Phase 3 ([ftp.1000genomes.ebi.ac.uk](http://ftp.1000genomes.ebi.ac.uk))
- **File:** `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- **Population:** ASW (African Ancestry in Southwest USA) — 61 individuals
- **Chromosome:** 22 only (biallelic SNPs, MAF ≥ 0.01)
- **Markers after filtering:** ~1,055,454 SNPs
- **Known parent-child pair:** NA20317 (mother) → NA20318 (child)

> **Note on phasing:** The 1000G Phase 3 VCFs are already phased (using SHAPEIT2). This phased information is used by both tools. GERMLINE is given phased haploid input via `vcf_to_germline_hap.py`, and Beagle uses the phased genotypes directly.

## Installation

### Compile GERMLINE

```bash
cd tools
tar -xzf germline.tar.gz
cd germline-master
make
```

### Download Beagle 4.1

```bash
mkdir -p tools
wget -O tools/beagle.21Jan17.6cc.jar \
  https://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar
```

### Python dependencies

```bash
pip install --user matplotlib
```

## Running the Analysis

```bash
# Full pipeline (requires PLINK in PATH and Beagle jar)
bash run_analysis.sh

# Or run individual steps:
# 1. Subset VCF to ASW
bcftools view -S data/processed/asw.samples -m2 -M2 -v snps \
    -Oz -o data/processed/asw.chr22.bial.vcf.gz \
    data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# 2. Convert phased VCF to GERMLINE haploid format
python3 scripts/vcf_to_germline_hap.py \
    data/processed/asw.chr22.bial.vcf.gz \
    --out data/processed/asw_hap

# 3. Run GERMLINE (phased haploid mode)
tools/germline-master/germline \
    -input data/processed/asw_hap.ped data/processed/asw_hap.map \
    -output results/asw_germline_phased \
    -min_m 1

# 4. Run Beagle IBD detection
java -Xmx4g -jar tools/beagle.21Jan17.6cc.jar \
    gt=data/processed/asw.chr22.bial.vcf.gz \
    out=results/asw_beagle \
    ibd=true ibdcm=0.5 ibdlod=1e-6 nthreads=4

# 5. PLINK relatedness baseline
plink --bfile data/processed/asw_bial --genome --out results/asw_plink

# 6. Summarize and plot
python3 scripts/summarize_ibd_results.py \
    --beagle-ibd results/asw_beagle.ibd.gz \
    --germline-match results/asw_germline_phased.match \
    --plink-genome results/asw_plink.genome \
    --outdir results/summary

python3 scripts/plot_analysis.py \
    --beagle-ibd results/asw_beagle.ibd.gz \
    --germline-match results/asw_germline_phased.match \
    --plink-genome results/asw_plink.genome \
    --overall-metrics results/summary/overall_metrics.tsv \
    --pairs "NA20317:NA20318" "NA20359:NA20362" "NA20320:NA20321" \
    --outdir results/figures
```

## Small Test Example (Sanity Check)

To verify the pipeline on a known parent-child pair (NA20317 / NA20318) using chromosome 22 only:

```bash
# Subset to just the parent-child pair
echo -e "NA20317\nNA20318" > data/processed/sanity_pair.samples
bcftools view -S data/processed/sanity_pair.samples \
    -Oz -o data/processed/sanity_pair.chr22.vcf.gz \
    data/processed/asw.chr22.bial.vcf.gz

# Run PLINK sanity check
plink --vcf data/processed/sanity_pair.chr22.vcf.gz \
    --genome --double-id --const-fid 0 \
    --out results/sanity/plink_sanity

# Run GERMLINE sanity check (phased)
python3 scripts/vcf_to_germline_hap.py \
    data/processed/sanity_pair.chr22.vcf.gz \
    --out data/processed/sanity_hap
tools/germline-master/germline \
    -input data/processed/sanity_hap.ped data/processed/sanity_hap.map \
    -output results/sanity/germline_sanity_phased -min_m 1

# Run Beagle sanity check
java -Xmx2g -jar tools/beagle.21Jan17.6cc.jar \
    gt=data/processed/sanity_pair.chr22.vcf.gz \
    out=results/sanity/beagle_sanity \
    ibd=true ibdcm=0.5 ibdlod=1e-6 nthreads=4
```

**Expected results for NA20317 / NA20318 (mother-child):**
- PLINK PI_HAT ≈ 0.5 (first-degree relative; note: may be inflated in ASW due to population structure)
- Both GERMLINE and Beagle should detect multiple IBD segments across chr22

## Key Results (Chr22, ASW Population)

| Metric | GERMLINE (phased) | Beagle |
|---|---|---|
| Total IBD segments | 151 | 27 |
| Total IBD bp | 212.6 Mb | 16.6 Mb |
| Parent-child segments detected | 4 | 2 |
| Jaccard overlap (bp) | 4.7% | — |

GERMLINE (phased haploid mode) detects more segments and more total IBD than Beagle. The low Jaccard overlap (4.7%) reflects differences in the algorithms: GERMLINE looks for exact haplotype matches while Beagle uses probabilistic HMM inference.

## Output Files

| File | Description |
|---|---|
| `results/asw_germline_phased.match` | GERMLINE haplotype-level IBD matches (tab-separated) |
| `results/asw_beagle.ibd.gz` | Beagle IBD segments (haplotype-level) |
| `results/asw_plink.genome` | PLINK pairwise IBD coefficients (Z0, Z1, Z2, PI_HAT) |
| `results/summary/overall_metrics.tsv` | Aggregate segment counts and Jaccard similarity |
| `results/summary/pair_metrics.tsv` | Per-pair segment counts, bp, and PI_HAT |
| `results/summary/parent_child_check.tsv` | Validation on the known parent-child pair |
| `results/figures/seg_length_hist.png` | Segment length distribution comparison |
| `results/figures/pihat_vs_ibd.png` | PI_HAT vs chr22 IBD sharing scatter |
| `results/figures/method_overlap.png` | Total IBD bp by method (Beagle-only, shared, GERMLINE-only) |
| `results/figures/karyogram_*.png` | Per-pair chromosome 22 IBD maps |

## LLM Usage

Claude (Anthropic) was used to assist with Python script development and documentation editing. The analysis design, interpretation, and core methods were developed by the project authors.
