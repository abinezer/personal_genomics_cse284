#!/usr/bin/env bash
# run_analysis.sh - Complete IBD comparison pipeline: GERMLINE vs Beagle
# Usage: bash run_analysis.sh [--skip-data-prep]
#
# Requirements:
#   plink   (PLINK 1.9, must be in PATH)
#   java    (Java 8+)
#   python3 (Python 3.6+, with matplotlib)
#   bcftools (via module load or in PATH)
#   tools/germline-master/germline (compiled from source in tools/)
#   beagle.jar (Beagle 4.1, path set by BEAGLE_JAR below)
#
# Data required in data/raw/:
#   ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#   integrated_call_samples_v3.20130502.ALL.panel
#   integrated_call_samples_v3.20250704.ALL.ped
#
# All paths are relative to the repo root.
set -euo pipefail

# ── Configurable paths ────────────────────────────────────────────────────────
BEAGLE_JAR="${BEAGLE_JAR:-tools/beagle.21Jan17.6cc.jar}"
GERMLINE="${GERMLINE:-tools/germline-master/germline}"
PLINK="${PLINK:-plink}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
PYTHON="${PYTHON:-python3}"

RAW_VCF="data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL="data/raw/integrated_call_samples_v3.20130502.ALL.panel"
PED_FILE="data/raw/integrated_call_samples_v3.20250704.ALL.ped"
CHR_LEN=51244237
POPULATION="ASW"
PARENT="NA20317"
CHILD="NA20318"

mkdir -p data/processed data/raw results/sanity results/summary results/figures logs

# ── Load modules if on cluster ────────────────────────────────────────────────
module load bcftools/1.19 2>/dev/null || true

# ═══════════════════════════════════════════════════════════════════
# STEP 0: Extract parent-child pairs from pedigree
# ═══════════════════════════════════════════════════════════════════
echo "[0] Extracting parent-child pairs..."
$PYTHON scripts/extract_parent_child_pairs.py \
    --ped "$PED_FILE" \
    --out data/processed/parent_child_pairs.tsv

# ═══════════════════════════════════════════════════════════════════
# STEP 1: Subset VCF to ASW population (biallelic SNPs only)
# ═══════════════════════════════════════════════════════════════════
echo "[1] Subsetting VCF to $POPULATION population..."
awk -v pop="$POPULATION" 'NR>1 && $2==pop {print $1}' "$PANEL" > data/processed/asw.samples

$BCFTOOLS view \
    -S data/processed/asw.samples \
    -m 2 -M 2 -v snps \
    -O z -o data/processed/asw.chr22.bial.vcf.gz \
    "$RAW_VCF"

$BCFTOOLS index -t data/processed/asw.chr22.bial.vcf.gz

# ═══════════════════════════════════════════════════════════════════
# STEP 2: Convert to PLINK format for PLINK and GERMLINE (unphased)
# ═══════════════════════════════════════════════════════════════════
echo "[2] Converting to PLINK format..."
$PLINK \
    --vcf data/processed/asw.chr22.bial.vcf.gz \
    --recode \
    --out data/processed/asw_bial \
    --allow-extra-chr --chr 22 \
    --biallelic-only strict \
    --geno 0.05 \
    --maf 0.01 \
    --const-fid 0 \
    --double-id \
    2>&1 | tee logs/plink_recode.log

# Fill missing genotypes (required by GERMLINE)
$PLINK \
    --file data/processed/asw_bial \
    --recode \
    --fill-missing-a2 \
    --out data/processed/asw_bial_filled \
    2>&1 | tee logs/plink_fill.log

# ═══════════════════════════════════════════════════════════════════
# STEP 3: PLINK genome (relatedness baseline)
# ═══════════════════════════════════════════════════════════════════
echo "[3] Running PLINK --genome..."
$PLINK \
    --bfile data/processed/asw_bial \
    --genome \
    --out results/asw_plink \
    2>&1 | tee logs/plink_genome.log

# ═══════════════════════════════════════════════════════════════════
# STEP 4: Create phased haplotype input for GERMLINE
# ═══════════════════════════════════════════════════════════════════
echo "[4] Converting phased VCF to GERMLINE haploid format..."
$PYTHON scripts/vcf_to_germline_hap.py \
    data/processed/asw.chr22.bial.vcf.gz \
    --out data/processed/asw_hap \
    2>&1 | tee logs/vcf_to_hap.log

# ═══════════════════════════════════════════════════════════════════
# STEP 5: Run GERMLINE (phased haploid mode)
# ═══════════════════════════════════════════════════════════════════
echo "[5] Running GERMLINE (phased)..."
$GERMLINE \
    -input data/processed/asw_hap.ped data/processed/asw_hap.map \
    -output results/asw_germline_phased \
    -min_m 1 \
    > results/asw_germline_phased.log 2>&1
echo "  GERMLINE segments: $(wc -l < results/asw_germline_phased.match)"

# ═══════════════════════════════════════════════════════════════════
# STEP 6: Run Beagle IBD detection
# ═══════════════════════════════════════════════════════════════════
echo "[6] Running Beagle IBD detection..."
java -Xmx4g -jar "$BEAGLE_JAR" \
    gt=data/processed/asw.chr22.bial.vcf.gz \
    out=results/asw_beagle \
    ibd=true \
    ibdcm=0.5 \
    ibdlod=1e-6 \
    nthreads=4 \
    > results/asw_beagle.log 2>&1
echo "  Beagle segments: $(python3 -c "import gzip; d=gzip.open('results/asw_beagle.ibd.gz','rb').read(); print(d.count(b'\n'))")"

# ═══════════════════════════════════════════════════════════════════
# STEP 7: Sanity check - validate parent-child pair on chr22
# ═══════════════════════════════════════════════════════════════════
echo "[7] Sanity check: parent-child ($PARENT vs $CHILD)..."
echo -e "$PARENT\n$CHILD" > data/processed/sanity_pair.samples
$BCFTOOLS view \
    -S data/processed/sanity_pair.samples \
    -O z -o data/processed/sanity_pair.chr22.vcf.gz \
    data/processed/asw.chr22.bial.vcf.gz
$BCFTOOLS index -t data/processed/sanity_pair.chr22.vcf.gz

# Sanity PLINK
$PLINK \
    --vcf data/processed/sanity_pair.chr22.vcf.gz \
    --genome --out results/sanity/plink_sanity \
    --double-id --const-fid 0 \
    2>&1 | tee logs/plink_sanity.log

# Sanity GERMLINE
$PYTHON scripts/vcf_to_germline_hap.py \
    data/processed/sanity_pair.chr22.vcf.gz \
    --out data/processed/sanity_hap
$GERMLINE \
    -input data/processed/sanity_hap.ped data/processed/sanity_hap.map \
    -output results/sanity/germline_sanity_phased \
    -min_m 1 \
    > results/sanity/germline_sanity_phased.log 2>&1

# Sanity Beagle
java -Xmx2g -jar "$BEAGLE_JAR" \
    gt=data/processed/sanity_pair.chr22.vcf.gz \
    out=results/sanity/beagle_sanity \
    ibd=true ibdcm=0.5 ibdlod=1e-6 nthreads=4 \
    > results/sanity/beagle_sanity.log 2>&1

echo "  Sanity PLINK PI_HAT:"
awk 'NR>1 {print "  Z0="$7,"Z1="$8,"Z2="$9,"PI_HAT="$10}' results/sanity/plink_sanity.genome
echo "  Sanity GERMLINE segments: $(wc -l < results/sanity/germline_sanity_phased.match)"

# ═══════════════════════════════════════════════════════════════════
# STEP 8: Summarize results
# ═══════════════════════════════════════════════════════════════════
echo "[8] Summarizing results..."
$PYTHON scripts/summarize_ibd_results.py \
    --beagle-ibd results/asw_beagle.ibd.gz \
    --germline-match results/asw_germline_phased.match \
    --plink-genome results/asw_plink.genome \
    --parent "$PARENT" \
    --child "$CHILD" \
    --outdir results/summary

# ═══════════════════════════════════════════════════════════════════
# STEP 9: Generate figures
# ═══════════════════════════════════════════════════════════════════
echo "[9] Generating figures..."
$PYTHON scripts/plot_analysis.py \
    --beagle-ibd results/asw_beagle.ibd.gz \
    --germline-match results/asw_germline_phased.match \
    --plink-genome results/asw_plink.genome \
    --overall-metrics results/summary/overall_metrics.tsv \
    --pairs "NA20317:NA20318" "NA20359:NA20362" "NA20320:NA20321" \
    --outdir results/figures

echo ""
echo "=== Pipeline complete ==="
echo "  Results: results/"
echo "  Figures: results/figures/"
echo "  Summary: results/summary/"
