#!/usr/bin/env bash
set -euo pipefail

# -- Config ----------------------------------------------------------
BCFTOOLS="${BCFTOOLS:-bcftools}"
PYTHON="${PYTHON:-python3}"

PANEL="data/raw/integrated_call_samples_v3.20130502.ALL.panel"
POPULATION="ASW"

CHR="${1:-22}"
if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be an integer from 1 to 22 (got: $CHR)"
    exit 1
fi

FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
MAP_BASE_URL="https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3"

mkdir -p data/processed data/raw logs

# -- Load modules ------------------------------------------------------
module load bcftools/1.19 2>/dev/null || true

echo "Building ASW sample list..."
awk -v pop="$POPULATION" 'NR>1 && $2==pop {print $1}' "$PANEL" \
    > data/processed/asw.samples
echo "  $(wc -l < data/processed/asw.samples) ASW samples"


RAW_VCF="data/raw/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
mkdir -p "data/processed/chr${CHR}"

# -- Download VCF if not present ------------------------------------------------------
if [[ ! -f "$RAW_VCF" ]]; then
    echo "  Downloading chr${CHR} VCF from 1000G FTP..."
    FNAME="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    wget -q --show-progress \
        -O "$RAW_VCF" \
        "${FTP_BASE}/${FNAME}"
    wget -q \
        -O "${RAW_VCF}.tbi" \
        "${FTP_BASE}/${FNAME}.tbi"
    echo "  Downloaded $RAW_VCF"
else
    echo "  VCF already present: $RAW_VCF"
fi

ASW_VCF="data/processed/chr${CHR}/asw.chr${CHR}.bial.vcf.gz"

# -- Subset to ASW, biallelic SNPs ------------------------------------------------------
if [[ ! -f "$ASW_VCF" ]]; then
    echo "  Subsetting to ASW biallelic SNPs..."
    $BCFTOOLS view \
        -S data/processed/asw.samples \
        -m 2 -M 2 -v snps \
        -O z -o "$ASW_VCF" \
        "$RAW_VCF"
    $BCFTOOLS index -t "$ASW_VCF"
    echo "  Variants: $ASW_VCF"
else
    echo "  ASW VCF already present: $ASW_VCF"
fi

GMAP_DIR="data/raw/maps"
GMAP_FILE="${GMAP_DIR}/genetic_map_chr${CHR}_combined_b37.txt"

# -- Download genetic map ( for germline2) ------------------------------------------------------
mkdir -p "$GMAP_DIR"
if [[ ! -f "$GMAP_FILE" ]]; then
    echo "  Downloading genetic map for chr${CHR}..."
    wget -q --show-progress \
        -O "$GMAP_FILE" \
        "${MAP_BASE_URL}/genetic_map_chr${CHR}_combined_b37.txt" || {
        echo "  ERROR: failed to download $GMAP_FILE"
        echo "  Tried: ${MAP_BASE_URL}/genetic_map_chr${CHR}_combined_b37.txt"
        echo "  Please place the map manually at: $GMAP_FILE"
        exit 1
    }
else
    echo "  Genetic map already exists: $GMAP_FILE"
fi