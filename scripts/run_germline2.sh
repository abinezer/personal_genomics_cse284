#!/usr/bin/env bash
set -euo pipefail

# -- Config --------------------------------
GERMLINE2="${GERMLINE2:-tools/germline2/g2}"
PYTHON="${PYTHON:-python3}"

CHR="${1:-22}"

# -- Prechecks ------------------------------------------------------
if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be an integer from 1 to 22 (got: $CHR)"
    exit 1
fi

if [[ ! -x "$GERMLINE2" ]]; then
    echo "ERROR: germline2 binary not found or not executable at: $GERMLINE2"
    echo "  Run setup.sh to compile germline2"
    exit 1
fi

ASW_VCF="data/processed/chr${CHR}/asw.chr${CHR}.bial.vcf.gz"
if [[ ! -f "$ASW_VCF" ]]; then
    echo "ERROR: ASW VCF not found: $ASW_VCF"
    echo "  Run scripts/prepare_data.sh $CHR first"
    exit 1
fi

GMAP_FILE="data/raw/maps/genetic_map_chr${CHR}_combined_b37.txt"
if [[ ! -f "$GMAP_FILE" ]]; then
    echo "ERROR: genetic map not found: $GMAP_FILE"
    echo "  Run scripts/prepare_data.sh $CHR first"
    exit 1
fi

# -- Setup paths ------------------------------------------------------
HAPS_PREFIX="data/processed/chr${CHR}/asw_shapeit"
OUTDIR="results/chr${CHR}"
GERM2_MATCH="${OUTDIR}/asw_germline2.match"

mkdir -p "$OUTDIR"

# -- Convert phased VCF to SHAPEIT haps/sample format ------------------------------------------------------
if [[ ! -f "${HAPS_PREFIX}.haps" || ! -f "${HAPS_PREFIX}.sample" ]]; then
    echo "Converting phased VCF to SHAPEIT haps/sample format for chr${CHR}..."
    $PYTHON scripts/vcf_to_haps_sample.py \
        "$ASW_VCF" \
        --out "$HAPS_PREFIX" \
        2>&1 | tee "logs/vcf_to_haps_chr${CHR}.log"
else
    echo "SHAPEIT files already present: ${HAPS_PREFIX}.haps/.sample"
fi

# -- Run germline2 ------------------------------------------------------
if [[ ! -f "$GERM2_MATCH" ]]; then
    echo "Running germline2 (haploid mode) for chr${CHR}..."
    "$GERMLINE2" \
        -h \
        -m 1 \
        "${HAPS_PREFIX}.haps" \
        "${HAPS_PREFIX}.sample" \
        "$GMAP_FILE" \
        "$GERM2_MATCH" \
        > "${OUTDIR}/asw_germline2.log" 2>&1
    echo "germline2 segments: $(wc -l < "$GERM2_MATCH")"
else
    echo "germline2 output already present for chr${CHR} ($(wc -l < "$GERM2_MATCH") segments)"
fi

echo "Done: $GERM2_MATCH"
