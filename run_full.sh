#!/usr/bin/env bash
# run_full.sh — Orchestrate full IBD pipeline: prep, germline2, Beagle, analysis
#
# Usage:
#   bash run_full.sh                          # default: chr22 only
#   bash run_full.sh --chromosomes "22"       # explicit single chromosome
#   bash run_full.sh --chromosomes "21 22"    # multiple chromosomes
#   bash run_full.sh --chromosomes all        # genome-wide (autosomes 1-22)
#
# Requires:
#   - setup.sh run once (for tools)
#   - Conda environment activated
set -euo pipefail

# -- Config ----------------------------------------------------------------------
PED_FILE="data/raw/integrated_call_samples_v3.20250704.ALL.ped"
PARENT="NA20317"
CHILD="NA20318"

# -- Parse arguments ----------------------------------------------------------------------
CHROMOSOMES="22"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --chromosomes)
            CHROMOSOMES="$2"
            shift 2
            ;;
        *)
            echo "ERROR: unknown argument: $1"
            echo "Usage: bash run_full.sh [--chromosomes 'all'|'21 22'|'22']"
            exit 1
            ;;
    esac
done

# -- Expand "all" to autosomes 1-22 ----------------------------------------------------------------------
if [[ "$CHROMOSOMES" == "all" ]]; then
    CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
fi

# -- Module setup ----------------------------------------------------------------------
module load bcftools/1.19 2>/dev/null || true

# -- Create output directories ----------------------------------------------------------------------
mkdir -p data/processed data/raw results/summary results/figures logs

# -- Info: extract parent-child pairs ----------------------------------------------------------------------
echo "[0] Extracting parent-child pairs from pedigree..."
python3 scripts/extract_parent_child_pairs.py \
    --ped "$PED_FILE" \
    --out data/processed/parent_child_pairs.tsv


# -- Steps 1-3: prep data, run germline2, run beagle on every chr ----------------------------------------------------------------------
for CHR in $CHROMOSOMES; do
    echo ""
    echo "══════════════════════════════════════════"
    echo "  Chromosome $CHR"
    echo "══════════════════════════════════════════"

    echo "[1.${CHR}] Preparing data..."
    bash scripts/prep_data.sh "$CHR"

    echo "[2.${CHR}] Running germline2..."
    bash scripts/run_germline2.sh "$CHR"

    echo "[3.${CHR}] Running Beagle..."
    bash scripts/run_beagle.sh "$CHR"
done

# -- Step 4: analyze and plot results ----------------------------------------------------------------------
echo ""
echo "[4] Analyzing results across chromosomes..."
bash scripts/analyze.sh "$CHROMOSOMES" "$PARENT" "$CHILD"

echo ""
echo "=== Pipeline complete ==="
echo "  Chromosomes: $CHROMOSOMES"
echo "  Parent-child pair: $PARENT vs $CHILD"
echo "  Results: results/"
echo "  Figures: results/figures/"
echo "  Summary: results/summary/"
