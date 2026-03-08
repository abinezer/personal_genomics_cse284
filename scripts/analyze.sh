#!/usr/bin/env bash
set -euo pipefail

# -- Config --------------------------------------------------------------------
PYTHON="${PYTHON:-python3}"
PARENT="NA20317"
CHILD="NA20318"

CHROMOSOMES="${1:-22}"

# -- Setup paths --------------------------------------------------------------------
mkdir -p results/summary results/figures logs

# -- Build argument lists for germline2 and Beagle files --------------------------------------------------------------------
GERM_ARGS=""
BEAGLE_ARGS=""
for CHR in $CHROMOSOMES; do
    GERM_ARGS="$GERM_ARGS results/chr${CHR}/asw_germline2.match"
    BEAGLE_ARGS="$BEAGLE_ARGS results/chr${CHR}/asw_beagle.ibd.gz"
done

# -- Summarize results across chromosomes --------------------------------------------------------------------
echo "Summarizing IBD results for chromosomes: $CHROMOSOMES"
$PYTHON scripts/summarize_ibd_results.py \
    --chromosomes $CHROMOSOMES \
    --germline-matches $GERM_ARGS \
    --beagle-ibds $BEAGLE_ARGS \
    --parent "$PARENT" \
    --child "$CHILD" \
    --outdir results/summary

# -- Generate figures --------------------------------------------------------------------
echo "Generating figures for chromosomes: $CHROMOSOMES"
$PYTHON scripts/plot_analysis.py \
    --chromosomes $CHROMOSOMES \
    --germline-matches $GERM_ARGS \
    --beagle-ibds $BEAGLE_ARGS \
    --parent "$PARENT" \
    --child "$CHILD" \
    --outdir results/figures

echo "Done: results/summary/ and results/figures/"
