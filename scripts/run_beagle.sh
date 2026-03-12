#!/usr/bin/env bash
set -euo pipefail

# -- Config --------------------------------------------------------------------
BEAGLE_JAR="${BEAGLE_JAR:-tools/beagle.21Jan17.6cc.jar}"
PYTHON="${PYTHON:-python3}"
MEMORY="${MEMORY:-4g}"
THREADS="${THREADS:-4}"

CHR="${1:-22}"

# -- Chromosome validation --------------------------------------------------------------------
if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: chromosome must be an integer from 1 to 22 (got: $CHR)"
    exit 1
fi

# -- Prechecks --------------------------------------------------------------------
if ! command -v java &>/dev/null; then
    echo "ERROR: java not found in PATH"
    echo "  Install Java 8+ or ensure it is in PATH"
    exit 1
fi

if [[ ! -f "$BEAGLE_JAR" ]]; then
    echo "ERROR: Beagle jar not found: $BEAGLE_JAR"
    echo "  Run setup.sh to download Beagle"
    exit 1
fi

ASW_VCF="data/processed/chr${CHR}/asw.chr${CHR}.bial.vcf.gz"
if [[ ! -f "$ASW_VCF" ]]; then
    echo "ERROR: ASW VCF not found: $ASW_VCF"
    echo "  Run scripts/prep_data.sh $CHR first"
    exit 1
fi

# -- Setup paths --------------------------------------------------------------------
OUTDIR="results/chr${CHR}"
BEAGLE_OUT="${OUTDIR}/asw_beagle"

mkdir -p "$OUTDIR"

# -- Run Beagle IBD detection ------------------------------------------------------
if [[ ! -f "${BEAGLE_OUT}.ibd.gz" ]]; then
    echo "Running Beagle IBD detection for chr${CHR}..."
    java -Xmx${MEMORY} -jar "$BEAGLE_JAR" \
        gt="$ASW_VCF" \
        out="$BEAGLE_OUT" \
        ibd=true \
        ibdcm=0.3 \
        ibdlod=2.0 \
        nthreads=${THREADS} \
        > "${BEAGLE_OUT}.log" 2>&1
    N_BEAGLE=$($PYTHON -c "import gzip; d=gzip.open('${BEAGLE_OUT}.ibd.gz','rb').read(); print(d.count(b'\n'))")
    echo "Beagle segments: ${N_BEAGLE}"
else
    N_BEAGLE=$($PYTHON -c "import gzip; d=gzip.open('${BEAGLE_OUT}.ibd.gz','rb').read(); print(d.count(b'\n'))")
    echo "Beagle output already present for chr${CHR} (${N_BEAGLE} segments)"
fi

echo "Done: ${BEAGLE_OUT}.ibd.gz"
