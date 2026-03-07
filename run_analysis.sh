#!/usr/bin/env bash
# run_analysis.sh — IBD comparison pipeline: GERMLINE vs Beagle
#
# Usage:
#   bash run_analysis.sh                        # run all chromosomes in CHROMOSOMES
#   CHROMOSOMES="22" bash run_analysis.sh       # single chromosome
#
# Requirements:
#   java         (Java 8+, used for Beagle)
#   python3      (Python 3.6+, with matplotlib)
#   bcftools     (loaded via module or in PATH)
#   tools/germline-master/germline
#   tools/beagle.21Jan17.6cc.jar
#
# Data in data/raw/:
#   ALL.chr{N}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#   integrated_call_samples_v3.20130502.ALL.panel
#   integrated_call_samples_v3.20250704.ALL.ped
set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
BEAGLE_JAR="${BEAGLE_JAR:-tools/beagle.21Jan17.6cc.jar}"
GERMLINE="${GERMLINE:-tools/germline-master/germline}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
PYTHON="${PYTHON:-python3}"

PANEL="data/raw/integrated_call_samples_v3.20130502.ALL.panel"
PED_FILE="data/raw/integrated_call_samples_v3.20250704.ALL.ped"
POPULATION="ASW"
PARENT="NA20317"
CHILD="NA20318"

# Chromosomes to process (space-separated)
CHROMOSOMES="${CHROMOSOMES:-13 22}"

# FTP base for 1000G Phase 3 VCFs
FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF_PATTERN="ALL.chr{CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

mkdir -p data/processed data/raw results/summary results/figures logs

# ── Load modules ──────────────────────────────────────────────────────────────
module load bcftools/1.19 2>/dev/null || true

# ── Locate java ───────────────────────────────────────────────────────────────
if command -v java &>/dev/null; then
    JAVA="java"
elif [[ -x /usr/bin/java ]]; then
    JAVA="/usr/bin/java"
else
    echo "ERROR: java not found. Install Java 8+ or load it via module." >&2
    exit 1
fi

# ── Compile GERMLINE if binary is missing ─────────────────────────────────────
if [[ ! -x "$GERMLINE" ]]; then
    echo "[setup] GERMLINE binary not found at $GERMLINE — compiling from source..."
    if [[ ! -d "tools/germline-master" ]]; then
        echo "  Downloading GERMLINE source from GitHub..."
        wget -q -O tools/germline.zip \
            "https://github.com/sgusev/GERMLINE/archive/refs/heads/master.zip"
        unzip -q tools/germline.zip -d tools/
        mv tools/GERMLINE-master tools/germline-master 2>/dev/null || true
        rm -f tools/germline.zip
    fi
    pushd tools/germline-master > /dev/null
    make clean 2>/dev/null || true
    make
    popd > /dev/null
    if [[ ! -x "$GERMLINE" ]]; then
        echo "ERROR: GERMLINE compilation failed. Check tools/germline-master/ for errors." >&2
        exit 1
    fi
    echo "  GERMLINE compiled: $GERMLINE"
fi

# ── Download Beagle jar if missing ───────────────────────────────────────────
if [[ ! -f "$BEAGLE_JAR" ]]; then
    echo "[setup] Beagle jar not found — downloading..."
    wget -q -O "$BEAGLE_JAR" \
        "https://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar"
    echo "  Downloaded: $BEAGLE_JAR"
fi

# ── Download required 1000G metadata files if missing ────────────────────────
BASE_1KG="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
if [[ ! -f "$PANEL" ]]; then
    echo "[setup] Downloading population panel..."
    wget -q -O "$PANEL" "${BASE_1KG}/integrated_call_samples_v3.20130502.ALL.panel"
    echo "  Downloaded: $PANEL"
fi
if [[ ! -f "$PED_FILE" ]]; then
    echo "[setup] Downloading pedigree file..."
    wget -q -O "$PED_FILE" "${BASE_1KG}/integrated_call_samples_v3.20250704.ALL.ped"
    echo "  Downloaded: $PED_FILE"
fi

# ═══════════════════════════════════════════════════════════════════
# STEP 0: Extract parent-child pairs from pedigree (informational)
# ═══════════════════════════════════════════════════════════════════
echo "[0] Extracting parent-child pairs from pedigree..."
$PYTHON scripts/extract_parent_child_pairs.py \
    --ped "$PED_FILE" \
    --out data/processed/parent_child_pairs.tsv

# ═══════════════════════════════════════════════════════════════════
# STEP 1: Extract ASW sample list
# ═══════════════════════════════════════════════════════════════════
echo "[1] Building ASW sample list..."
awk -v pop="$POPULATION" 'NR>1 && $2==pop {print $1}' "$PANEL" \
    > data/processed/asw.samples
echo "  $(wc -l < data/processed/asw.samples) ASW samples"

# ═══════════════════════════════════════════════════════════════════
# STEP 2: Per-chromosome processing
# ═══════════════════════════════════════════════════════════════════
for CHR in $CHROMOSOMES; do
    echo ""
    echo "══════════════════════════════════════════"
    echo "  Chromosome $CHR"
    echo "══════════════════════════════════════════"

    RAW_VCF="data/raw/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    OUTDIR="results/chr${CHR}"
    mkdir -p "$OUTDIR" "data/processed/chr${CHR}"

    # ── 2a: Download VCF if not present ─────────────────────────────
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

    # ── 2b: Subset to ASW, biallelic SNPs ───────────────────────────
    if [[ ! -f "$ASW_VCF" ]]; then
        echo "  Subsetting to ASW biallelic SNPs..."
        $BCFTOOLS view \
            -S data/processed/asw.samples \
            -m 2 -M 2 -v snps \
            -O z -o "$ASW_VCF" \
            "$RAW_VCF"
        $BCFTOOLS index -t "$ASW_VCF"
        echo "  Variants: $($BCFTOOLS stats $ASW_VCF | grep '^SN' | grep 'SNPs' | awk '{print $4}')"
    else
        echo "  ASW VCF already present: $ASW_VCF"
    fi

    HAP_PREFIX="data/processed/chr${CHR}/asw_hap"
    GERM_OUT="${OUTDIR}/asw_germline"
    BEAGLE_OUT="${OUTDIR}/asw_beagle"

    # ── 2c: Convert phased VCF to GERMLINE haploid PED ──────────────
    if [[ ! -f "${HAP_PREFIX}.ped" ]]; then
        echo "  Converting phased VCF to GERMLINE haploid format..."
        $PYTHON scripts/vcf_to_germline_hap.py \
            "$ASW_VCF" \
            --out "$HAP_PREFIX" \
            2>&1 | tee "logs/vcf_to_hap_chr${CHR}.log"
    else
        echo "  Haploid PED already present: ${HAP_PREFIX}.ped"
    fi

    # ── 2d: Run GERMLINE ─────────────────────────────────────────────
    if [[ ! -f "${GERM_OUT}.match" ]]; then
        echo "  Running GERMLINE..."
        printf "1\n${HAP_PREFIX}.map\n${HAP_PREFIX}.ped\n${GERM_OUT}\n" | \
            $GERMLINE -min_m 1 > "${GERM_OUT}.log" 2>&1
        echo "  GERMLINE segments: $(wc -l < ${GERM_OUT}.match)"
    else
        echo "  GERMLINE output already present ($(wc -l < ${GERM_OUT}.match) segments)"
    fi

    # ── 2e: Run Beagle IBD detection ────────────────────────────────
    if [[ ! -f "${BEAGLE_OUT}.ibd.gz" ]]; then
        echo "  Running Beagle IBD detection..."
        $JAVA -Xmx4g -jar "$BEAGLE_JAR" \
            gt="$ASW_VCF" \
            out="$BEAGLE_OUT" \
            ibd=true \
            ibdcm=0.5 \
            ibdlod=1e-6 \
            nthreads=4 \
            > "${BEAGLE_OUT}.log" 2>&1
        N_BEAGLE=$($PYTHON -c "import gzip; d=gzip.open('${BEAGLE_OUT}.ibd.gz','rb').read(); print(d.count(b'\n'))")
        echo "  Beagle segments: ${N_BEAGLE}"
    else
        N_BEAGLE=$($PYTHON -c "import gzip; d=gzip.open('${BEAGLE_OUT}.ibd.gz','rb').read(); print(d.count(b'\n'))")
        echo "  Beagle output already present (${N_BEAGLE} segments)"
    fi

done

# ═══════════════════════════════════════════════════════════════════
# STEP 3: Summarize results across chromosomes
# ═══════════════════════════════════════════════════════════════════
echo ""
echo "[3] Summarizing results..."

# Build argument lists for the summarize script
GERM_ARGS=""
BEAGLE_ARGS=""
for CHR in $CHROMOSOMES; do
    GERM_ARGS="$GERM_ARGS results/chr${CHR}/asw_germline.match"
    BEAGLE_ARGS="$BEAGLE_ARGS results/chr${CHR}/asw_beagle.ibd.gz"
done

$PYTHON scripts/summarize_ibd_results.py \
    --chromosomes $CHROMOSOMES \
    --germline-matches $GERM_ARGS \
    --beagle-ibds $BEAGLE_ARGS \
    --parent "$PARENT" \
    --child "$CHILD" \
    --outdir results/summary

# ═══════════════════════════════════════════════════════════════════
# STEP 4: Generate figures
# ═══════════════════════════════════════════════════════════════════
echo ""
echo "[4] Generating figures..."

$PYTHON scripts/plot_analysis.py \
    --chromosomes $CHROMOSOMES \
    --germline-matches $GERM_ARGS \
    --beagle-ibds $BEAGLE_ARGS \
    --parent "$PARENT" \
    --child "$CHILD" \
    --outdir results/figures

echo ""
echo "=== Pipeline complete ==="
echo "  Chromosomes: $CHROMOSOMES"
echo "  Parent-child pair: $PARENT vs $CHILD"
echo "  Results: results/"
echo "  Figures: results/figures/"
echo "  Summary: results/summary/"
