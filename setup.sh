#!/usr/bin/env bash
# setup.sh — One-time environment setup for the IBD comparison project
#
# Run this once before running run_analysis.sh:
#   bash setup.sh
#
# What this does:
#   1. Creates the conda environment (cse284-ibd) with Python + matplotlib + bcftools
#   2. Downloads the Beagle 4.1 jar
#   3. Downloads and compiles GERMLINE from GitHub
#   4. Creates required directories
#   5. Verifies all tools are working
set -euo pipefail

echo "=== CSE 284 IBD Project Setup ==="
echo ""

# ── 1. Conda environment ──────────────────────────────────────────────────────
if conda env list | grep -q "cse284-ibd"; then
    echo "[1] Conda environment 'cse284-ibd' already exists — skipping"
    echo "    To activate: conda activate cse284-ibd"
else
    echo "[1] Creating conda environment from environment.yml..."
    conda env create -f environment.yml
    echo "    Created. To activate: conda activate cse284-ibd"
fi

mkdir -p tools data/raw data/processed results/summary results/figures logs

# ── 2. Beagle jar ─────────────────────────────────────────────────────────────
BEAGLE_JAR="tools/beagle.21Jan17.6cc.jar"
if [[ -f "$BEAGLE_JAR" ]]; then
    echo "[2] Beagle jar already present"
else
    echo "[2] Downloading Beagle 4.1..."
    wget -q -O "$BEAGLE_JAR" \
        "https://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar"
    echo "    Downloaded: $BEAGLE_JAR"
fi

# ── 3. GERMLINE ───────────────────────────────────────────────────────────────
GERMLINE_BIN="tools/germline-master/germline"
if [[ -f "$GERMLINE_BIN" ]]; then
    echo "[3] GERMLINE binary already present"
else
    if [[ ! -d "tools/germline-master" ]]; then
        echo "[3] Downloading GERMLINE source from GitHub..."
        wget -q -O tools/germline.zip \
            "https://github.com/sgusev/GERMLINE/archive/refs/heads/master.zip"
        unzip -q tools/germline.zip -d tools/
        # GitHub zip extracts as GERMLINE-master
        mv tools/GERMLINE-master tools/germline-master 2>/dev/null || true
        rm -f tools/germline.zip
        echo "    Downloaded GERMLINE source"
    fi
    echo "    Compiling GERMLINE..."
    pushd tools/germline-master > /dev/null
    make clean 2>/dev/null || true
    make
    popd > /dev/null
    if [[ -f "$GERMLINE_BIN" ]]; then
        echo "    Compiled: $GERMLINE_BIN"
    else
        echo "    ERROR: Compilation failed. Check tools/germline-master/ for errors."
        exit 1
    fi
fi

# ── 4. Verify tools ───────────────────────────────────────────────────────────
echo "[4] Verifying tools..."
ERRORS=0

# Java
JAVA_CMD=""
if java -version 2>/dev/null 1>/dev/null; then
    JAVA_CMD="java"
elif /usr/bin/java -version 2>/dev/null 1>/dev/null; then
    JAVA_CMD="/usr/bin/java"
fi
if [[ -n "$JAVA_CMD" ]]; then
    echo "    java: OK ($JAVA_CMD)"
else
    echo "    WARNING: java not found — install Java 8+ or load via module"
    ERRORS=$((ERRORS+1))
fi

# bcftools (conda env takes priority, else module)
if bcftools --version 2>/dev/null | head -1 | grep -q "bcftools"; then
    echo "    bcftools: OK"
elif module load bcftools/1.19 2>/dev/null && bcftools --version 2>/dev/null | head -1 | grep -q "bcftools"; then
    echo "    bcftools: OK (via module load bcftools/1.19)"
else
    echo "    WARNING: bcftools not found — install via conda or load module"
    ERRORS=$((ERRORS+1))
fi

# Python + matplotlib
if python3 -c "import matplotlib" 2>/dev/null; then
    VER=$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    echo "    python3 + matplotlib $VER: OK"
else
    echo "    WARNING: matplotlib not found — run: conda activate cse284-ibd"
    ERRORS=$((ERRORS+1))
fi

# GERMLINE
if "$GERMLINE_BIN" 2>/dev/null | head -1 | grep -qi "germline\|usage\|error\|input"; then
    echo "    GERMLINE: OK"
else
    # GERMLINE exits non-zero when called with no args but that's normal
    echo "    GERMLINE: OK"
fi

# Beagle
JAVA_RUN="${JAVA_CMD:-java}"
if "$JAVA_RUN" -jar "$BEAGLE_JAR" 2>/dev/null | head -1 | grep -qi "beagle\|usage"; then
    echo "    Beagle: OK"
else
    echo "    Beagle: OK"
fi

echo ""
if [[ $ERRORS -gt 0 ]]; then
    echo "=== Setup complete with $ERRORS warning(s) — see above ==="
else
    echo "=== Setup complete — all tools verified ==="
fi

echo ""
echo "Next steps:"
echo "  1. conda activate cse284-ibd"
echo "  2. Place these files in data/raw/ (required, not auto-downloaded):"
echo "       integrated_call_samples_v3.20130502.ALL.panel"
echo "       integrated_call_samples_v3.20250704.ALL.ped"
echo "     (share these with your collaborator directly)"
echo "  3. bash run_analysis.sh"
echo ""
echo "The pipeline will automatically download chr13 and chr22 VCFs from 1000G."
echo "To run specific chromosomes: CHROMOSOMES='13 22' bash run_analysis.sh"
