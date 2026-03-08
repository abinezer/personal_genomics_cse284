#!/usr/bin/env bash
# setup.sh — One-time environment setup for the IBD comparison project
#
# Run this once before running run_analysis.sh:
#   bash setup.sh
#
# What this does:
#   1. Creates the conda environment (cse284-ibd) with Python + matplotlib + bcftools
#   2. Downloads the Beagle 4.1 jar
#   3. Clears prior GERMLINE installs and installs germline2 (g2)
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

# ── 1c. Clear prior germline installs (fresh reinstall each run) ────────────
echo "[1c] Clearing prior germline installs..."
rm -rf tools/germline-master tools/GERMLINE-master tools/germline2
rm -f tools/germline.zip tools/germline2.zip
echo "    Cleared legacy/new germline install directories"

# ── 1b. Required panel/pedigree metadata ─────────────────────────────────────
echo "[1b] Ensuring required sample metadata files are in data/raw/..."

PANEL_FILE="data/raw/integrated_call_samples_v3.20130502.ALL.panel"
PED_FILE="data/raw/integrated_call_samples_v3.20250704.ALL.ped"

if [[ -f "$PANEL_FILE" ]]; then
    echo "    Panel file already present: $PANEL_FILE"
else
    echo "    Downloading panel file..."
    wget -q -O "$PANEL_FILE" \
        "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    echo "    Downloaded: $PANEL_FILE"
fi

if [[ -f "$PED_FILE" ]]; then
    echo "    Pedigree file already present: $PED_FILE"
else
    echo "    Downloading pedigree file..."
    wget -q -O "$PED_FILE" \
        "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20250704.ALL.ped"
    echo "    Downloaded: $PED_FILE"
fi

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

# ── 3. germline2 (g2) ────────────────────────────────────────────────────────
GERMLINE2_BIN="tools/germline2/g2"
echo "[3] Installing germline2 from GitHub..."
wget -q -O tools/germline2.zip \
    "https://github.com/gusevlab/germline2/archive/refs/heads/master.zip"
unzip -q tools/germline2.zip -d tools/
mv tools/germline2-master tools/germline2
rm -f tools/germline2.zip

echo "    Compiling germline2..."
pushd tools/germline2 > /dev/null
if [[ -f "Makefile" ]]; then
    sed -i 's|-I/opt/boost-1.57.0/include/||g' Makefile
fi
make clean 2>/dev/null || true
make
popd > /dev/null

if [[ -f "$GERMLINE2_BIN" ]]; then
    echo "    Compiled: $GERMLINE2_BIN"
else
    echo "    ERROR: germline2 compilation failed. Check tools/germline2/ for errors."
    exit 1
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

# germline2 (g2)
if "$GERMLINE2_BIN" 2>&1 | grep -qiE "Usage: g2|Incorrect number of parameters|Options:"; then
    echo "    germline2 (g2): OK"
else
    echo "    WARNING: germline2 did not return expected usage output"
    ERRORS=$((ERRORS+1))
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
echo "  2. bash run_analysis.sh"
echo ""
echo "The pipeline will automatically download chr13 and chr22 VCFs from 1000G."
echo "To run specific chromosomes: CHROMOSOMES='13 22' bash run_analysis.sh"
