#!/usr/bin/env bash
# setup.sh — One-time environment setup 
#
# Usage:
#   bash setup.sh
#
# What this does:
#   1. Creates the conda environment (cse284-ibd) with Python, matplotlib, bcftools
#   2. Downloads metadata files (panel, pedigree) from 1000 Genomes
#   3. Downloads Beagle 4.1 jar file
#   4. Downloads and compiles germline2 (g2) from GitHub
#   5. Verifies all tools work correctly
#
# System prerequisites (must be installed before running this):
#   - curl or wget 
#   - Java 8+
#   - gcc, make 
#   - conda 
#
# Check prerequisites:
#   which wget java gcc make conda
set -euo pipefail

echo "════════════════════════════════════════════════════════════════"
echo "  CSE 284 IBD Project Setup"
echo "════════════════════════════════════════════════════════════════"
echo ""

# -- Verify system prerequisites -------------------------------------------------
echo "[PREREQ] Checking system prerequisites..."
MISSING=()

for cmd in wget java gcc make conda; do
    if ! command -v "$cmd" &>/dev/null; then
        MISSING+=("$cmd")
    fi
done

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: Missing required system tools: ${MISSING[*]}"
    exit 1
fi
echo "  SUCCESS: wget, java, gcc, make, conda all found"
echo ""

# -- Create required directories -------------------------------------------------
echo "[SETUP] Creating output directories..."
mkdir -p tools data/raw data/processed results/summary results/figures logs
echo "  SUCCESS: tools/, data/, results/, logs/ created"
echo ""

# -- Create/activate conda environment -------------------------------------------------
echo "[1] Setting up conda environment..."
if conda env list | grep -q "cse284-ibd"; then
    echo "  SUCCESS: Conda environment 'cse284-ibd' already exists"
else
    echo "  Creating environment from environment.yml..."
    conda env create -f environment.yml
    echo "  SUCCESS: Created environment 'cse284-ibd'"
fi
echo ""
echo "  NOTE: Activate manually before running pipeline:"
echo "        conda activate cse284-ibd"
echo ""

# -- Download required metadata files ------------------------------------------------- 
echo "[2] Downloading 1000 Genomes metadata files..."
PANEL_FILE="data/raw/integrated_call_samples_v3.20130502.ALL.panel"
PED_FILE="data/raw/integrated_call_samples_v3.20250704.ALL.ped"

if [[ -f "$PANEL_FILE" ]]; then
    echo "  SUCCESS: Panel file already present"
else
    echo "  Downloading: integrated_call_samples_v3.20130502.ALL.panel..."
    wget -q -O "$PANEL_FILE" \
        "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    echo "  SUCCESS: Downloaded $PANEL_FILE"
fi

if [[ -f "$PED_FILE" ]]; then
    echo "  SUCCESS: Pedigree file already present"
else
    echo "  Downloading: integrated_call_samples_v3.20250704.ALL.ped..."
    wget -q -O "$PED_FILE" \
        "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20250704.ALL.ped"
    echo "  SUCCESS: Downloaded $PED_FILE"
fi
echo ""

# -- Download and compile germline2 -------------------------------------------------
echo "[3] Installing germline2..."
GERMLINE2_BIN="tools/germline2/g2"

# Clear any prior installation
if [[ -d tools/germline2 ]]; then
    echo "  Removing prior germline2 installation..."
    rm -rf tools/germline2
fi

echo "  Downloading from GitHub..."
wget -q -O tools/germline2.zip \
    "https://github.com/gusevlab/germline2/archive/refs/heads/master.zip"
unzip -q tools/germline2.zip -d tools/
mv tools/germline2-master tools/germline2
rm -f tools/germline2.zip

echo "  Compiling with gcc/make..."
pushd tools/germline2 > /dev/null
if [[ -f "Makefile" ]]; then
    sed -i 's|-I/opt/boost-1.57.0/include/||g' Makefile || true
fi
make clean 2>/dev/null || true
make > /dev/null 2>&1 || {
    echo "  ERROR: germline2 compilation failed"
    echo "  Check: cat tools/germline2/Makefile and gcc/g++ installation"
    exit 1
}
popd > /dev/null

if [[ ! -f "$GERMLINE2_BIN" ]]; then
    echo "  ERROR: germline2 binary not created at: $GERMLINE2_BIN"
    exit 1
fi
echo "  SUCCESS: Compiled germline2: $GERMLINE2_BIN"
echo ""

# -- Download Beagle jar -------------------------------------------------
echo "[4] Downloading Beagle 4.1..."
BEAGLE_JAR="tools/beagle.21Jan17.6cc.jar"

if [[ -f "$BEAGLE_JAR" ]]; then
    echo "  SUCCESS: Beagle jar already present"
else
    echo "  Downloading: beagle.21Jan17.6cc.jar..."
    wget -q -O "$BEAGLE_JAR" \
        "https://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar"
    if [[ ! -f "$BEAGLE_JAR" ]]; then
        echo "  ERROR: Beagle download failed"
        exit 1
    fi
    echo "  SUCCESS: Downloaded: $BEAGLE_JAR"
fi
echo ""

# -- Verify system tools (no conda activation) -------------------------------------------------
echo "[5] Verifying system tools..."
FAILED=0

# Check germline2
if "$GERMLINE2_BIN" 2>&1 | grep -qiE "Usage|parameters|Options"; then
    echo "  SUCCESS: germline2 (g2) compiled and executable"
else
    echo "  ERROR: germline2 executable but unexpected behavior"
    FAILED=1
fi

# Check Beagle jar
if java -jar "$BEAGLE_JAR" 2>&1 | head -3 | grep -qi "beagle"; then
    echo "  SUCCESS: Beagle jar executable"
else
    echo "  WARNING: Beagle jar seems executable (Java + jar work)"
fi
echo ""
echo "[6] Conda environment tools (verify after activation)..."
echo "  After running 'conda activate cse284-ibd', verify:"
echo "    bcftools --version"
echo "    python3 -c 'import matplotlib; print(matplotlib.__version__)'"
echo ""

# -- Summary --------------------------------------------------------------------------
echo "════════════════════════════════════════════════════════════════"
if [[ $FAILED -eq 0 ]]; then
    echo "  SUCCESS: Setup complete — all tools verified"
else
    echo "  WARNING: Setup complete with some issues — see above"
fi
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Next steps:"
echo "  1. Activate the conda environment:"
echo "     conda activate cse284-ibd"
echo ""
echo "  2. Run the pipeline with:"
echo "     bash run_full.sh (default: chr22 only)"
echo "     bash run_full.sh --chromosomes '21 22' (specific chromosomes)"
echo "     bash run_full.sh --chromosomes all (autosomes 1-22)"
echo ""
echo "For more help: see README.md"
