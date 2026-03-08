#!/usr/bin/env bash
# Cleans up generated files from the pipeline (by run_full.sh)
# DOES NOT delete files generated in setup.sh 
#
# Usage:
#   bash cleanup.sh 

set -euo pipefail

echo "[CLEANUP] Removing generated pipeline outputs..."

rm -rf data/processed/*
rm -rf data/raw/maps/*
rm -rf results/summary results/chr*
rm -rf logs/*
rm -rf data/raw/ALL.chr*.vcf.gz*

echo "[CLEANUP] Done."