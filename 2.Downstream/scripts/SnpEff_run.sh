#!/usr/bin/env bash
set -euo pipefail

# ---- Config ----
SNPEFF_HOME="${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/snpEff"
SNPEFF_JAR="$SNPEFF_HOME/snpEff.jar"
# Genome version — passed from downstream pipeline via SNPEFF_GENOME env
# Default: GRCh38.99 (hg38). For hg19: GRCh37.75
GENOME="${SNPEFF_GENOME:-GRCh38.99}"

OUT_DIR="$1"

cd ${OUT_DIR}/SnpEff/input/

for f in *_snpeff_input.vcf; do
trait=$(basename "$f" _snpeff_input.vcf)
# ---- Run SnpEff ----
  java -Xmx8g -jar "$SNPEFF_JAR" \
    -v \
    -stats "$OUT_DIR/SnpEff/${trait}.html" \
    "$GENOME" \
    "$f" \
    > "$OUT_DIR/SnpEff/${trait}.ann.vcf"
done
