#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=nextflow

# =============================================================================
# run_pipeline.sh — GWASphere external version
# Launches nf-gwas (REGENIE) via Nextflow on SLURM
#
# Prerequisites:
#   - Java 11+
#   - Nextflow (https://www.nextflow.io/docs/latest/install.html)
#   - Singularity
# =============================================================================

# -----------------------------------------------------------------------------
# Check prerequisites
# -----------------------------------------------------------------------------
if ! command -v java &>/dev/null; then
  echo "ERROR: Java not found. Install from https://www.nextflow.io/docs/latest/install.html" >&2
  exit 1
fi

if ! command -v "$NEXTFLOW_BIN" &>/dev/null; then
  echo "ERROR: Nextflow not found at '${NEXTFLOW_BIN}'." >&2
  echo "       Install from https://www.nextflow.io/docs/latest/install.html" >&2
  exit 1
fi

if ! command -v singularity &>/dev/null; then
  echo "ERROR: Singularity not found. Please install Singularity on your cluster." >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
CONFFILE=$(realpath "$1")
OUTDIR=$(sed -n 's/^params\.outdir *= *"\(.*\)"/\1/p' "$CONFFILE" | head -n1)
OUTDIR=${OUTDIR//$'\r'/}
OTHER_ARGS=("${@:2}")

echo "Running on:        ${SLURM_JOB_NODELIST}"
echo "Config:            ${CONFFILE}"
echo "Output dir:        ${OUTDIR}"
echo "Nextflow:          $(which $NEXTFLOW_BIN)"

# -----------------------------------------------------------------------------
# Java — auto-detect from PATH
# -----------------------------------------------------------------------------
export NXF_JAVA_HOME="${NXF_JAVA_HOME:-$(dirname $(dirname $(readlink -f $(which java))))}"
export PATH="$NXF_JAVA_HOME/bin:$PATH"
export NXF_OPTS="-Djdk.lang.Process.launchMechanism=VFORK"

# Singularity cache for nf-gwas containers
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/nxf_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
echo "Singularity cache: ${NXF_SINGULARITY_CACHEDIR}"

# -----------------------------------------------------------------------------
# Launch
# -----------------------------------------------------------------------------
mkdir -p "$OUTDIR"
export NXF_LOG_FILE="$OUTDIR/.nextflow.log"

cd "$OUTDIR" || exit 1

echo "Start nf-gwas pipeline with: ${CONFFILE}"

"$NEXTFLOW_BIN" run genepi/gwas-regenie -r v1.0.11 -profile slurm \
  -work-dir "$OUTDIR/work" \
  -c "$CONFFILE" \
  "${OTHER_ARGS[@]}"