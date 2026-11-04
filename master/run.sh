#!/bin/bash
#SBATCH --job-name=gwasphere
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodelist=YOUR_COMPUTE_NODE
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/your/path/gwasphere_%j.out
#SBATCH --error=/your/path/gwasphere_%j.err

# =============================================================================
# run.sh — Launch the GWASphere master Shiny app
#
# Usage:
#   sbatch run.sh              # submit to SLURM (recommended on HPC)
#   bash run.sh                # run interactively (local/debug)
#
# To access from your local machine, open an SSH tunnel:
#   ssh -L 3838:COMPUTE_NODE:3838 username@login_node
# then open http://localhost:3838 in your browser.
# =============================================================================

# =============================================================================
# CONFIGURATION
# =============================================================================

# Path to the GWASphere repository — edit this
export GWASPHERE_BASE_DIR="/path/to/gwasphere/"

# Java — edit if 'java' is not in your PATH
# e.g. export NXF_JAVA_HOME="/usr/lib/jvm/java-21-openjdk-amd64"
export NXF_JAVA_HOME="${NXF_JAVA_HOME:-}"
[[ -n "$NXF_JAVA_HOME" ]] && export PATH="$NXF_JAVA_HOME/bin:$PATH"

# Nextflow — edit if 'nextflow' is not in your PATH
# e.g. export NEXTFLOW_BIN="$HOME/bin/nextflow"
export NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"

# Path to the tools folder — edit if you extracted tools to a custom location
export GWASPHERE_TOOLS_DIR="${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}"

# Ollama — edit if you set up Ollama in a custom location
export OLLAMA_DIR="${OLLAMA_DIR:-$HOME/gwasphere_tools/ollama}"
export OLLAMA_MODELS="${OLLAMA_DIR}/models"
# =============================================================================

SCRIPT_DIR="${GWASPHERE_BASE_DIR}/master"

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

PORT=3838 
HOST="0.0.0.0"

echo "HOST value: '$HOST'"
while [[ $# -gt 0 ]]; do
  case $1 in
    --port) PORT="$2"; shift 2 ;;
    --host) HOST="$2"; shift 2 ;;
    --help)
      echo "Usage: sbatch run.sh  or  bash run.sh [--port N] [--host H]"
      exit 0 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

# Activate conda environment
if command -v conda &>/dev/null; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda env list | grep -q "gwasphere"; then
    conda activate gwasphere
    echo "[run.sh] Conda env: gwasphere"
  else
    echo "[run.sh] WARNING: gwasphere env not found. Run: conda env create -f environment.yml"
  fi
fi

echo ""
echo "=============================================="
echo "  GWASphere"
echo "=============================================="
echo "  Base dir:    ${GWASPHERE_BASE_DIR}"
echo "  Tools dir:   ${GWASPHERE_TOOLS_DIR}"
echo "  Nextflow:    ${NEXTFLOW_BIN}"
echo "  Java home:   ${NXF_JAVA_HOME:-system default}"
echo "  URL:         http://localhost:${PORT}"
echo "=============================================="
echo ""

cd "${SCRIPT_DIR}"
Rscript -e "shiny::runApp('${SCRIPT_DIR}', host='0.0.0.0', port=${PORT}, launch.browser=FALSE)"