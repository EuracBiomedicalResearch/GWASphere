#!/bin/bash
#SBATCH --job-name=gwasphere-ollama-setup
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --partition=YOUR_PARTITION

# =============================================================================
# setup_ollama.sh
# Downloads Ollama and pulls the AI model for GWASphere
#
# Run once before using AI reports:
#   sbatch setup_ollama.sh
#
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directory where Ollama .sif and models will be stored
OLLAMA_DIR="${OLLAMA_DIR:-$HOME/gwasphere_tools/ollama}"

# Model to pull — change if you want a different model
# Recommended: qwen3:8b (fast), qwen3:32b (more accurate, needs more GPU memory)
OLLAMA_MODEL="${OLLAMA_MODEL:-qwen3:8b}"

# =============================================================================
# SETUP
# =============================================================================

mkdir -p "$OLLAMA_DIR/models"
cd "$OLLAMA_DIR"

# Step 1 — Build Ollama Singularity image
if [[ -f "$OLLAMA_DIR/ollama.sif" ]]; then
  echo "[INFO] ollama.sif already exists — skipping build"
else
  echo "[INFO] Building ollama.sif from Docker Hub..."
  singularity build ollama.sif docker://ollama/ollama
  echo "[INFO] ollama.sif built successfully"
fi

# Step 2 — Start Ollama server
export OLLAMA_MODELS="$OLLAMA_DIR/models"

singularity exec \
  --bind "$OLLAMA_MODELS:$OLLAMA_MODELS" \
  "$OLLAMA_DIR/ollama.sif" \
  ollama serve > "$OLLAMA_DIR/ollama_serve.log" 2>&1 &

SERVE_PID=$!

# Wait for server to be ready
echo "[INFO] Waiting for Ollama server..."
for i in {1..60}; do
  sleep 2
  if curl -s http://127.0.0.1:11434 >/dev/null 2>&1; then
    echo "[INFO] Ollama server is up"
    break
  fi
  if [[ $i -eq 60 ]]; then
    echo "ERROR: Ollama server did not start" >&2
    kill $SERVE_PID || true
    exit 1
  fi
done

# Step 3 — Pull model
echo "[INFO] Pulling model: $OLLAMA_MODEL"
singularity exec \
  --bind "$OLLAMA_MODELS:$OLLAMA_MODELS" \
  "$OLLAMA_DIR/ollama.sif" \
  ollama pull "$OLLAMA_MODEL"

# Step 4 — Verify
echo "[INFO] Available models:"
singularity exec \
  --bind "$OLLAMA_MODELS:$OLLAMA_MODELS" \
  "$OLLAMA_DIR/ollama.sif" \
  ollama list

# Stop server
kill $SERVE_PID || true
wait $SERVE_PID 2>/dev/null || true

echo ""
echo "=============================================="
echo "  Ollama setup complete!"
echo "  Models dir: $OLLAMA_MODELS"
echo "  SIF:        $OLLAMA_DIR/ollama.sif"
echo ""
echo "  Make sure OLLAMA_DIR in master/run.sh is set to:"
echo "  export OLLAMA_DIR=\"$OLLAMA_DIR\""
echo "=============================================="