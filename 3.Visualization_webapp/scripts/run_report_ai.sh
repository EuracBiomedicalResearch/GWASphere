#!/usr/bin/env bash
#SBATCH --job-name=gwas_ai_insights
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=ai_insights.%j.out
#SBATCH --error=ai_insights.%j.err

set -Eeuo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <evidence.json> <assessment.json> <status.json>" >&2
  exit 2
fi

EVIDENCE_FILE="$(realpath "$1")"
ASSESSMENT_FILE="$2"
STATUS_FILE="$3"
STATUS_DIR="$(dirname "$STATUS_FILE")"
mkdir -p "$STATUS_DIR"

export OLLAMA_CONTEXT_LENGTH="${OLLAMA_CONTEXT_LENGTH:-16384}"
export OLLAMA_NUM_THREADS="${OLLAMA_NUM_THREADS:-8}"
export OLLAMA_NUM_PARALLEL="${OLLAMA_NUM_PARALLEL:-1}"
export OLLAMA_MODEL="${OLLAMA_MODEL:-qwen3:8b}"
export OLLAMA_NUM_CTX="${OLLAMA_NUM_CTX:-16384}"
export OLLAMA_NUM_PREDICT="${OLLAMA_NUM_PREDICT:-2600}"

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  OLLAMA_PORT="${OLLAMA_PORT:-$((20000 + SLURM_JOB_ID % 20000))}"
else
  OLLAMA_PORT="${OLLAMA_PORT:-11434}"
fi
export OLLAMA_HOST="127.0.0.1:${OLLAMA_PORT}"
export OLLAMA_URL="http://${OLLAMA_HOST}/api/generate"
# The image defines OLLAMA_HOST internally. Singularity gives precedence to
# SINGULARITYENV_* variables, so bind Ollama inside the container to the same
# per-job port that the host-side readiness check and Python client use.
export SINGULARITYENV_OLLAMA_HOST="0.0.0.0:${OLLAMA_PORT}"

OLLAMA_SIF="${OLLAMA_SIF:-${OLLAMA_DIR:-$HOME/gwasphere_tools/ollama}/ollama.sif}"
export OLLAMA_MODELS="${OLLAMA_MODELS:-${OLLAMA_DIR:-$HOME/gwasphere_tools/ollama}/models}"
OLLAMA_PID=""
VIEWER_DIR="${VIEWER_APP_DIR:-$(dirname "$(realpath "$0")")/..}"
WORKER="${GWAS_AI_WORKER:-${VIEWER_DIR}/scripts/generate_report_insights.py}"

write_shell_status() {
  local state="$1"
  local error_text="${2:-}"
  python - "$STATUS_FILE" "$state" "$error_text" "${SLURM_JOB_ID:-}" <<'PY'
import json, os, sys, tempfile
path, state, error, job_id = sys.argv[1:]
payload = {"status": state, "job_id": job_id or None, "error": error or None}
directory = os.path.dirname(path)
fd, tmp = tempfile.mkstemp(prefix="status_", suffix=".json", dir=directory)
with os.fdopen(fd, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2)
os.replace(tmp, path)
PY
}

cleanup() {
  rc=$?
  if [[ -n "$OLLAMA_PID" ]] && kill -0 "$OLLAMA_PID" 2>/dev/null; then
    kill "$OLLAMA_PID" 2>/dev/null || true
    wait "$OLLAMA_PID" 2>/dev/null || true
  fi
  if [[ $rc -ne 0 ]]; then
    write_shell_status "failed" "GPU job exited with code $rc" || true
  fi
  exit "$rc"
}
trap cleanup EXIT INT TERM

[[ -f "$EVIDENCE_FILE" ]] || { write_shell_status failed "Evidence file not found"; exit 2; }
[[ -f "$WORKER" ]] || { write_shell_status failed "AI worker not found: $WORKER"; exit 2; }
[[ -f "$OLLAMA_SIF" ]] || { write_shell_status failed "Ollama container not found"; exit 2; }
[[ -d "$OLLAMA_MODELS" ]] || { write_shell_status failed "Ollama model directory not found"; exit 2; }


write_shell_status "starting"
singularity exec --nv --bind "$OLLAMA_MODELS:$OLLAMA_MODELS" \
  "$OLLAMA_SIF" ollama serve >"$STATUS_DIR/ollama_${SLURM_JOB_ID:-local}.log" 2>&1 &
OLLAMA_PID=$!

ready=0
for _ in $(seq 1 90); do
  if ! kill -0 "$OLLAMA_PID" 2>/dev/null; then
    write_shell_status failed "Ollama exited before readiness"
    exit 1
  fi
  if curl --fail --silent --max-time 2 "http://${OLLAMA_HOST}/api/tags" >/dev/null; then
    ready=1
    break
  fi
  sleep 2
done
[[ "$ready" -eq 1 ]] || { write_shell_status failed "Ollama readiness timeout"; exit 1; }

python "$WORKER" "$EVIDENCE_FILE" "$ASSESSMENT_FILE" "$STATUS_FILE"