#!/usr/bin/env bash
#SBATCH --job-name=gwas_ai_traits
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <manifest.tsv>" >&2
  exit 2
fi

MANIFEST_FILE="$(realpath "$1")"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
ROW="$(sed -n "${TASK_ID}p" "$MANIFEST_FILE")"
[[ -n "$ROW" ]] || { echo "No manifest row for array task $TASK_ID" >&2; exit 2; }

IFS=$'\t' read -r EVIDENCE_FILE ASSESSMENT_FILE STATUS_FILE <<<"$ROW"
[[ -n "$EVIDENCE_FILE" && -n "$ASSESSMENT_FILE" && -n "$STATUS_FILE" ]] || {
  echo "Malformed manifest row $TASK_ID" >&2
  exit 2
}

# Use VIEWER_APP_DIR (set by master app) to find the single runner
VIEWER_DIR="${VIEWER_APP_DIR:-$(dirname "$(realpath "$0")")/..}"
SINGLE_RUNNER="${GWAS_AI_SINGLE_RUNNER:-${VIEWER_DIR}/scripts/run_report_ai.sh}"
[[ -f "$SINGLE_RUNNER" ]] || { echo "Single-task runner not found: $SINGLE_RUNNER" >&2; exit 2; }

exec bash "$SINGLE_RUNNER" "$EVIDENCE_FILE" "$ASSESSMENT_FILE" "$STATUS_FILE"
