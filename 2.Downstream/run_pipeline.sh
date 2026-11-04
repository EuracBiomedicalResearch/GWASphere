#!/bin/bash

#SBATCH --mem=128G
#SBATCH --job-name=post_gwas_pip
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH -o pipeline.%j.out
#SBATCH -e pipeline.%j.err

set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

# =============================================================================
# PATHS — overridable via environment variables passed with sbatch --export
# Syntax ${VAR:-default}: uses the env variable if set, otherwise the default.
# =============================================================================
OUT_DIR="${OUT_DIR:-/downstream}"
INPUT_DIR="${INPUT_DIR:-/results}"
mkdir -p "${OUT_DIR}"

# Save REGENIE input directory path so the viewer can find tophits
echo "${INPUT_DIR}" > "${OUT_DIR}/regenie_input_dir.txt"

# Save all run parameters for provenance tracking in the viewer
python3 - << PYEOF
import json, os
from datetime import datetime

params = {
    "run_date":            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "out_dir":             os.environ.get("OUT_DIR", ""),
    "input_dir":           os.environ.get("INPUT_DIR", ""),
    "ref_prefix":          os.environ.get("REF_PREFIX", ""),
    "clumping": {
        "P1":              os.environ.get("P1", "5e-8"),
        "PTAG":            os.environ.get("PTAG", "0.05"),
        "STEP1_R2":        os.environ.get("STEP1_R2", "0.6"),
        "STEP2_R2":        os.environ.get("STEP2_R2", "0.1"),
        "MERGE_KB":        os.environ.get("MERGE_KB", "250"),
        "MAF":             os.environ.get("MAF", "0.001"),
    },
    "novelty": {
        "trait_map_xlsx":         os.environ.get("TRAIT_MAP_XLSX", ""),
        "catalog_tsv":            os.environ.get("CATALOG_TSV", ""),
        "default_trait_string":   os.environ.get("DEFAULT_TRAIT_STRING", ""),
        "window_kb":              "500",
        "sig_level":              "5e-8",
    },
    "magma": {
        "gene_loc":        os.environ.get("GENE_LOC", ""),
        "geneset":         os.environ.get("GENESET", ""),
        "gtex":            os.environ.get("GTEX", ""),
        "parallel":        os.environ.get("MAGMA_PARALLEL", "8"),
    },
    "steps_run": {
        "magma":           os.environ.get("RUN_MAGMA", "true"),
        "novelty_pubmed":  os.environ.get("RUN_NOVELTY", "true"),
    },
    "ai_reports_dir": os.path.join(os.environ.get("OUT_DIR", "."), "ai_reports")
}

out_file = os.path.join(os.environ.get("OUT_DIR", "."), "run_params.json")
with open(out_file, "w") as f:
    json.dump(params, f, indent=2)
print(f"Saved run parameters to {out_file}")
PYEOF

REF_PREFIX="${REF_PREFIX:-}"
if [[ -z "$REF_PREFIX" ]]; then
  echo "ERROR: REF_PREFIX is not set. Please provide the path to the PLINK reference files via --export REF_PREFIX=/path/to/imputed_chr"
  exit 1
fi
# Auto-detect: single merged file or per-chromosome files
if [[ -f "${REF_PREFIX}.bed" ]]; then
  REF_MODE="single"
  echo "PLINK reference: single merged file detected (${REF_PREFIX}.bed)"
else
  REF_MODE="split"
  echo "PLINK reference: per-chromosome files detected (${REF_PREFIX}1.bed, ...)"
fi
export REF_MODE

# GWAS Catalog
TRAIT_MAP_XLSX="${TRAIT_MAP_XLSX:-path/to/file/}"  # Excel file mapping trait to EFO ID (columns: trait, EFO_ID)
CATALOG_TSV="${CATALOG_TSV:-${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/GWAS_catalog/gwas-catalog-download-associations-alt-full.tsv}"
# Restore spaces encoded as underscores (sbatch --export cannot handle spaces in values)
DEFAULT_TRAIT_STRING="${DEFAULT_TRAIT_STRING:-}"
DEFAULT_TRAIT_STRING="${DEFAULT_TRAIT_STRING//_/ }"

# MAGMA reference files
GENE_LOC="${GENE_LOC:-${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/MAGMA/NCBI38/NCBI38.gene.loc}"
# REF_FILE_CHR: with #CHR# for split files, without for single file
if [[ "$REF_MODE" == "single" ]]; then
  REF_FILE_CHR="${REF_FILE_CHR:-${REF_PREFIX}}"
else
  REF_FILE_CHR="${REF_FILE_CHR:-${REF_PREFIX}#CHR#}"
fi
export REF_FILE_CHR
GENESET="${GENESET:-${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/MAGMA/MSigDB/c5.all.v2025.1.Hs.entrez.gmt}"
GTEX="${GTEX:-${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/MAGMA/GTEx/gtex_v8_ts_general_avg_log2TPM_entrez.txt}"
MAGMA_PARALLEL="${MAGMA_PARALLEL:-8}"

# Optional steps — set to "false" to skip
RUN_MAGMA="${RUN_MAGMA:-true}"
RUN_NOVELTY="${RUN_NOVELTY:-true}"

# SLURM partition — passed from the Shiny app via --export
SLURM_PARTITION="${SLURM_PARTITION:-computational}"

# =============================================================================
# STEP 1: MAGMA gene-based analysis (submitted as SLURM array)
# =============================================================================
if [[ "$RUN_MAGMA" == "true" ]]; then

echo "STEP 8: MAGMA gene-based analysis (submitted as SLURM array)"

ls "${INPUT_DIR}"/*.regenie.gz > "${OUT_DIR}/traits.txt"
mkdir -p "${OUT_DIR}/logs"
N_TRAITS=$(wc -l < "${OUT_DIR}/traits.txt")

if [[ "$N_TRAITS" -eq 0 ]]; then
  echo "ERROR: no .regenie.gz files found in ${INPUT_DIR}"
  exit 1
fi

MAGMA_JOB_ID=$(sbatch --parsable \
  --array=1-${N_TRAITS}%${MAGMA_PARALLEL} \
  --partition="${SLURM_PARTITION}" \
  --cpus-per-task="${SLURM_CPUS_PER_TASK:-4}" \
  --mem=64G \
  --output="${OUT_DIR}/logs/magma_%A_%a.out" \
  --error="${OUT_DIR}/logs/magma_%A_%a.err" \
  scripts/magma_one_trait.sh \
    "${OUT_DIR}/traits.txt" \
    "$OUT_DIR" \
    "$INPUT_DIR" \
    "$GENE_LOC" \
    "$REF_FILE_CHR" \
    "$GENESET" \
    "$GTEX")

echo "Submitted MAGMA array job: $MAGMA_JOB_ID"

# =============================================================================
# STEP 2: MAGMA gene symbol conversion (depends on STEP 1)
# =============================================================================
echo "STEP 8b: MAGMA gene symbol conversion"

MAGMA_CONVERT_JOB_ID=$(sbatch --parsable \
  --dependency=afterok:${MAGMA_JOB_ID} \
  --partition="${SLURM_PARTITION}" \
  --mem=5G \
  --output="${OUT_DIR}/logs/magma_convert_symbols_%j.out" \
  --error="${OUT_DIR}/logs/magma_convert_symbols_%j.err" \
  --wrap="Rscript scripts/magma_convert_symbols.R \"$OUT_DIR\"")

echo "Submitted MAGMA conversion job: $MAGMA_CONVERT_JOB_ID"

else
  echo "STEP 8: MAGMA skipped (RUN_MAGMA=false)"
fi

# =============================================================================
# STEP 3 & 4: get_mdata and QQ-plots (run in parallel)
# =============================================================================
echo "STEP 1: GWAS summary statistics (get_mdata)"
bash scripts/get_mdata.sh "$INPUT_DIR" &
PID_GET_MDATA=$!

echo "STEP 1b: QQ-plots"
Rscript scripts/qqplots_faster.R "$OUT_DIR" "$INPUT_DIR" &
PID_QQPLOTS=$!

wait "$PID_GET_MDATA" || {
  echo "ERROR: get_mdata failed"
  exit 1
}
echo "STEP 1 completed: GWAS summary statistics"

# =============================================================================
# STEP 3b: Manhattan plots (runs after get_mdata)
# =============================================================================
echo "STEP 1c: Manhattan plots"
Rscript scripts/MH_plots.R "$OUT_DIR" "$INPUT_DIR" || {
  echo "ERROR: Manhattan plots failed"
  exit 1
}
echo "STEP 1c completed: Manhattan plots"

# =============================================================================
# STEPS 6-13: main analysis branch (clumping → locus zoom), runs in background
# =============================================================================
echo "STEP 2-7: starting main analysis branch"
(
  set -e

  echo "STEP 2: LD clumping and locus definition"
  bash scripts/clump_merge_locus.sh \
    "$OUT_DIR" \
    "$INPUT_DIR" \
    "$REF_PREFIX"

  echo "STEP 3: candidate SNP extraction"
  bash scripts/candidate_extraction.sh "$OUT_DIR"

  if [[ "$RUN_NOVELTY" == "true" ]]; then
    echo "STEP 4: novelty detection (GWAS Catalog)"
    python scripts/novel_detection.py \
      "$OUT_DIR" \
      "$INPUT_DIR" \
      "$TRAIT_MAP_XLSX" \
      "$CATALOG_TSV" \
      "$DEFAULT_TRAIT_STRING"

    echo "STEP 5: PubMed literature search"
    Rscript scripts/PubMed.R \
      "$OUT_DIR" \
      "$TRAIT_MAP_XLSX" \
      "$DEFAULT_TRAIT_STRING" &
    PID_PUBMED=$!
  else
    echo "STEP 4-5: Novelty/GWAS Catalog + PubMed skipped (RUN_NOVELTY=false)"
    PID_PUBMED=""
  fi

  echo "STEP 6: SnpEff annotation"
  Rscript scripts/SnpEff_input_prep.R \
    "$INPUT_DIR" \
    "$OUT_DIR"

  echo "STEP 6b: SnpEff running"
  bash scripts/SnpEff_run.sh "$OUT_DIR"

  echo "STEP 6c: SnpEff output refinement"
  Rscript scripts/SnpEff_refine_output.R "$OUT_DIR"

  if [[ -n "${PID_PUBMED:-}" ]]; then
    wait "$PID_PUBMED" || {
      echo "ERROR: PubMed search failed"
      exit 1
    }
  fi

  echo "STEP 7: LocusZoom regional plots"
  Rscript scripts/locuszoom.R \
    "$OUT_DIR" \
    "$INPUT_DIR" \
    "$GENE_LOC"

  echo "STEP 2-7 completed"
) &
PID_MAIN_BRANCH=$!

# Wait for both parallel branches
FAILED=0

wait "$PID_QQPLOTS" || {
  echo "ERROR: QQ-plots failed"
  FAILED=1
}

wait "$PID_MAIN_BRANCH" || {
  echo "ERROR: main analysis branch failed"
  FAILED=1
}

if [[ "$FAILED" -ne 0 ]]; then
  echo "ERROR: one or more pipeline branches failed"
  exit 1
fi

# Gene-based Manhattan plots (after MAGMA conversion completes)
if [[ "$RUN_MAGMA" == "true" ]]; then
  echo "STEP 8c: Gene-based Manhattan plots"
  Rscript scripts/magma_gene_manhattan.R "$OUT_DIR" || {
    echo "WARNING: Gene-based Manhattan plots failed (non-fatal)"
  }
fi

echo "All pipeline steps completed."
echo "AI reports can be generated on-demand from the Explore Results tab."

echo "Pipeline completed — local steps done, submitted jobs running in background"