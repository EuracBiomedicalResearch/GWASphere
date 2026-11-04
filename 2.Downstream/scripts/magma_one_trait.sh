#!/usr/bin/env bash
set -euo pipefail

export PATH="${GWASPHERE_TOOLS_DIR:-$HOME/gwasphere_tools}/MAGMA:$PATH"

# =========================================================
# INPUT
# =========================================================
TRAIT_LIST="$1"

f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TRAIT_LIST")

if [[ -z "${f:-}" ]]; then
  echo "ERROR: no file for array index ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi


# =========================================================
# CONFIG
# =========================================================
OUT_DIR="${2}/magma"
INPUT_DIR="$3"

GENE_LOC="$4"
REF_FILE_CHR="$5"
GENESET="$6"
GTEX="$7"

# LOG10P_COL auto-detected from header
# N_COL auto-detected from header

MAX_WORKERS="${SLURM_CPUS_PER_TASK:-4}"
# MAGMA batch mode requires at least 2 batches
[[ "$MAX_WORKERS" -lt 2 ]] && MAX_WORKERS=2

# REF_MODE inherited from run_pipeline.sh via --export ALL
# Auto-detect as fallback if not set
REF_MODE="${REF_MODE:-}"
if [[ -z "$REF_MODE" ]]; then
  [[ -f "${REF_FILE_CHR}.bed" ]] && REF_MODE="single" || REF_MODE="split"
  echo "[INFO] REF_MODE auto-detected: ${REF_MODE}"
fi

mkdir -p \
  "$OUT_DIR/annotation" \
  "$OUT_DIR/gene_analysis" \
  "$OUT_DIR/gene_set_analysis" \
  "$OUT_DIR/tissue_analysis" \
  "$OUT_DIR/tmp" \
  logs

cd "$INPUT_DIR"

trait="$(basename "$f" .regenie.gz)"

echo "========================================"
echo "[TRAIT] $trait"
echo "[FILE]  $f"

# =========================================================
# FILES
# =========================================================
TMP_SNPS_LOC="$OUT_DIR/tmp/${trait}.snploc.tsv"
SNP_PVAL_FILE="$OUT_DIR/tmp/${trait}.pval.tsv"

ANNO_PREFIX="$OUT_DIR/annotation/${trait}"
ANNO_LOC="${ANNO_PREFIX}.genes.annot"

BATCH_PREFIX="$OUT_DIR/gene_analysis/${trait}"
GENE_PREFIX="$OUT_DIR/gene_analysis/${trait}"
GENE_RAW="${GENE_PREFIX}.genes.raw"

# =========================================================
# 1) Create SNP location + p-value files
# =========================================================
echo "[INFO] Creating MAGMA input files"

# Auto-detect PLINK .bim ID format (chr-prefixed or numeric)
if [[ "$REF_MODE" == "single" ]]; then
  _bim_check="${REF_FILE_CHR}.bim"
else
  _bim_check="${REF_FILE_CHR//#CHR#/1}.bim"
fi
BIM_HAS_CHR=0
if [[ -f "$_bim_check" ]]; then
  _sample_id=$(awk 'NR==2{print $2; exit}' "$_bim_check")
  [[ "$_sample_id" == chr* ]] && BIM_HAS_CHR=1
  echo "[INFO] PLINK .bim ID format: ${_sample_id} (bim_has_chr=${BIM_HAS_CHR})"
fi

zcat "$f" | awk -F"\t" \
  -v snploc="$TMP_SNPS_LOC" \
  -v pval="$SNP_PVAL_FILE" \
  -v bim_has_chr="$BIM_HAS_CHR" '
  BEGIN{
    OFS="\t"
    print "SNP","CHR","BP" > snploc
    print "ID","PVAL","N" > pval
  }
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="LOG10P") pcol=i
      if($i=="N") ncol=i
    }
    next
  }
  pcol && NF>=pcol && $pcol!="" && $pcol!="NA" {

    p = exp(log(10)*(-$pcol))

    # Match ID format to PLINK .bim:
    # - If .bim uses chr-prefixed IDs (e.g. chr1:pos:A:T) → keep chr prefix
    # - If .bim uses numeric IDs (e.g. 1:pos:A:T) → strip chr prefix
    chrom = $1
    sub(/^chr/, "", chrom)
    if (bim_has_chr) {
      id = "chr" chrom ":" $2 ":" $4 ":" $5
    } else {
      id = chrom ":" $2 ":" $4 ":" $5
    }

    print id,chrom,$2 > snploc
    print id,p,(ncol?$ncol:"NA") > pval
  }
'

# =========================================================
# 2) Annotation
# =========================================================
if [[ -s "$ANNO_LOC" ]]; then
  echo "[INFO] Annotation exists: $ANNO_LOC"
else
  echo "[INFO] Running MAGMA annotation"
  magma \
    --annotate window=35,10 \
    --snp-loc "$TMP_SNPS_LOC" \
    --gene-loc "$GENE_LOC" \
    --out "$ANNO_PREFIX"
fi

# =========================================================
# 3) Gene analysis
# =========================================================
if [[ -s "$GENE_RAW" ]]; then
  echo "[INFO] Gene results already exist: $GENE_RAW"
else
  echo "[INFO] Running MAGMA gene analysis"

  rm -f "${BATCH_PREFIX}".batch* "${BATCH_PREFIX}".log

  mkdir -p "${OUT_DIR}/gene_analysis/tmp_batches"
  BATCH_TMP="${OUT_DIR}/gene_analysis/tmp_batches/${trait}"
  PARALLEL_LOG="${OUT_DIR}/gene_analysis/tmp_batches/${trait}_parallel.log"

  # Export variables needed inside the parallel functions
  export BATCH_TMP ANNO_LOC SNP_PVAL_FILE MAX_WORKERS REF_FILE_CHR

  if [[ -f "${REF_FILE_CHR}.bed" ]]; then
    # ---- Mode A: single merged file — numeric batch mode ----
    echo "[INFO] MAGMA reference: single merged file — numeric batch mode"

    seq 1 "$MAX_WORKERS" | parallel \
      --jobs "$MAX_WORKERS" \
      --joblog "$PARALLEL_LOG" \
      --halt now,fail=1 \
      magma \
        --batch {} "$MAX_WORKERS" \
        --bfile "$REF_FILE_CHR" synonyms=0 \
        --gene-annot "$ANNO_LOC" \
        --gene-model snp-wise=mean \
        --pval "$SNP_PVAL_FILE" use=ID,PVAL ncol=N \
        --out "$BATCH_TMP"

  else
    # ---- Mode B: per-chromosome files — chromosome batch mode ----
    declare -a CHR_LABELS=()
    shopt -s nullglob
    for bim in "${REF_FILE_CHR//#CHR#/}"*.bim; do
      prefix="${bim%.bim}"
      suffix="${prefix#${REF_FILE_CHR//#CHR#/}}"
      if [[ "$suffix" =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]; then
        [[ -f "${prefix}.bed" && -f "${prefix}.fam" ]] && CHR_LABELS+=( "$suffix" )
      fi
    done
    shopt -u nullglob

    if [[ ${#CHR_LABELS[@]} -eq 0 ]]; then
      echo "ERROR: no PLINK files found for prefix: ${REF_FILE_CHR//#CHR#/}" >&2
      exit 1
    fi

    echo "[INFO] MAGMA reference: per-chromosome — chromosome batch mode"
    echo "[INFO] Detected chromosomes: ${CHR_LABELS[*]}"

    printf '%s\n' "${CHR_LABELS[@]}" | parallel \
      --jobs "$MAX_WORKERS" \
      --joblog "$PARALLEL_LOG" \
      --halt now,fail=1 \
      magma \
        --batch {} chr \
        --bfile "$REF_FILE_CHR" synonyms=0 \
        --gene-annot "$ANNO_LOC" \
        --gene-model snp-wise=mean \
        --pval "$SNP_PVAL_FILE" use=ID,PVAL ncol=N \
        --out "$BATCH_TMP"
  fi

  echo "[INFO] Parallel job log:"
  cat "$PARALLEL_LOG"

  echo "[INFO] Merging MAGMA batches"
  magma \
    --merge "$BATCH_TMP" \
    --out "$GENE_PREFIX"

  rm -f "${BATCH_TMP}".batch* "${BATCH_TMP}"*.genes.out "${BATCH_TMP}"*.genes.raw
fi

# =========================================================
# 4) Check merged output
# =========================================================
if [[ ! -s "$GENE_RAW" ]]; then
  echo "ERROR: missing merged genes.raw: $GENE_RAW"
  exit 1
fi

# =========================================================
# 5) Gene-set analysis
# =========================================================
if [[ -s "$GENESET" ]]; then
  echo "[INFO] Running gene-set analysis"
  magma \
    --gene-results "$GENE_RAW" \
    --set-annot "$GENESET" \
    --out "$OUT_DIR/gene_set_analysis/${trait}"
fi

# =========================================================
# 6) Tissue analysis
# =========================================================
if [[ -s "$GTEX" ]]; then
  echo "[INFO] Running tissue analysis"
  magma \
    --gene-results "$GENE_RAW" \
    --gene-covar "$GTEX" \
    --model direction=greater \
    --out "$OUT_DIR/tissue_analysis/${trait}"
fi

# =========================================================
# CLEANUP
# =========================================================
rm -f "$TMP_SNPS_LOC" "$SNP_PVAL_FILE"
rm -f "${BATCH_PREFIX}".batch* "${BATCH_PREFIX}".log "${PARALLEL_LOG}"

echo "[DONE] $trait"