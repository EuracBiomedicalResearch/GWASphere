#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# DOUBLE CLUMP + IND-SIG LD EXPANSION + LD-BOUNDARY MERGE
###############################################################################
# Genomic risk loci identified using PLINK following a FUMA-inspired approach
# (Watanabe et al., Nat Commun 2017).
# Independent significant SNPs: P < 5e-8, MAF >= 0.001, r2 < 0.6.
# Lead SNPs identified with secondary clumping (r2 < 0.1).
# Regions within 250 kb merged into genomic risk loci.
###############################################################################

# ---------- CONFIG ----------
OUT_DIR="$1"
INPUT_DIR="$2"
REF_PREFIX="$3"

echo "OUT_DIR: $OUT_DIR"
echo "INPUT_DIR: $INPUT_DIR"
echo "REF_PREFIX: $REF_PREFIX"

WINDOW_KB=10000      # LD expansion window
STEP1_KB=10000       # LD search window for clumping step 1
STEP2_KB=10000       # LD search window for clumping step 2

# Overridable via environment variables (passed from the Shiny app via sbatch --export)
MERGE_KB="${MERGE_KB:-250}"         # merge distance for risk loci (kb)
PTAG="${PTAG:-0.05}"                # GWAS-tagged threshold for candidate SNPs
P1="${P1:-5e-8}"                    # genome-wide significance for indsig
STEP1_R2="${STEP1_R2:-0.6}"         # LD r2 threshold for independent significant SNPs
STEP2_R2="${STEP2_R2:-0.1}"         # LD r2 threshold for lead SNPs
MAF="${MAF:-0.001}"                 # minimum minor allele frequency

mkdir -p "$OUT_DIR"/{tmp,clump,ld,blocks,loci}

TMP="$OUT_DIR/tmp"

# ---------- helpers ----------

# Build per-chr sumstats: ID + PVAL (and optionally BP/CHR if you want)
build_sumstats_chr () {
  local regenie="$1"
  local chr="$2"
  local out="$3"

  zcat "$regenie" | awk -F'\t' -v chr="$chr" -v OFS="\t" '
    NR==1{
      for(i=1;i<=NF;i++) if($i=="LOG10P") col=i
      print "ID","PVAL"
      next
    }
    $1==chr && col {
      p = exp(log(10)*(-$col))
      print $3, p
    }
  ' > "$out"
}

# Build per-chr p-map: ID -> PVAL (for fast filtering P<0.05)
build_pmap_chr () {
  local sumstats="$1"
  local out="$2"
  # skip header
  awk 'NR>1{print $1"\t"$2}' "$sumstats" > "$out"
}

# Extract clumped index SNPs: CHR SNP BP P (header-robust)
extract_index_from_clumped () {
  local clumped="$1"
  local out="$2"
  awk -v OFS="\t" '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="CHR") chr=i
        if($i=="SNP") snp=i
        if($i=="BP")  bp=i
        if($i=="P")   p=i
        if($i=="SP2") sp2=i
      }
      next
    }
    $1=="CHR"{next}
    NF==0{next}
    {
      if(!chr||!snp||!bp||!p){
        print "ERROR: missing CHR/SNP/BP/P in "FILENAME > "/dev/stderr"; exit 1
      }
      print $chr,$snp,$bp,$p
    }
  ' "$clumped" > "$out"
}

# Build mapping indsig -> lead from step2 clumped:
# Each line: lead \t member_indsig
build_lead_to_indsig_map () {
  local step2_clumped="$1"
  local out="$2"

  awk -v OFS="\t" '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP") snp=i
        if($i=="SP2") sp2=i
      }
      next
    }
    $1=="CHR"{next}
    NF==0{next}
    {
      lead=$snp
      # lead itself is an indsig member
      print lead, lead
      # SP2 contains "snp(p),snp(p)" or "NONE"
      if(sp2 && $sp2!="NONE"){
        n=split($sp2,a,",")
        for(i=1;i<=n;i++){
          x=a[i]
          sub(/\(.*/,"",x)
          if(x!="") print lead, x
        }
      }
    }
  ' "$step2_clumped" > "$out"
}

# Run clump
run_clump () {
  local bfile="$1"
  local sumstats="$2"
  local out="$3"
  local p1="$4"
  local p2="$5"
  local r2="$6"
  local kb="$7"

  plink --bfile "$bfile" \
    --maf "$MAF"  \
    --clump "$sumstats" \
    --clump-snp-field ID \
    --clump-field PVAL \
    --clump-p1 "$p1" \
    --clump-p2 "$p2" \
    --clump-r2 "$r2" \
    --clump-kb "$kb" \
    --out "$out" >/dev/null
}

# Run LD expansion for a list (indsig) on a chromosome
run_ld_indsig () {
  local bfile="$1"
  local snplist="$2"
  local out="$3"
  plink --bfile "$bfile" \
    --maf "$MAF"  \
    --r2 \
    --ld-snp-list "$snplist" \
    --ld-window-kb "$WINDOW_KB" \
    --ld-window 999999 \
    --ld-window-r2 "$STEP1_R2" \
    --out "$out" >/dev/null
}

# Filter LD pairs to candidate SNPs (SNP_B) with GWAS P<PTAG, then compute boundaries per indsig
# Output: chr start end indsig
ld_to_indsig_blocks () {
  local ldfile="$1"
  local pmap="$2"         # ID \t P
  local indsig_list="$3"  # lista SNP indsig, uno per riga
  local outbed="$4"

  awk -v OFS="\t" -v PTAG="$PTAG" -v pmap="$pmap" -v ilist="$indsig_list" '
    BEGIN{
      while((getline < pmap)>0){ p[$1]=$2+0 }
      close(pmap)

      while((getline < ilist)>0){ indsig[$1]=1 }
      close(ilist)
    }
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="CHR_A") ca=i
        if($i=="CHR_B") cb=i
        if($i=="SNP_A") sa=i
        if($i=="SNP_B") sb=i
        if($i=="BP_A")  ba=i
        if($i=="BP_B")  bb=i
      }
      if(!ca||!cb||!sa||!sb||!ba||!bb){
        print "ERROR: missing columns in LD header "FILENAME > "/dev/stderr"
        exit 1
      }
      next
    }
    NF==0{next}
    {
      # qui entra il blocco che ti ho mandato:
      # scegli quale SNP della coppia è l independent significant SNP
      if($(sa) in indsig){
        chr=$(ca); a=$(sa); b=$(sb); bpA=$(ba)+0; bpB=$(bb)+0
      } else if($(sb) in indsig){
        chr=$(cb); a=$(sb); b=$(sa); bpA=$(bb)+0; bpB=$(ba)+0
      } else {
        next
      }

      # Always include the indsig SNP position
      if(!(a in minbp) || bpA<minbp[a]) minbp[a]=bpA
      if(!(a in maxbp) || bpA>maxbp[a]) maxbp[a]=bpA
      chrmap[a]=chr

      # Candidate SNPs: proxy b must have GWAS P < 0.05
      if(b in p && p[b] < PTAG){
        if(bpB<minbp[a]) minbp[a]=bpB
        if(bpB>maxbp[a]) maxbp[a]=bpB
      }
    }
    END{
      for(a in minbp){
        s=minbp[a]-1
        if(s<0) s=0
        e=maxbp[a]
        print chrmap[a], s, e, a
      }
    }
  ' "$ldfile" > "$outbed"
}

# Extract all LD pairs (indsig -> proxy) from PLINK .ld
# Output appended: trait chr indsig proxy bpA bpB r2
extract_ld_pairs () {
  local ldfile="$1"
  local trait="$2"
  local chr="$3"
  local indsig_list="$4"
  local pmap="$5"
  local out="$6"

  awk -v OFS="\t" \
      -v trait="$trait" \
      -v chr="$chr" \
      -v ilist="$indsig_list" \
      -v pmap="$pmap" '
    BEGIN{
      while((getline < ilist)>0){ indsig[$1]=1 }
      close(ilist)

      while((getline < pmap)>0){ p[$1]=$2+0 }
      close(pmap)
    }
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP_A") sa=i
        if($i=="SNP_B") sb=i
        if($i=="BP_A")  ba=i
        if($i=="BP_B")  bb=i
        if($i=="R2")    r2=i
      }
      if(!sa || !sb || !ba || !bb || !r2){
        print "ERROR: missing columns in LD file " FILENAME > "/dev/stderr"
        exit 1
      }
      next
    }
    NF==0{next}
    {
      if($(sa) in indsig){
        a=$(sa)
        b=$(sb)
        bpA=$(ba)+0
        bpB=$(bb)+0
      } else if($(sb) in indsig){
        a=$(sb)
        b=$(sa)
        bpA=$(bb)+0
        bpB=$(ba)+0
      } else {
        next
      }

      indsig_p = (a in p ? p[a] : "NA")
      proxy_p  = (b in p ? p[b] : "NA")

      print trait, chr, a, b, bpA, bpB, $(r2)+0, indsig_p, proxy_p
    }
  ' "$ldfile" >> "$out"
}



###############################################################################
# MAIN
###############################################################################

combined="$OUT_DIR/loci/ALL_TRAITS.loci_table.tsv"
echo -e "trait\tlocus_id\tchr\tstart\tend\trepLead\tminP\tleadSNPs\tnLead\tnIndSig" > "$combined"

shopt -s nullglob
regenies=( "$INPUT_DIR"/*.regenie.gz )

[[ ${#regenies[@]} -gt 0 ]] || { 
  echo "No *.regenie.gz found in $INPUT_DIR"
  exit 1
}

for f in "${regenies[@]}"; do
  trait=$(basename "$f" .regenie.gz)
  echo "=============================="
  echo "[TRAIT] $trait"
  # Per-trait accumulators
  lead_p_map="$TMP/${trait}.lead_p.tsv"; : > "$lead_p_map"
  lead_to_indsig_all="$TMP/${trait}.lead_to_indsig.tsv"; : > "$lead_to_indsig_all"
  indsig_blocks_all="$TMP/${trait}.indsig_blocks.bed"; : > "$indsig_blocks_all"

  # LD pairs table: all SNPs in LD with an independent significant SNP (FUMA-like)
  ld_pairs_all="$OUT_DIR/ld/${trait}.indsig_ld_pairs.tsv"
  echo -e "trait\tchr\tindsig\tproxy\tbp_indsig\tbp_proxy\tr2\tindsig_pval\tproxy_pval" > "$ld_pairs_all"

  # ---------- per chr: build sumstats, step1, step2, LD expansion on indsig ----------
  # Auto-detect reference mode if not already set by run_pipeline.sh
  if [[ -z "${REF_MODE:-}" ]]; then
    if [[ -f "${REF_PREFIX}.bed" ]]; then
      REF_MODE="single"
    else
      REF_MODE="split"
    fi
  fi

  for chr in $(seq 1 22); do
    if [[ "$REF_MODE" == "single" ]]; then
      bfile="${REF_PREFIX}"
    else
      bfile="${REF_PREFIX}${chr}"
    fi

    ss="$TMP/${trait}_chr${chr}.sumstats.tsv"
    pmap="$TMP/${trait}_chr${chr}.pmap.tsv"
    build_sumstats_chr "$f" "$chr" "$ss"
    build_pmap_chr "$ss" "$pmap"

    # STEP1 (indsig)
    out1="$OUT_DIR/clump/${trait}_chr${chr}.step1"
    run_clump "$bfile" "$ss" "$out1" "$P1" "$PTAG" "$STEP1_R2" "$STEP1_KB"

    indsig="$TMP/${trait}_chr${chr}.indsig.tsv"  # CHR SNP BP P
    if [[ -s "${out1}.clumped" ]]; then
      extract_index_from_clumped "${out1}.clumped" "$indsig" || true
    else
      : > "$indsig"
    fi
    [[ -s "$indsig" ]] || continue

    # STEP2 (lead) on indsig only
    indsig_ss="$TMP/${trait}_chr${chr}.indsig_sumstats.tsv"
    { echo -e "ID\tPVAL"; awk -v OFS="\t" '{print $2,$4}' "$indsig"; } > "$indsig_ss"

    out2="$OUT_DIR/clump/${trait}_chr${chr}.step2"
    run_clump "$bfile" "$indsig_ss" "$out2" "1" "1" "$STEP2_R2" "$STEP2_KB"

    # Extract leads and P
    leads="$TMP/${trait}_chr${chr}.leads.tsv"    # CHR SNP BP P
    if [[ -s "${out2}.clumped" ]]; then
      extract_index_from_clumped "${out2}.clumped" "$leads" || true
    else
      : > "$leads"
    fi
    [[ -s "$leads" ]] || continue

    awk -v OFS="\t" '{print $2,$4}' "$leads" >> "$lead_p_map"

    # lead -> indsig membership mapping (from step2 clumped)
    l2i="$TMP/${trait}_chr${chr}.lead_to_indsig.tsv"
    build_lead_to_indsig_map "${out2}.clumped" "$l2i"
    cat "$l2i" >> "$lead_to_indsig_all"

    # LD expansion around INDSIG SNPs (r2>=0.6) within WINDOW_KB, then boundaries using P<0.05 filter
    indsig_list="$TMP/${trait}_chr${chr}.indsig_snps.txt"
    awk '{print $2}' "$indsig" | sort -u > "$indsig_list"

    outld="$OUT_DIR/ld/${trait}_chr${chr}.indsigLD"
    run_ld_indsig "$bfile" "$indsig_list" "$outld"

    ldfile="${outld}.ld"
    [[ -s "$ldfile" ]] || continue

    extract_ld_pairs "$ldfile" "$trait" "$chr" "$indsig_list" "$pmap" "$ld_pairs_all"

    blockbed="$TMP/${trait}_chr${chr}.indsig_blocks.bed"
    ld_to_indsig_blocks "$ldfile" "$pmap" "$indsig_list" "$blockbed"
    cat "$blockbed" >> "$indsig_blocks_all"
  done

  # If no leads / no blocks
  if [[ ! -s "$lead_p_map" || ! -s "$indsig_blocks_all" || ! -s "$lead_to_indsig_all" ]]; then
    echo "  [skip] insufficient data (no leads or no indsig blocks) for $trait"
    continue
  fi

  # Dedup lead->P with min P
  lead_p_min="$TMP/${trait}.lead_p.min.tsv"
  awk -v OFS="\t" '{s=$1;p=$2+0; if(!(s in best)||p<best[s]) best[s]=p} END{for(s in best) print s,best[s]}' \
    "$lead_p_map" > "$lead_p_min"
  mv "$lead_p_min" "$lead_p_map"

  # ---------- Merge indsig blocks by lead membership ----------
  # Build quick lookup indsig->(chr,start,end)
  # Then for each lead, union boundaries: start=min, end=max over assigned indsig blocks.
  lead_blocks="$TMP/${trait}.lead_blocks.bed"  # chr start end lead
  awk -v OFS="\t" '
    FNR==NR{
      # indsig blocks bed: chr start end indsig
      chr[$4]=$1; st[$4]=$2+0; en[$4]=$3+0
      next
    }
    {
      lead=$1; ind=$2
      if(ind in chr){
        if(!(lead in min) || st[ind]<min[lead]) min[lead]=st[ind]
        if(!(lead in max) || en[ind]>max[lead]) max[lead]=en[ind]
        leadchr[lead]=chr[ind]
        nind[lead]++
      }
    }
    END{
      for(l in min){
        print leadchr[l], min[l], max[l], l
      }
    }
  ' "$indsig_blocks_all" "$lead_to_indsig_all" > "$lead_blocks"

  if [[ ! -s "$lead_blocks" ]]; then
    echo "  [skip] no lead blocks created for $trait (ID mismatches?)"
    continue
  fi

  # ---------- Merge lead blocks into final risk loci by overlap or within 250kb ----------
  merged="$TMP/${trait}.merged.bed"
  sort -k1,1 -k2,2n "$lead_blocks" \
    | bedtools merge -d $((MERGE_KB*1000)) -c 4 -o distinct \
    > "$merged"

  out_loci="$OUT_DIR/loci/${trait}.loci.tsv"
  echo -e "trait\tlocus_id\tchr\tstart\tend\trepLead\tminP\tleadSNPs\tnLead\tnIndSig" > "$out_loci"

  # Representative lead = minP among leads in locus; also compute nIndSig by summing unique indsig across included leads.
  # For nIndSig we use lead_to_indsig mapping and a per-locus set.
  awk -v OFS="\t" -v trait="$trait" -v pmap="$lead_p_map" -v l2i="$lead_to_indsig_all" '
    BEGIN{
      while((getline < pmap)>0){ p[$1]=$2+0 } close(pmap)
      while((getline < l2i)>0){
        lead=$1; ind=$2
        key=lead SUBSEP ind
        if(!(key in seen)){ seen[key]=1; indlist[lead]=indlist[lead] "," ind }
      }
      close(l2i)
    }
    {
      chr=$1; start=$2; end=$3; leads=$4
      n=split(leads,a,",")
      minp=1e99; rep=""
      # build indsig set for this locus
      delete locInd
      nind=0

      for(i=1;i<=n;i++){
        L=a[i]
        if(L in p && p[L]<minp){ minp=p[L]; rep=L }
        if(rep=="" && !(L in p)) rep=L

        # add indsig from this lead
        s=indlist[L]
        if(s!=""){
          m=split(s,b,",")
          for(j=1;j<=m;j++){
            if(b[j]=="" ) continue
            if(!(b[j] in locInd)){ locInd[b[j]]=1; nind++ }
          }
        }
      }

      locus_id=trait "_locus" NR
      if(minp==1e99) minp="NA"
      print trait,locus_id,chr,start,end,rep,minp,leads,n,nind
    }
  ' "$merged" >> "$out_loci"

  awk 'NR>1' "$out_loci" >> "$combined"

  echo "  [done] loci=$(($(wc -l < "$out_loci")-1)) -> $out_loci"

  tmpdedup="$TMP/${trait}.indsig_ld_pairs.dedup.tsv"
  {
    head -n 1 "$ld_pairs_all"
    tail -n +2 "$ld_pairs_all" \
      | awk -v OFS="\t" '
        {
          key=$1 FS $2 FS $3 FS $4
          if(!(key in best) || $7>best[key]) { best[key]=$7; row[key]=$0 }
        }
        END{ for(k in row) print row[k] }
      ' | sort -k2,2n -k5,5n
  } > "$tmpdedup"
  mv "$tmpdedup" "$ld_pairs_all"
done


#MAP ID to RSID
cd "$INPUT_DIR"

out="$OUT_DIR/loci/all_traits_id_rsid_unique.tsv"
tmp="$(mktemp)"

# header output
printf "ID\tRSID\n" > "$tmp"

for f in *.regenie.gz; do
  trait="$(basename "$f" .regenie.gz)"

  echo "Processing $trait from $f"

  zcat "$f" | awk -v OFS="\t" '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="ID") id=i
        if($i=="RSID") rs=i
      }
      if(!id || !rs){
        print "ERROR: missing ID/RSID columns" > "/dev/stderr"
        exit 1
      }
      next
    }
    { print $id, $rs }
  ' >> "$tmp"
done

{
  head -n 1 "$tmp"
  tail -n +2 "$tmp" | sort -u
} > "$out"

rm -f "$tmp"
echo "Wrote: $out"


LOCI_TSV="$OUT_DIR/loci/ALL_TRAITS.loci_table.tsv"
MAP_TSV="$OUT_DIR/loci/all_traits_id_rsid_unique.tsv"
OUT_TSV="$OUT_DIR/loci/ALL_TRAITS.loci_table_with_rsid.tsv"

awk -v OFS="\t" '
  BEGIN { FS="\t" }

  # 1) carica mapping ID -> RSID (solo primo file)
  NR==FNR {
    if (FNR==1) next
    id=$1; rs=$2
    gsub(/\r/, "", id); gsub(/\r/, "", rs)
    gsub(/^[ \t]+|[ \t]+$/, "", id)
    gsub(/^[ \t]+|[ \t]+$/, "", rs)
    map[id]=rs
    next
  }

  # 2) header del secondo file (LOCI_TSV)
  FNR==1 {
    for (i=1; i<=NF; i++) {
      gsub(/\r/, "", $i)
      if ($i=="leadSNPs") lead_col=i
    }
    if (!lead_col) {
      print "ERROR: leadSNPs column not found in loci table" > "/dev/stderr"
      exit 1
    }
    print $0, "leadRSIDs"
    next
  }

  # 3) per ogni riga: split leadSNPs e mappa
  {
    lead=$lead_col
    gsub(/\r/, "", lead)

    n = split(lead, arr, ",")
    rsids=""
    for (i=1; i<=n; i++) {
      gsub(/\r/, "", arr[i])
      gsub(/^[ \t]+|[ \t]+$/, "", arr[i])
      r = (arr[i] in map ? map[arr[i]] : "NA")
      rsids = rsids (i==1 ? r : "," r)
    }
    print $0, rsids
  }
' "$MAP_TSV" "$LOCI_TSV" > "$OUT_TSV"

echo "Wrote: $OUT_TSV"