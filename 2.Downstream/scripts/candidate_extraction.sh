#!/usr/bin/env bash

# Candidate SNPs were defined as variants in linkage disequilibrium with an independent significant SNP (r2≥0.6) and showing nominal association evidence (P<0.05). For each candidate SNP, the independent significant SNP with the highest LD was retained, while independent significant SNPs themselves were included as candidate variants with r2=1.

OUT_DIR="$1"

set -euo pipefail

PTAG=0.05
R2THR=0.6

mkdir -p "${OUT_DIR}/ld/candidate_snps"

for LDPAIRS in "$OUT_DIR"/ld/*.indsig_ld_pairs.tsv; do

  trait=$(basename "$LDPAIRS" .indsig_ld_pairs.tsv)
  echo "=== Processing trait: $trait ==="

  OUT="$OUT_DIR/ld/candidate_snps/${trait}.candidate_snps.tsv"

  awk -F'\t' -v OFS="\t" -v PTAG="$PTAG" -v R2="$R2THR" '
    NR==1 { next }

 {
    indsig=$3
    cand=$4
    bp_indsig=$5
    bp_proxy=$6
    r2=$7
    indsig_p=$8
    proxy_p=$9

    seen_indsig[indsig]=1
    indsig_bp[indsig]=bp_indsig
    indsig_pval[indsig]=indsig_p

    if(proxy_p != "NA" && r2 >= R2 && proxy_p < PTAG){
      if (!(cand in best) || r2 > best[cand]) {
        best[cand]=r2
        bestindsig[cand]=indsig
        bestbp[cand]=bp_proxy
        bestp[cand]=proxy_p
        type[cand]="proxy"
      }
    }
  }

  END{
    print "candidateSNP","bestIndSigSNP","bp_candidate","maxR2","candidateP","type"

    for(c in best){
      print c,bestindsig[c],bestbp[c],best[c],bestp[c],type[c]
    }

    for(s in seen_indsig){
      if(!(s in best)){
        print s,s,indsig_bp[s],1,indsig_pval[s],"independent_significant"
      }
    }
  }
  ' "$LDPAIRS" > "$OUT"

  echo "  wrote: $OUT"

done