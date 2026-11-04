#!/bin/bash

DIR=$1

cd "$DIR"
mkdir -p mplot

## Filter out variants with -LOG10P <= 2.3
for f in *.regenie.gz; do
    echo "Processing: $f"
    base=$(basename "$f" .regenie.gz)
    out="mplot/${base}_mplot.txt"

    gunzip -c "$f" | awk -F"\t" -v base="$base" '
        NR==1{
            for(i=1;i<=NF;i++) col[$i]=i
            print "trait","chrom","pos","variant","gene_name","freqA1","b","se","p"
        }
        NR>1{
            l10p=$(col["LOG10P"])
            if(l10p~/^[0-9.]+$/ && l10p+0>2.3)
                print base,$1,$2,$3,$(col["GENE_NAME"]),$6,$(col["BETA"]),$(col["SE"]),10^(-l10p)
        }
    ' > "$out"

    # Remove empty file
    nlines=$(wc -l < "$out")
    if [ "$nlines" -le 1 ]; then
        rm -f "$out"
    fi
done