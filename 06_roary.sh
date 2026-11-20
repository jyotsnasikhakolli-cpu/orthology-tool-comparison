#!/usr/bin/env bash
# Roary: CD-HIT pre-clustering plus gene-neighbourhood splitting.
# Run twice to measure how much the identity threshold alone drives the result.
set -euo pipefail

roary -f roary_result  -e -n -v -p 8        gffs/*.gff   # default 95% identity
roary -f roary_relaxed -e -n -v -p 8 -i 50  gffs/*.gff   # relaxed 50% identity

for d in roary_result roary_relaxed; do
    echo "== $d =="
    cat "$d/summary_statistics.txt"
    echo -n "Clusters in all 10 genomes: "
    awk -F',' 'NR>1 && $4==10' "$d/gene_presence_absence.csv" | wc -l
done
