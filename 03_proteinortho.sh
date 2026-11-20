#!/usr/bin/env bash
# Proteinortho: reciprocal best hits with adaptive similarity, then spectral clustering.
set -euo pipefail

mkdir -p Proteinortho
cd Proteinortho
proteinortho -project=proteinortho_result -p=blastp+ -cpus=10 -singles -selfblast ../CDS_clean/*.txt

# Clusters present in all 10 genomes (column 1 = number of species)
echo "Core clusters:"
awk -F'\t' '$1==10' proteinortho_result.proteinortho.tsv | wc -l
