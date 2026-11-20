#!/usr/bin/env bash
# Panaroo: graph-based pangenome with error correction (merge, refind, clean).
set -euo pipefail

source ~/panaroo_env/bin/activate
rm -rf panaroo_result
panaroo -i gffs/*.gff -o panaroo_result --clean-mode sensitive --remove-invalid-genes -a core -t 8
deactivate

cat panaroo_result/summary_statistics.txt

# Core-gene tree from the Panaroo core alignment
iqtree2 -s panaroo_result/core_gene_alignment.aln -m GTR+G -B 1000 -T AUTO --prefix panaroo_result/core_tree
