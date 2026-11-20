#!/usr/bin/env bash
# PPanGGOLiN: partitioned pangenome graph (persistent / shell / cloud).
set -euo pipefail

source ~/ppanggolin_env/bin/activate

# PPanGGOLiN takes a tab-separated list of genome name and GFF path
for f in gffs/*.gff; do
    printf "%s\t%s\n" "$(basename "$f" .gff)" "$(realpath "$f")"
done > gff_list.tsv

rm -rf ppanggolin_result
ppanggolin workflow --anno gff_list.tsv -o ppanggolin_result -c 8
ppanggolin draw -p ppanggolin_result/pangenome.h5 --tile_plot --ucurve -o ppanggolin_result --force

deactivate
ls ppanggolin_result/
