#!/usr/bin/env bash
# Fetch GenBank records, protein FASTA and whole-genome FASTA for each accession.
# Requires NCBI Entrez Direct (sudo apt install ncbi-entrez-direct).
set -euo pipefail

ACCESSIONS=${1:-accessions.txt}
mkdir -p gbks CDS genomes

while read -r acc; do
    [[ -z "$acc" || "$acc" == \#* ]] && continue
    echo "Fetching $acc"
    efetch -db nucleotide -id "$acc" -format gb          > "gbks/$acc.gbk"
    efetch -db nucleotide -id "$acc" -format fasta_cds_aa > "CDS/$acc.txt"
    efetch -db nucleotide -id "$acc" -format fasta        > "genomes/$acc.fasta"
    sleep 1
done < "$ACCESSIONS"

# Sanity check: protein count and CDS feature count should match per genome
echo -e "\naccession\tproteins\tCDS_features"
for f in gbks/*.gbk; do
    acc=$(basename "$f" .gbk)
    echo -e "$acc\t$(grep -c '>' "CDS/$acc.txt")\t$(grep -c '^     CDS ' "$f")"
done
