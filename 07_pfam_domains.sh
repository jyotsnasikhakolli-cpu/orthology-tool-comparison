#!/usr/bin/env bash
# Pfam domain annotation with HMMER hmmscan, one output per genome.
set -euo pipefail

if [[ ! -f Pfam-A.hmm.h3m ]]; then
    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
fi

mkdir -p hmmscan_out
for f in CDS_clean/*.txt; do
    acc=$(basename "$f" .txt)
    hmmscan --cpu 8 --domtblout "hmmscan_out/$acc.domtbl" Pfam-A.hmm "$f" > /dev/null
done

# Keep domain hits with independent E-value < 1e-5 (column 13 of domtblout)
cat hmmscan_out/*.domtbl | grep -v '^#' | awk '$13 < 1e-5' > hmmscan_out/all_hits_filtered.tsv
wc -l hmmscan_out/all_hits_filtered.tsv
