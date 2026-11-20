#!/usr/bin/env bash
# Build tool inputs from the downloaded NCBI records.
set -euo pipefail

# GenBank -> GFF3 (Panaroo, PPanGGOLiN, Roary). Needs biopython.
source ~/panaroo_env/bin/activate
python scripts/gbk2gff.py gbks gffs
deactivate

# Check every genome has psbA after conversion (should print 1 for each)
grep -c "gene=psbA" gffs/*.gff

# Protein FASTA: strip gap characters from sequence lines (they break Proteinortho)
mkdir -p CDS_clean
for f in CDS/*.txt; do
    sed '/^>/! s/-//g' "$f" > "CDS_clean/$(basename "$f")"
done

# Remove Windows Zone.Identifier files if any were copied in from Windows
find . -name "*Zone.Identifier" -delete
