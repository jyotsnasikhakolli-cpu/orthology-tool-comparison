#!/usr/bin/env bash
# Whole-plastome alignment and maximum likelihood trees in IQ-TREE 2 and RAxML-NG.
set -euo pipefail

mkdir -p phylogeny
cat genomes/*.fasta > phylogeny/cp10.fasta
mafft --auto --thread 8 phylogeny/cp10.fasta > phylogeny/cp10.aln

# IQ-TREE 2: ModelFinder (selected TVM+F+G4 under BIC), ultrafast bootstrap + SH-aLRT
iqtree2 -s phylogeny/cp10.aln -m MFP -B 1000 -alrt 1000 -T AUTO --prefix phylogeny/cp10_iqtree

# RAxML-NG with the equivalent model and standard bootstrap
raxml-ng --all --msa phylogeny/cp10.aln --model TVM+FO+G4 --bs-trees 1000 \
         --threads 8 --prefix phylogeny/cp10_raxml
