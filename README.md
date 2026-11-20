# Orthology tool comparison on plant chloroplast genomes

Four orthology and pangenome tools (Proteinortho, Panaroo, PPanGGOLiN and Roary) run on the same set of ten chloroplast genomes from three angiosperm families. The question is simple: which tool recovers the right gene families when the genomes are distantly related?

Chloroplast genomes make a useful test case because the answer is largely known. Genes such as *psbA*, *rbcL* and *matK* carry the same names in basil, wheat and potato, so the annotations act as a ground truth. A good tool should put all ten copies of *psbA* into one cluster.

## Dataset

| Group | Species | Accession |
|---|---|---|
| *Ocimum* (Lamiaceae) | *O. basilicum* | NC_035143.1 |
| | *O. tenuiflorum* | NC_043873.1 |
| | *O. gratissimum* | MW348919.1 |
| *Aconitum* (Ranunculaceae) | *A. carmichaelii* | KY006977 |
| | *A. kusnezoffii* | NC_031422 |
| | *A. ciliare* | NC_031420 |
| Poaceae | wheat | NC_002762 |
| | rice | NC_001320 |
| | maize | NC_001666 |
| Solanaceae | potato | NC_008096.2 |

A second panel of 16 *Ocimum* plastomes was run through Panaroo, Proteinortho and PPanGGOLiN as a within-genus comparison (accessions in `accessions_ocimum16.txt`).

## Pipeline

```
NCBI accessions
 ├── GenBank (.gbk) ──> gbk2gff.py ──> GFF3 ──> Panaroo, PPanGGOLiN, Roary
 ├── protein FASTA ──> gap characters stripped ──> Proteinortho, hmmscan (Pfam-A)
 └── whole plastome FASTA ──> MAFFT ──> IQ-TREE 2 / RAxML-NG
```

All tools read from the same NCBI records, so differences in results come from the tools rather than from file preparation.

| Script | Step |
|---|---|
| `scripts/01_download.sh` | Fetch GenBank, protein FASTA and genome FASTA for every accession |
| `scripts/02_prepare_inputs.sh` | Convert GenBank to GFF3, strip gap characters from protein FASTA |
| `scripts/gbk2gff.py` | GenBank to GFF3 converter with a product-name fallback (see below) |
| `scripts/03_proteinortho.sh` | Proteinortho with BLASTP+ |
| `scripts/04_panaroo.sh` | Panaroo, sensitive clean mode |
| `scripts/05_ppanggolin.sh` | PPanGGOLiN full workflow |
| `scripts/06_roary.sh` | Roary at default (95%) and relaxed (50%) identity |
| `scripts/07_pfam_domains.sh` | hmmscan against Pfam-A, filtered at E < 1e-5 |
| `scripts/08_phylogeny.sh` | Whole-plastome alignment and maximum likelihood trees |

## Results (ten-genome panel)

| Metric | Proteinortho | Panaroo | PPanGGOLiN | Roary |
|---|---|---|---|---|
| Clusters or gene families | 120 | 152 | 205 | 317 |
| Genes placed | 873 | 768 | 873 | 749 |
| Present in all ten genomes | 51 | 37 | 34 | 13 |
| Singleton clusters | 32 | 27 | 78 | 155 |
| Clusters containing paralogs | 25 | 0 | 39 | 0 |
| Mean cluster size | 7.28 | 5.05 | 4.26 | 2.36 |

**Main findings**

- **Proteinortho performs best across families.** Its adaptive identity threshold (down to 25%) recovered 51 core families, closest to the 55 to 70 shared protein-coding genes expected in angiosperm plastomes.
- **Bacterial pangenome tools fragment cross-family data.** Panaroo (98% identity), Roary (95%) and PPanGGOLiN were built for strains of one species. Across monocots and eudicots they split single genes into many clusters, and Roary left almost half its clusters as single-genome orphans. The PPanGGOLiN frequency plot peaks at exactly three genomes, matching the taxonomic group sizes rather than the biology.
- **The core sets are nested.** Every gene in the Roary core appears in the Panaroo and PPanGGOLiN cores, so the tools differ in how much they recover, not in what they consider conserved.
- **Panaroo's core genes still carry correct phylogenetic signal.** The IQ-TREE tree from its core alignment recovers all four groups at 100% bootstrap.
- **Annotation quality matters as much as the tool.** KY006977 is a direct GenBank submission with only 61 CDS (others carry 83 to 111) and no gene name qualifiers. It is missing 33 Pfam families, including RuBisCO and ATP synthase domains, which points to under-annotation rather than gene loss.
- **Domain architecture is almost invariant** across the panel, apart from whole-gene losses in Poaceae (*accD*, *ycf1*, *ycf2*) and an *Ocimum*-specific difference in *ycf1* pseudogene copy number at the inverted repeat boundary.

**Which tool to use when**

- Orthology across genera or families: Proteinortho (or OrthoFinder).
- Pangenome within one species, uncertain annotation quality: Panaroo.
- Large populations, fragmented assemblies, or interest in the core to accessory gradient: PPanGGOLiN.
- Roary: historical reference only; Panaroo supersedes it.

## Problems solved along the way

- **Missing gene names.** KY006977 has no `gene=` qualifiers, so Panaroo treated its genes as unnamed. `gbk2gff.py` falls back to the product description and maps common plastid products to gene symbols (for example, "photosystem II protein D1" to *psbA*). It also warns when CDS counts vary by more than 30% between genomes.
- **Gap characters** in the NCBI protein FASTA broke Proteinortho; they are stripped in `02_prepare_inputs.sh`.
- **panX** was dropped because its codebase requires Python 2, which is no longer packaged on current Ubuntu. Roary replaced it, which gives a cleaner contrast anyway: Panaroo is effectively Roary plus an error-correction layer.
- **Prokka was deliberately not used.** It is built for bacteria and would relabel plastid genes as generic or hypothetical, destroying the gene-name ground truth.

## Environment

Ubuntu under WSL. Panaroo and PPanGGOLiN in separate Python virtual environments; Proteinortho, BLAST+, Roary, HMMER, MAFFT and IQ-TREE 2 installed system-wide; RAxML-NG from the precompiled binary.

```bash
sudo apt install proteinortho ncbi-blast+ roary hmmer mafft iqtree ncbi-entrez-direct mmseqs2
python3 -m venv ~/panaroo_env && source ~/panaroo_env/bin/activate && pip install panaroo biopython && deactivate
python3 -m venv ~/ppanggolin_env && source ~/ppanggolin_env/bin/activate && pip install ppanggolin && deactivate
```

Tool versions used: Proteinortho 6.3.1, BLAST+ 2.17.0, PPanGGOLiN 2.3.1, RAxML-NG 1.2.2.

## Running it

```bash
bash scripts/01_download.sh
bash scripts/02_prepare_inputs.sh
bash scripts/03_proteinortho.sh
bash scripts/04_panaroo.sh
bash scripts/05_ppanggolin.sh
bash scripts/06_roary.sh
bash scripts/07_pfam_domains.sh     # needs Pfam-A.hmm, see script
bash scripts/08_phylogeny.sh
```
