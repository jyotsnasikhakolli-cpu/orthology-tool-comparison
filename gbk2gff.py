#!/usr/bin/env python3
"""Convert GenBank records to Prokka-style GFF3 for Panaroo, PPanGGOLiN and Roary.

When a CDS has no gene= qualifier, the gene name is taken from the product
description, with common plastid products mapped to standard gene symbols.

Usage: python gbk2gff.py <genbank_dir> <output_dir>
"""
import sys, os, glob
from Bio import SeqIO

PRODUCT_MAP = {
    "photosystem ii protein d1": "psbA",
    "photosystem ii protein d2": "psbD",
    "maturase k": "matK",
    "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit": "rbcL",
    "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit": "rbcL",
    "cytochrome f": "petA",
    "atp synthase cf1 alpha subunit": "atpA",
    "atp synthase cf1 beta subunit": "atpB",
    "atp synthase cf1 epsilon subunit": "atpE",
    "atp synthase cf0 subunit i": "atpF",
    "atp synthase cf0 subunit iii": "atpH",
    "atp synthase cf0 subunit iv": "atpI",
}


def normalise(product):
    key = product.lower().strip()
    if key in PRODUCT_MAP:
        return PRODUCT_MAP[key]
    token = product.split(";")[0].split(",")[0].strip()
    token = token.replace(" ", "_").replace("/", "-").replace("'", "")
    token = "".join(c for c in token if c.isalnum() or c in "_-")
    return token[:40] if token else "unknown"


def convert(gbk_path, out_dir, stats):
    base = os.path.basename(gbk_path)
    for ext in (".gbk", ".gbff", ".gb", ".genbank"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    safe = base.replace(".", "_")
    out_path = os.path.join(out_dir, safe + ".gff")
    records = list(SeqIO.parse(gbk_path, "genbank"))
    if not records:
        print(f"  [SKIP] {gbk_path}: no GenBank records found")
        return
    n_cds = 0
    n_from_product = 0
    with open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for rec in records:
            out.write(f"##sequence-region {rec.id} 1 {len(rec.seq)}\n")
        for rec in records:
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                if "pseudo" in feat.qualifiers or "pseudogene" in feat.qualifiers:
                    continue
                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                strand = "+" if feat.location.strand == 1 else "-"
                n_cds += 1
                product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
                gene = feat.qualifiers.get("gene", [""])[0].strip()
                if not gene:
                    gene = normalise(product)
                    n_from_product += 1
                locus = feat.qualifiers.get("locus_tag", [f"{safe}_{n_cds:05d}"])[0]
                product_clean = product.replace(";", ",").replace("=", "-")
                gene = gene.replace(";", ",").replace("=", "-")
                attrs = (f"ID={safe}_{n_cds:05d};locus_tag={locus}"
                         f";gene={gene};product={product_clean}")
                out.write("\t".join([rec.id, "GenBank", "CDS", str(start), str(end),
                                     ".", strand, "0", attrs]) + "\n")
        out.write("##FASTA\n")
        for rec in records:
            out.write(f">{rec.id}\n")
            seq = str(rec.seq)
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    note = f"  ({n_from_product} named from product)" if n_from_product else ""
    print(f"  {os.path.basename(gbk_path)} -> {os.path.basename(out_path)}  {n_cds} CDS{note}")
    stats.append((safe, n_cds))


def main():
    if len(sys.argv) != 3:
        print("Usage: python gbk2gff.py <genbank_dir> <output_dir>")
        sys.exit(1)
    in_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    files = []
    for ext in ("*.gbk", "*.gbff", "*.gb", "*.genbank"):
        files.extend(glob.glob(os.path.join(in_dir, ext)))
    files.sort()
    if not files:
        print(f"No GenBank files found in {in_dir}")
        sys.exit(1)
    print(f"Converting {len(files)} files...\n")
    stats = []
    for f in files:
        convert(f, out_dir, stats)
    total = sum(n for _, n in stats)
    print(f"\nDone. {total} CDS features across {len(stats)} genomes -> {out_dir}/")
    if stats:
        counts = [n for _, n in stats]
        if max(counts) > min(counts) * 1.3:
            print("\n[!] CDS counts vary a lot between genomes:")
            for name, n in sorted(stats, key=lambda x: x[1]):
                print(f"      {name:<20} {n}")
            print("    Sparse records may be less thoroughly annotated than the rest.")


if __name__ == "__main__":
    main()
