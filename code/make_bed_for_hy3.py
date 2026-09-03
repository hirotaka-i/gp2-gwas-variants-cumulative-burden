#!/usr/bin/env python3

import csv
import re
import argparse

GENE_STRAND = {
    "APOE": "+",
    "GBA1": "-",
    "LRRK2": "+",
}


def parse_coordinate(coord: str):
    """'Chr19:44908684' -> ('chr19', 44908684)  1-based"""
    m = re.match(r"\s*chr?([0-9XYMT]+)\s*:\s*([\d,]+)\s*", coord, re.IGNORECASE)
    if not m:
        raise ValueError(f"Cannot detect the position: {coord!r}")
    chrom = m.group(1).upper().replace("MT", "M")
    pos = int(m.group(2).replace(",", ""))
    return chrom, pos


def parse_alleles(allele: str):
    """'T > C' -> ('T', 'C')"""
    parts = re.split(r"\s*>\s*", allele.strip())
    if len(parts) != 2:
        raise ValueError(f"Cannot parse Ref/Alt alleles: {allele!r}")
    return parts[0].strip(), parts[1].strip()


def convert(in_path: str, out_path: str):
    rows_out = []
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        # Strip whitespace from header fields
        reader.fieldnames = [h.strip() for h in reader.fieldnames]

        for row in reader:
            row = {k.strip(): (v.strip() if v else v) for k, v in row.items()}

            gene = row["Gene"]
            variant = row["Variant"]
            rsid = row["dbSNP ID"]
            chrom, pos_1based = parse_coordinate(row["GRCh38 Coordinate"])
            ref, alt = parse_alleles(row["Ref/Alt Allele (Genomic)"])

            start = pos_1based    # 1-based
            end = pos_1based          
            name = f"{gene}|{variant}|{rsid}|{ref}>{alt}"
            score = "0"
            strand = GENE_STRAND.get(gene, ".")

            rows_out.append([chrom, str(start), str(end), name, score, strand])

    # Sort by coordinates (chrom, start)
    def chrom_key(c):
        num = c[3:]
        return (0, int(num)) if num.isdigit() else (1, num)

    rows_out.sort(key=lambda r: (chrom_key(r[0]), int(r[1])))

    with open(out_path, "w", encoding="utf-8") as f:
        # Enable trackline if needed:
        # f.write('track name="genetic_variants" description="GRCh38 variants"\n')
        for r in rows_out:
            f.write("\t".join(r) + "\n")

    print(f"Wrote {len(rows_out)} rows -> {out_path}")


convert("variant_list/curated_variants_hy3.csv", "temp/curated_variants_hy3.bed")

# combine multiple BED files (keep the first 3 columns)
## remove duplicate entries and sort
bed1 = "temp/curated_variants_hy3.bed"
bed2 = "temp/Leonard_risk.bed"
out_bed = "temp/gp2_r12_hy3_combined.bed"

unique_rows = set()
for bed in [bed1, bed2]:
    with open(bed, "r", encoding="utf-8") as fin:
        for line in fin:
            cols = line.strip().split("\t")
            unique_rows.add(tuple(cols[:3]))

with open(out_bed, "w", encoding="utf-8") as fout:
    for row in sorted(unique_rows):
        fout.write("\t".join(row) + "\n")

print(f"Combined BED files -> {out_bed}")