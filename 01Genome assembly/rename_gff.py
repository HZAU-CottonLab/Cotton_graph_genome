#!/usr/bin/env python3
import re
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python rename.py input_file output_prefix")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

    with open(input_file, 'r') as gff:
        count = 0
        mRNA_count = 0
        cds_count = 0
        exon_count = 0

        for line in gff:
            if not line.startswith("\n"):
                records = line.split("\t")
                records[1] = "EVM"

            if re.search(r"\tgene\t", line):
                count += 10
                mRNA_count = 0
                gene_id = output_prefix + str(count).zfill(6)
                records[8] = f"ID={gene_id}"
            elif re.search(r"\tmRNA\t", line):
                cds_count = 0
                exon_count = 0
                mRNA_count += 1
                mRNA_id = f"{gene_id}.{mRNA_count}"
                records[8] = f"ID={mRNA_id};Parent={gene_id}"
            elif re.search(r"\texon\t", line):
                exon_count += 1
                exon_id = f"{mRNA_id}_exon_{exon_count}"
                records[8] = f"ID={exon_id};Parent={mRNA_id}"
            elif re.search(r"\tCDS\t", line):
                cds_count += 1
                cds_id = f"{mRNA_id}_cds_{cds_count}"
                records[8] = f"ID={cds_id};Parent={mRNA_id}"
            else:
                continue

            print("\t".join(records), end='')

if __name__ == "__main__":
    main()
