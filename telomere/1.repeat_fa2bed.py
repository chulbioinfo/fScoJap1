#!/usr/bin/env python
# This code modified from the following link
# https://www.danielecook.com/generate-a-bedfile-of-masked-ranges-a-fasta-file/

import sys
fNAME = "GCF_027409825.1_fScoJap1.pri_genomic.fna"
oNAME = "GCF_027409825.1_fScoJap1.pri_genomic_repeat.bed"
chrom = ""
pos = -1
start = -1
in_masked_region = False

fpout = open(oNAME,'w')
with open(fNAME, "r") as fh:
    for line in fh:
        if line.startswith(">"):
            if in_masked_region:  # last masked region from previous chrom
                fpout.write(f"{chrom}\t{start}\t{pos}")
                fpout.write("\n")
                start = -1  # not needed actually
                in_masked_region = False
            pos = 0
            chrom = line.split(" ")[0].replace(">", "").strip()
        else:
            for c in line.strip():
                if not in_masked_region and c.islower() == True:
                    in_masked_region = True
                    start = pos
                elif in_masked_region and c.islower() == False:
                    in_masked_region = False
                    fpout.write(f"{chrom}\t{start}\t{pos}")
                    fpout.write("\n")

                pos += 1

if in_masked_region:  # last masked region in last chrom
    fpout.write(f"{chrom}\t{start}\t{pos}")
    fpout.write("\n")
fpout.close()
