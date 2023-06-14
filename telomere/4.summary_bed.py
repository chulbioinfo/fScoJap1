#!/usr/bin/env python

import sys
import os

flist = ["GCF_027409825.1_fScoJap1.pri_genomic_repeat.bed","fScoJap1_repeat_telomere_AACCCT.bed","GCF_027409825.1_fScoJap1.pri_genomic.AACCCTn_nOver1_loci.bed"]

nScf_iLen_dic = {}
fpin = open("GCF_027409825.1_fScoJap1.pri_genomic.len_scaffold.txt",'r')
for line in fpin:
    part = line.strip().split("\t")
    nScf = part[0].split(" ")[0]
    iLen = int(part[1])
    nScf_iLen_dic.setdefault(nScf,iLen)
fpin.close()


def bed2table(fNAME):

    scafID_len_dic = {}
    scafID_len_dic_5prime_30kbp = {}
    scafID_len_dic_3prime_30kbp = {}
    
    fpin = open(fNAME,'r')
    for line in fpin:
        part = line.strip().split('\t')
        nScf = part[0]
        iStart = int(part[1])
        iEnd = int(part[2])
        iLen = iEnd - iStart
        scafID_len_dic.setdefault(nScf,0)
        scafID_len_dic_5prime_30kbp.setdefault(nScf,0)
        scafID_len_dic_3prime_30kbp.setdefault(nScf,0)
        
        scafID_len_dic[nScf]+=iLen
        if iEnd <= 30000:
            scafID_len_dic_5prime_30kbp[nScf]+=iLen
        if iStart >= (nScf_iLen_dic[nScf]-30000):
            scafID_len_dic_3prime_30kbp[nScf]+=iLen
            
    fpin.close()
    
    oNAME = os.path.basename(fNAME)+".txt"
    fpout = open(oNAME,'w')
    for scafID in scafID_len_dic.keys():
        fpout.write(scafID +'\t'+ str(scafID_len_dic[scafID])+'\n')
    fpout.close()

    oNAME = os.path.basename(fNAME)+"5prime_30kbp.txt"
    fpout = open(oNAME,'w')
    for scafID in scafID_len_dic_5prime_30kbp.keys():
        fpout.write(scafID +'\t'+ str(scafID_len_dic_5prime_30kbp[scafID])+'\n')
    fpout.close()

    oNAME = os.path.basename(fNAME)+"3prime_30kbp.txt"
    fpout = open(oNAME,'w')
    for scafID in scafID_len_dic_3prime_30kbp.keys():
        fpout.write(scafID +'\t'+ str(scafID_len_dic_3prime_30kbp[scafID])+'\n')
    fpout.close()

def main():
    for fNAME in flist:
        bed2table(fNAME)
        
if __name__ == "__main__":
    main()
