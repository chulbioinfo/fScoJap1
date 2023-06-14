import sys
fpin = open("GCF_027409825.1_fScoJap1.pri_genomic.fna",'r')
fpout = open("GCF_027409825.1_fScoJap1.pri_genomic.len_scaffold.txt",'w')
scafID_len_dic = {}
for line in fpin:
    if line[0]==">":
        scafID = line[1:].strip()#.split(" ")[0]
        scafID_len_dic.setdefault(scafID,0)
    else:
        scafID_len_dic[scafID]+=len(line.strip())
fpin.close()

for scafID in scafID_len_dic.keys():
    tmpline = scafID +'\t'+ str(scafID_len_dic[scafID])+'\n'
    fpout.write(tmpline)
fpout.close()
