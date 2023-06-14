import sys

fpin = open("GCF_027409825.1_fScoJap1.pri_genomic.fna",'r')
fpout = open("GCF_027409825.1_fScoJap1.pri_genomic.upper.fna",'w')
iFlag_start = 0
for line in fpin:
    if line[0]==">":
        if iFlag_start == 0:
            fpout.write(line)
            iFlag_start = 1
        else:
            fpout.write("\n")
            fpout.write(line)
    else:
        fpout.write(line.strip().upper())
fpout.write('\n')
fpin.close()
fpout.close()
print("End of upper")

nScaff_nSeq_dic = {}
nScaff_list = []
fpin = open("GCF_027409825.1_fScoJap1.pri_genomic.upper.fna","r")
fpout= open("GCF_027409825.1_fScoJap1.pri_genomic.AACCCT.bed",'w')
for line in fpin:
  if line[0]==">":
    nScaff  = line[1:].strip().split(" ")[0]
    #nScaff_nSeq_dic.setdefault(nScaff,'')
    #nScaff_list.append(nScaff)
    print("Reading: ",nScaff)
  else:
    nSeq = line.strip()
    for i in range(0,len(nSeq)-6):
      tmpSeq = nSeq[i:i+6]
      if tmpSeq == "AACCCT":
        tmpline = nScaff + '\t' +str(i) +'\t'+ str(i+6) +'\t'+ "AACCCT" +'\t'+"1"+"\t"+"+"+'\n'
        fpout.write(tmpline)
      elif tmpSeq == "AGGGTT":
        tmpline = nScaff + '\t' +str(i) +'\t'+ str(i+6) +'\t'+ "AGGGTT" +'\t'+"1"+"\t"+"-"+'\n'
        fpout.write(tmpline)
    #nScaff_nSeq_dic[nScaff]+= nSeq
fpin.close()
fpout.close()

fpin = open("GCF_027409825.1_fScoJap1.pri_genomic.AACCCT.bed",'r')
tmp_nScaff = ""
tmp_nStart = ""
tmp_nEnd = ""
tmp_nStrand = ""
tmp_NAME = ""
tmp_lenBlock = 1
nScaff_nStart_info_dic = {}
nScaff_nStartList_dic = {}
for line in fpin:
    part= line.strip().split("\t")
    nScaff = part[0]
    nStart = part[1]
    nEnd   = part[2]
    nNAME  = part[3]
    nStrand = part[5]

    if tmp_nScaff == nScaff:
        if tmp_nEnd == nStart:
            tmp_nEnd = nEnd
            tmp_nStrand = nStrand
            tmp_NAME = nNAME
            tmp_lenBlock += 1

        else:
            if tmp_lenBlock >= 1:
                nScaff_nStartList_dic.setdefault(nScaff,[])
                nScaff_nStartList_dic[nScaff].append(tmp_nStart)
                nScaff_nStart_info_dic.setdefault(nScaff,{})
                nScaff_nStart_info_dic[nScaff].setdefault(tmp_nStart,[tmp_nEnd,tmp_lenBlock,tmp_NAME,tmp_nStrand])
            tmp_nScaff = nScaff
            tmp_nStart = nStart
            tmp_nEnd   = nEnd
            tmp_nStrand = nStrand
            tmp_NAME = nNAME
            tmp_lenBlock = 1
    else:
        if tmp_lenBlock >= 1:
            nScaff_nStartList_dic.setdefault(nScaff,[])
            nScaff_nStartList_dic[nScaff].append(tmp_nStart)
            nScaff_nStart_info_dic.setdefault(nScaff,{})
            nScaff_nStart_info_dic[nScaff].setdefault(tmp_nStart,[tmp_nEnd,tmp_lenBlock,tmp_NAME,tmp_nStrand])
        tmp_nScaff = nScaff
        tmp_nStart = nStart
        tmp_nEnd   = nEnd
        tmp_nStrand = nStrand
        tmp_NAME = nNAME
        tmp_lenBlock = 1


fpout = open('GCF_027409825.1_fScoJap1.pri_genomic.AACCCTn_nOver1.bed','w')
for nScaff in nScaff_nStartList_dic.keys():
    for nStart in nScaff_nStartList_dic[nScaff]:
        nEnd  = nScaff_nStart_info_dic[nScaff][nStart][0]
        nNAME = "(AACCCT)"+str(nScaff_nStart_info_dic[nScaff][nStart][1])
        nScore = str(nScaff_nStart_info_dic[nScaff][nStart][1])
        nStrand = str(nScaff_nStart_info_dic[nScaff][nStart][3])
        tmpline = nScaff+"\t"+nStart+"\t"+nEnd+"\t"+nNAME+'\t'+nScore+'\t'+nStrand+'\t'+"\n"
        fpout.write(tmpline)
fpout.close()
