import math
import os

def loadInds(file2):
    data = open(file2).read().split("\n")
    c1 = data[0].split("\t")[1:]
    c2 = data[1].split("\t")[1:]
    return(c1,c2)


def loadGenepop(file1,c1,c2):
    data = open(file1).read().split("\n")
    loci = data[1].split(",")
    locs = {}
    alld = {}
    for d in data[2:]:
        d = d.split("\t")
        if len(d) > 1:
            ind = d[0].split(",")[0].strip(" ")
            d[0] = d[0].split(",")[1].strip(" ")
            if ind in c1 or ind in c2:
                locs[ind] = d
            alld[ind] = d
    return(loci,locs,alld)

def getGtLocus(locs,cx,pos):
    gt = []
    for ind in cx:
        g = locs.get(ind)[pos]
        if g != "0000":
            gt.append(g)
    gt = list(set(gt))
    return(gt)


def filterLoci(locs,c1,c2,loci):
    alho1he2 = []
    alho2he1 = []
    for pos in range(len(loci)):
        gt1 = getGtLocus(locs,c1,pos)
        gt2 = getGtLocus(locs,c2,pos)
        if len(gt1) == 1 and len(gt2) == 1:
            if gt1 != gt2:
                if gt1[0][0:2] != gt1[0][2:] and gt2[0][0:2] == gt2[0][2:]:
                    alho2he1.append(loci[pos])
                elif gt1[0][0:2] == gt1[0][2:] and gt2[0][0:2] != gt2[0][2:]:
                    alho1he2.append(loci[pos])

    return(alho1he2,alho2he1)

def makeGenepop(loci,alhohe,alld,output):
    out = open(output,"w")
    out.write("Genotypes fixed between LL and LR\n")
    out.write(",".join(alhohe))
    out.write("\nPOP")
    pos = []
    for lo in alhohe:
        pos.append(loci.index(lo))
    for ind,vals in alld.items():
        out.write("\n{} ,  ".format(ind))
        out.write("\t".join([vals[p] for p in pos]))
    out.close()


file1 = "genepop\m3M3n3_snps_c8_maf05_mds3_mdi03_fi_si.pop"
file2 = "listSamplesGenomeContent.txt"
folder = "fixed"
output1 = "fixed/temptest.pop"
output2 = "fixed/temptest.reflocs"

try:
    os.stat(folder)
except:
    os.mkdir(folder)

c1,c2 = loadInds(file2)
loci,locs,alld = loadGenepop(file1,c1,c2)
alho1he2,alho2he1 = filterLoci(locs,c1,c2,loci)
print(len(alho1he2))
print(len(alho2he1))
makeGenepop(loci,alho2he1,alld,output)
