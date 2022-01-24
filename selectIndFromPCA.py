import math

def calcMissing(file2):
    data = open(file2).read().split("\n")
    miss = {}
    for d in data[2:]:
        d = d.split("\t")
        if len(d) > 1:
            ind = d[0].split(",")[0].strip(" ")
            d[0] = d[0].split(",")[1].strip(" ")
            dmis = [da for da in d if da == "0000"]
            miss[ind] = len(dmis)/len(d)
    return(miss)

def loadPCAInds(file3):
    data = open(file3).read().split("\n")
    return(data)

def loadPCA(file1,inds):
    data = open(file1).read().split("\n")
    pca = {}
    for p in range(len(data)):
        if len(data[p]) > 3:
            ind = inds[p]
            pca1 = float(data[p].split(" ")[0])
            pca2 = float(data[p].split(" ")[1])
            pca[ind] = [pca1,pca2]
    return(pca)

def findCenter(pca,p1gx,p2gx):
    inInds = []
    ax1 = []
    ax2 = []
    for ind,vals in pca.items():
        if p1gx[0] < vals[0] < p1gx[1]:
            if p2gx[0] < vals[1] < p2gx[1]:
                inInds.append(ind)
                ax1.append(vals[0])
                ax2.append(vals[1])
    center = [sum(ax1)/len(ax1),sum(ax2)/len(ax2)]
    return(center,inInds)

def calcDistC(pos,center):
        dist = math.sqrt(pow(center[0] - pos[0],2) + pow(center[1] - pos[1],2))
        return(dist)

def calcIndsDistAcc(inInds,center,pca,propOfInd):
    distall = []
    indall = []
    for ind in inInds:
        pos = pca.get(ind)
        dist = calcDistC(pos,center)
        distall.append(dist)
        indall.append(ind)

    indin = [x for _,x in sorted(zip(distall,indall))]
    thres = int(round(len(distall)*propOfInd,0))
    indin = indin[:thres]
    return(indin)

def indWithLowestMiss(indin,miss):
    misind = []
    for i in indin:
        misind.append(miss.get(i))
    indacc = [x for _,x in sorted(zip(misind,indin))][0:10]
    return(indacc)

file1 =  "lea\pca\m3M3n3_snps_c8_maf05_mds3_mdi03_fi_si.pca\m3M3n3_snps_c8_maf05_mds3_mdi03_fi_si.projections"
file2 = "genepop\m3M3n3_snps_c8_maf05_mds3_mdi03_fi_si.pop"
file3 = "lfmm\m3M3n3_snps_c8_maf05_mds3_mdi03_fi_si.lfmm_ind"
output = "listSamplesGenomeContent.txt"

p1g1 = [-20,0]
p2g1 = [-5,10]
p1g2 = [10,25]
p2g2 = [-5,10]

propOfInd = 0.6

miss = calcMissing(file2)
inds = loadPCAInds(file3)
pca = loadPCA(file1,inds)
center1,inInds1 = findCenter(pca,p1g1,p2g1)
center2,inInds2 = findCenter(pca,p1g2,p2g2)
indin1 = calcIndsDistAcc(inInds1,center1,pca,propOfInd)
indin2 = calcIndsDistAcc(inInds2,center2,pca,propOfInd)

indacc1 = indWithLowestMiss(indin1,miss)
indacc2 = indWithLowestMiss(indin2,miss)

out = open(output,"w")
out.write("Cluster1\t" + "\t".join(indacc1))
out.write("\nCluster2\t" + "\t".join(indacc2))
out.close()
