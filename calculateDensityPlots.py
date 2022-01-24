import seaborn as sns
import matplotlib.pyplot as plt
import os

def loadInds(file):
    inds = []
    data = open(file).read().replace("\r\n","\n").split("\n")
    for line in data:
        line = line.strip(" ")
        if line != "":
            inds.append(line)
    return(inds)

def calc_cov_dist(val,ind,dens_all,dens_L,dens_R,refa,fix_all,fix_L,fix_R,covlim):
    if len(val) > 1:
        cov = int(val[1])
        if cov >= covlim:
            gt = [int(i) for i in val[0].split("/")]
            if gt[0] == refa and gt[1] != refa:
                pl = 0
                pr = 1
            elif gt[0] != refa and gt[1] == refa:
                pl = 1
                pr = 0
            if sum(gt) == 1:
                ac = val[2].split(",")
                old = dens_all.get(ind,[])
                old.append(int(ac[pl])/cov)
                old.append(int(ac[pr])/cov)
                dens_all[ind] = old
                old = dens_L.get(ind,[])
                old.append(int(ac[pl])/cov)
                dens_L[ind] = old
                old = dens_R.get(ind,[])
                old.append(int(ac[pr])/cov)
                dens_R[ind] = old
            elif sum(gt) == 0:
                fix_all[ind] = fix_all.get(ind,0) + 1
                if refa == 0:
                    fix_L[ind] = fix_L.get(ind,0) + 1
                elif refa == 1:
                    fix_R[ind] = fix_R.get(ind,0) + 1
            elif sum(gt) == 2:
                fix_all[ind] = fix_all.get(ind,0) + 1
                if refa == 0:
                    fix_R[ind] = fix_R.get(ind,0) + 1
                elif refa == 1:
                    fix_L[ind] = fix_L.get(ind,0) + 1

    return(dens_all,dens_L,dens_R,fix_all,fix_L,fix_R)

def setRef(a1,a2,ref_allele):
    if ref_allele == a1:
        return(0)
    elif ref_allele == a2:
        return(1)
    else:
        print("wrong reference")

def get_data(file1,rflocs,inds,covlim):
    dens_all = {}
    dens_L = {}
    dens_R = {}
    fix_all = {}
    fix_L = {}
    fix_R = {}
    data = open(file1).read().split("\n")
    for line in data:
        if len(line) > 5:
            if line.startswith("#CHROM"):
                inds2 = line.rstrip("\n").split("\t")[9:]
            elif line.startswith("#") == False:
                line2 = line.rstrip("\n").split("\t")
                loc = line2[0]
                pos = line2[1]
                na = str(loc) + "_" + str(pos)
                if na in rflocs.keys():
                    refa = setRef(line2[3],line2[4],rflocs.get(na))
                    for i in range(0,len(inds2)):
                        ind = inds2[i]
                        if ind in inds:
                            val = line2[9+i].split(":")
                            dens_all,dens_L,dens_R,fix_all,fix_L,fix_R = calc_cov_dist(val,ind,dens_all,dens_L,dens_R,refa,fix_all,fix_L,fix_R,covlim)
    return(dens_all,dens_L,dens_R,fix_all,fix_L,fix_R)

def make_density_plots(inds,dens_all,dens_L,dens_R,fix_all,fix_L,fix_R,folder,name,name2):
    sns.set_style('darkgrid')
    out2 = open(folder + "/" + name2, "w")
    out2.write("Specimen\tVariable\tFixed")
    for ind in sorted(inds):
        va = dens_all.get(ind,[])
        vl = dens_L.get(ind,[])
        vr = dens_R.get(ind,[])
        fa = fix_all.get(ind,0)
        fl = fix_L.get(ind,0)
        fr = fix_R.get(ind,0)
        fix,ax1 = plt.subplots()
        ax1 = sns.kdeplot(va,bw=0.1,color="black")
        ax1 = sns.kdeplot(vl,bw=0.1,color="green")
        ax1 = sns.kdeplot(vr,bw=0.1,color="brown")
        ax1.set_xlim(0,1)
        ax1.yaxis.set_label_position("left")
        ax1.set_ylabel("Density")

        ax2 = ax1.twinx()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel("Proportion of homozygous loci",rotation=270,labelpad=10)
        if (fl+fr+len(va)) > 0:
            ax2.bar(0.95,fl/(fl+fr+len(va)/2),width=0.1,color="green")
            ax2.bar(0.95,fr/(fl+fr+len(va)/2),bottom = fl/(fl+fr+len(va)/2),width=0.1,color="brown")
        ax2.set_ylim(0,1)
        ax2.set_xlim(0,1)

        plt.savefig(folder + "/" + ind + "_" + name + "_" + str(len(va)/2) + "_" + str(fa) + ".svg")
        plt.savefig(folder + "/" + ind + "_" + name + "_" + str(len(va)/2) + "_" + str(fa) + ".png")
        plt.close()
        out2.write("\n{}\t{}\t{}".format(ind,len(va)/2,fa))
    out2.close()

def load_ref(file):
    rflocs = {}
    for line in open(file).read().split("\n"):
        line2 = line.split("\t")
        if len(line2) == 2:
            rflocs[line2[0]] = line2[1]
    return(rflocs)


indfile = "remaining_inds_final.txt"
reflocs =  "temptest.reflocs"

folder = "density_temp2"

vcffile = "vcf/m3M3n3_snps.vcf"

name = "remaining_inds"
name2 = "remaining_inds.txt"

covlim = 8

try:
    os.stat(folder)
except:
    os.mkdir(folder)

inds = loadInds(indfile)
rflocs = load_ref(reflocs)
print(inds)
print(rflocs)

dens_all,dens_L,dens_R,fix_all,fix_L,fix_R = get_data(vcffile,rflocs,inds,covlim)
print(dens_all)
print(fix_all)
make_density_plots(inds,dens_all,dens_L,dens_R,fix_all,fix_L,fix_R,folder,name,name2)
