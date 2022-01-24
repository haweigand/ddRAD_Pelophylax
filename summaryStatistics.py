import os

def openFile(file):
    data = {}
    l = 0
    for line in open(file):
        line2 = line.replace("\r","").strip("\n")
        l += 1
        if l == 2:
            loci = line2.split(",")
        elif l > 2:
            if line2 != "POP":
                line3 = line2.replace(",","\t").split("\t")
                line3 = [l.strip(" ") for l in line3]
                data[line3[0]] = line3[1:]
    return(data,loci)

def calcStat(gt):
    if gt[0:2] == gt[2:]:
        return("ho")
    else:
        return("he")


def calcIndividualGenotypes(data):
    gtyps = {}
    for ind,locs in data.items():
        locs2 = [l for l in locs if l != "0000"]
        locshe = list(filter(lambda x: calcStat(x) == "he",locs2))
        locsho = list(filter(lambda x: calcStat(x) == "ho",locs2))
        gtyps[ind] = [len(locs),len(locs2),len(locsho),len(locshe)]
    return(gtyps)

def testCases(vals):
    idvals = list(set(vals))
    if len(idvals) == 1:
        if idvals[0] == "0000":
            return("amis")
        elif idvals[0][0:2] == idvals[0][2:]:
            return("homo")
        else:
            return("hetero")
    elif len(idvals) == 2:
        if "0000" in idvals:
            return("pmis")
        else:
            g1 = calcStat(vals[0])
            g2 = calcStat(vals[1])
            if g1 == "ho" and g2 == "ho":
                return("homhom")
            elif (g1 == "ho" and g2 == "he") or (g1 == "he" and g2 == "ho"):
                return("hethom")
            else:
                if (g1[0:2] == g2[2:]) and (g1[2:] == g2[0:2]):
                    return("hetero")
                else:
                    return("hethet")
    else:
        if "0000" in idvals:
            return("pmis")
        else:
            return("comp")

def calcDifferentAlleles(data):
    names = data.keys()
    inds = set([i[:-1] for i in names])
    types = {}
    for i in inds:
        testinds = list(filter(lambda x: x.startswith(i),names))
        ti_data = []
        if len(testinds) > 1:
            for ti in testinds:
                ti_data.append(data.get(ti))
            cases = []
            for p in range(0,len(ti_data[0])):
                vals = []
                for v in range(0,len(ti_data)):
                    vals.append(ti_data[v][p])
                cases.append(testCases(vals))
            res = {}
            for c in set(cases):
                res[c] = cases.count(c)
            types[i] = res
    return(types)

def makeHeaderSummary(out,all_inds,all_pairs,thresh):
    out.write("Name\tmds\tmdi\tLoci\tIndividuals")
    out.write("\tAvMissing\tAvHetero")
    out.write("\tAvWrong\tThreshWrongInds_" + str(thresh))
    return(out)

def makeHeaderInds(out,all_inds,all_pairs,thresh):
    out.write("Dataset\t\t\t% missing")
    out.write("\t".join(["" for i in all_inds]))
    out.write("\t% homo")
    out.write("\t".join(["" for i in all_inds]))
    out.write("\t% hetero")
    out.write("\t".join(["" for i in all_inds]))
    out.write("\t%wrong")
    out.write("\t".join(["" for i in all_pairs]))
    out.write("\t%He-He")
    out.write("\t".join(["" for i in all_pairs]))
    out.write("\t%Ho-Ho")
    out.write("\t".join(["" for i in all_pairs]))
    out.write("\t%He-Ho")
    out.write("\t".join(["" for i in all_pairs]))

    out.write("\nName\tmds\tmdi")
    out.write("\t" + "\t".join(all_inds))
    out.write("\t" + "\t".join(all_inds))
    out.write("\t" + "\t".join(all_inds))
    out.write("\t" + "\t".join(all_pairs))
    out.write("\t" + "\t".join(all_pairs))
    out.write("\t" + "\t".join(all_pairs))
    out.write("\t" + "\t".join(all_pairs))
    return(out)


def write_ind_data(all_inds,gtypes,t,out):
    for ind in all_inds:
        all = gtypes.get(ind,"NA")
        if all == "NA":
            out.write("\tNA")
        else:
            if t == "m":
                out.write("\t" + str(round((all[0]-all[1])/all[0],4)))
            elif t == "ho":
                out.write("\t" + str(round(all[2]/all[0],4)))
            elif t == "he":
                out.write("\t" + str(round(all[3]/all[0],4)))

def write_type_data(all_pairs,types,t,out,thresh,othres,co):
    for ind in all_pairs:
        all = types.get(ind,"NA")
        if all == "NA":
            out.write("\tNA")
        else:
            co += 1
            tot = sum([v for v in all.values()])
            if t == "w":
                wr = 0
                wr += all.get("hethom",0)
                wr += all.get("hethet",0)
                wr += all.get("homhom",0)
                wr += all.get("comp",0)
                out.write("\t" + str(round(wr/tot,4)))
                if wr/tot > thresh:
                    othres += 1
            elif t == "homhom":
                out.write("\t" + str(round(all.get("homhom",0)/tot,4)))
            elif t == "hethet":
                out.write("\t" + str(round(all.get("hethet",0)/tot,4)))
            elif t == "hethom":
                out.write("\t" + str(round(all.get("hethom",0)/tot,4)))
    return(othres,co)

def getWrong(types,thresh):
    wrong = 0
    thr = 0
    all = 0
    mis = 0
    for ind,cat in types.items():
        w = 0
        a = 0
        for c,v in cat.items():
            all += v
            a += v
            if c == "homhom" or c == "hethom" or c == "hethet" or c == "comp":
                wrong += v
                w += v
            elif c == "amis" or c == "pmis":
                mis += v
        if w/a >= thresh:
            thr += 1
    return(thr,all,wrong,mis)

def generateOutputSummary(results,output,all_inds,all_pairs,thresh):
    out = open(output,"w")
    out = makeHeaderSummary(out,all_inds,all_pairs,thresh)

    names = sorted(results.keys())
    for na in names:
        nab = na.replace("p_m","pm")
        na2 = nab.split("_")
        out.write("\n" + na2[0] + "\t" + na2[4][3:] + "\t" + na2[5][3:])
        gtypes = results.get(na)[0]
        types = results.get(na)[1]
        inds = results.get(na)[2]
        out.write("\t{}\t{}".format(gtypes.get(list(gtypes.keys())[0])[0],inds))
        totalloci = sum([dat[0] for ind,dat in gtypes.items()])
        totalmis = sum([dat[0] - dat[1] for ind,dat in gtypes.items()])
        totalhet = sum([dat[2] for ind,dat in gtypes.items()])
        out.write("\t{}\t{}".format(round(totalmis/totalloci,4),round(totalhet/(totalloci-totalmis),4)))
        thr,all,wrong,mis = getWrong(types,thresh)
        out.write("\t{}\t{}".format(round(wrong/all,4),thr))

def generateOutputInds(results,output,all_inds,all_pairs,thresh):
    out = open(output,"w")
    out = makeHeaderInds(out,all_inds,all_pairs,thresh)

    names = sorted(results.keys())
    for na in names:
        nab = na.replace("p_m","pm")
        na2 = nab.split("_")
        out.write("\n" + na2[0] + "\t" + na2[4][3:] + "\t" + na2[5][3:])
        gtypes = results.get(na)[0]
        types = results.get(na)[1]
        inds = results.get(na)[2]
        othres = 0
        co = 0
        write_ind_data(all_inds,gtypes,"m",out)
        write_ind_data(all_inds,gtypes,"ho",out)
        write_ind_data(all_inds,gtypes,"he",out)
        othres,co = write_type_data(all_pairs,types,"w",out,thresh,othres,co)
        write_type_data(all_pairs,types,"hethet",out,thresh,othres,co)
        write_type_data(all_pairs,types,"homhom",out,thresh,othres,co)
        write_type_data(all_pairs,types,"hethom",out,thresh,othres,co)
    out.close()

def openFiles(folder,ending,output,output2,tresh):
    results = {}
    all_inds = []
    all_pairs = []
    dircont = os.listdir(folder)
    dircont = [f for f in dircont if f.endswith(ending)]
    for file in dircont:
        data,loci = openFile(folder + "/" + file)
        inds = len(list(data.keys()))
        all_inds = list(set(all_inds + list(data.keys())))
        all_pairs = list(set(all_pairs + [i[:-1] for i in list(data.keys())]))
        gtyps = calcIndividualGenotypes(data)
        types = calcDifferentAlleles(data)
        results[file] = [gtyps,types,inds]
    generateOutputSummary(results,output,sorted(all_inds),sorted(all_pairs),thresh)
    generateOutputInds(results,output2,sorted(all_inds),sorted(all_pairs),thresh)

folder = "genepop"
ending = "_fi_si.pop"
output = "statistics_genepop_summary.txt"
output2 = "statistics_genepop_individuals.txt"
thresh = 0.03
openFiles(folder,ending,output,output2,thresh)
