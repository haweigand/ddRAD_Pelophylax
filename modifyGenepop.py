#!/usr/bin/env python

### Written by Hannah Weigand, University of Duisburg-Essen

### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.

### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import numpy as np
import sys
import argparse
from argparse import RawTextHelpFormatter
if sys.version_info.major == 2:
    from nose.tools import assert_equal
    from nose.tools import assert_items_equal
elif sys.version_info.major == 3:
    from nose.tools import assert_equal as assert_equal
    from nose.tools import assert_count_equal as assert_items_equal

#Define Error messages
class AllIndividualsLost(ValueError):
    pass

class AllLociLost(ValueError):
    pass

class WrongNumberLoci(ValueError):
    pass

class MissingLine(ValueError):
    pass

class WrongTypeLoci(ValueError):
    pass

class NoOverlapWithData(ValueError):
    pass

class WrongLocusListFormat(ValueError):
    pass

class WrongPopSepFormat(ValueError):
    pass

class NoInputDefind(ValueError):
    pass

class NoOutputDefind(ValueError):
    pass

#Define Locus class
class Locus(object):

    #Generate general structure of the Locus class
    def __init__(self, name, position, individuals,order = []):
        self.name = name
        self.position = position
        self.individuals = individuals
        self.order = order


    #Calculate allele counts
    def alleles(self,ind_list):
        counts_alleles = {}
        counts_all = 0
        for ind in ind_list:
            gt0 = self.individuals.get(ind)
            if gt0[0] != "00":
                a1 = gt0[0]
                a2 = gt0[1]
                if self.order == []:
                    self.order = [a1]
                    if a1 != a2:
                        self.order.append(a2)
                else:
                    if self.order.count(a1) < 1:
                        self.order.append(a1)
                    if self.order.count(a2) < 1:
                        self.order.append(a2)
                count_a1 = counts_alleles.get(a1, 0)
                count_a2 = counts_alleles.get(a2, 0)
                counts_alleles[a1] = count_a1 + 1
                counts_alleles[a2] = count_a2 + 1
                counts_all += 2

        return(counts_alleles, counts_all)

    #Calculate basic data statistics for a set of individuals
    def allele_frequency_statistics(self,ind_list):
        counts = {}
        counts_gt_u = {}
        counts_af = {}
        freq = {}
        total_only_data = 0
        total_with_missing = 0
        missing = 0

        for ind in ind_list:
            total_with_missing += 1
            gt0 = self.individuals.get(ind)
            gt = gt0[0] + gt0[1]
            count = counts.get(gt, 0)
            counts[gt] = count + 1

        if ["0000"] == list(counts.keys()):
            return(freq,[],0,1)

        for gt,count in counts.items():
            a1 = gt[0:2]
            a2 = gt[2:]
            gt2 = a2 + a1
            key = list(counts_gt_u.keys())
            if key.count(gt2)==1:
                count_o = counts_gt_u.get(gt2)
                counts_gt_u[gt2] = count_o + count
            else:
                counts_gt_u[gt] = count
            if a1 != "00":
                total_only_data += (2 * count)
                count_af = counts_af.get(a1, 0)
                counts_af[a1] = count_af + count
                count_af = counts_af.get(a2, 0)
                counts_af[a2] = count_af + count

        for allele,count in counts_af.items():
            freq[allele] = count/total_only_data

        gts = []
        gts = [gt for gt in list(counts_gt_u.keys()) if gt != "0000"]
        maf = sorted([fre for allele,fre in freq.items()])[0]
        missing_frequency = counts_gt_u.get("0000", 0)/total_with_missing

        #frequencies, maf, missing_frequency
        return(freq,gts,maf,missing_frequency)


#Define Individual class
class Individual(object):

    #Generate general structure of the Individual class
    def __init__(self, name, population = 0):
        self.name = name
        self.population = population

    #Calculate missing data frequency for a set of loci
    def missing_frequency(self, loc_list):
        missing = 0
        for loc in loc_list:
            gt = loc.individuals.get(self)
            if gt[0] == "00":
                missing += 1

        miss = missing/len(loc_list)
        return(miss)

    #Generate fasta sequence for a set of loci
    def fasta(self, loc_l):
        seq = ""
        for lo in loc_l:
            alleles = lo.individuals.get(self)
            a1 = alleles[0]
            a2 = alleles[1]
            if a1 == a2:
                if a1 == a:
                    seq = seq + "A"
                elif a1 == c:
                    seq = seq + "C"
                elif a1 == g:
                    seq = seq + "G"
                elif a1 == t:
                    seq = seq + "T"
                elif a1 == "00":
                    seq = seq + "N"
            else:
                if (a1 == a and a2 == c) or (a1 == c and a2 == a):
                    seq = seq + "M"
                elif (a1 == a and a2 == g) or (a1 == g and a2 == a):
                    seq = seq + "R"
                elif (a1 == a and a2 == t) or (a1 == t and a2 == a):
                    seq = seq + "W"
                elif (a1 == c and a2 == g) or (a1 == g and a2 == c):
                    seq = seq + "S"
                elif (a1 == c and a2 == t) or (a1 == t and a2 == c):
                    seq = seq + "Y"
                elif (a1 == g and a2 == t) or (a1 == t and a2 == g):
                    seq = seq + "K"

        return(seq)

    #Generate genotypes for a set of loci
    def get_genotype(self, locus):
        alleles = locus.individuals.get(self)
        gt = alleles[0] + alleles[1]
        return(gt)


#Inner Functions
#Filter SNPs by minor allele frequency and missing data frequency (threshold = maf & mds)
#Filter fixed (f = y) SNPs and SNPs with more than two alleles (mta = y)
def filter_loci(ind_l, loc_l, maf, mds,f,mta):
    loc_l_in = []
    maf_ex = 0
    mds_ex = 0
    both_ex = 0
    fixed = 0
    mtas = 0

    #Calculation if fixed loci and loci with more than two alleles are excluded
    if f == True and mta == True:
        for lo in loc_l:
            af, gts, fmaf, fmds = lo.allele_frequency_statistics(ind_l)
            if len(list(af.keys())) == 1:
                fixed += 1
            elif len(list(af.keys())) == 2:
                if len(gts) == 1:
                    fixed += 1
                else:
                    if fmaf > maf:
                        if fmds < mds:
                            loc_l_in.append(lo)
                        else:
                            mds_ex += 1
                    else:
                            if fmds < mds:
                                maf_ex += 1
                            else:
                                both_ex += 1

            else:
                mtas += 1

    #Calculation if fixed loci are excluded
    elif f == True and mta == False:
        for lo in loc_l:
            af, gts, fmaf, fmds = lo.allele_frequency_statistics(ind_l)
            if len(list(af.keys())) == 1:
                fixed += 1
            else:
                if len(gts) == 1:
                    fixed += 1
                else:
                    if fmaf > maf:
                        if fmds < mds:
                            loc_l_in.append(lo)
                        else:
                            mds_ex += 1
                    else:
                            if fmds < mds:
                                maf_ex += 1
                            else:
                                both_ex += 1

    #Calculation loci with more than two alleles are excluded
    elif f == False and mta == True:
        for lo in loc_l:
            af,gts, fmaf, fmds = lo.allele_frequency_statistics(ind_l)
            if len(list(af.keys())) > 2:
                mtas += 1
            else:
                if fmaf > maf:
                    if fmds < mds:
                        loc_l_in.append(lo)
                    else:
                        mds_ex += 1
                else:
                        if fmds < mds:
                            maf_ex += 1
                        else:
                            both_ex += 1

    #Calculation if none of the loci is excluded due to the number of alleles
    else:
        for lo in loc_l:
            af, gts, fmaf, fmds = lo.allele_frequency_statistics(ind_l)
            if fmaf > maf:
                if fmds < mds:
                    loc_l_in.append(lo)
                else:
                    mds_ex += 1
            else:
                    if fmds < mds:
                        maf_ex += 1
                    else:
                        both_ex += 1

    #Error if all loci are excluded
    if len(loc_l_in) == 0:
        raise AllLociLost("Try different filter parameters")

    return(loc_l_in, maf_ex, mds_ex, both_ex,fixed, mtas)

#Filter Individuals by missing data frequency (threshold = mdi)
def filter_individuals(ind_l,loc_l,mdi):
    ind_l_in = []

    for ind in ind_l:
        if ind.missing_frequency(loc_l) < mdi:
            ind_l_in.append(ind)

    #Error if all individuals are excluded
    if len(ind_l_in) == 0:
        raise AllIndividualsLost("Try different filter parameters")
    return(ind_l_in)

#Data loading
def load_data(my_data,slp,mess):
    with open(my_data) as file1:
        pop = 0
        l = 0
        loc_list = []
        ind_list = []
        empty_line = 0
        i = 0
        for line in file1:
            l += 1
            #Get title line
            if l == 1:
                title = line.rstrip("\n")

            else:
                line2 = line.replace(" ","\t").replace(",", "\t").replace("\t\t","\t").rstrip("\n").split("\t")
                #Check for empty lines within the file
                if line2 == [""] or line2 == [" "]:
                    empty_line == l
                else:
                    if empty_line != 0 and l > empty_line:
                        raise MissingLine("A unexcepted empty line is in the file at line " + str(empty_line))

                    #Exclude empty spaces
                    while line2[0] == "":
                        line2 = line2[1:]
                    while line2[-1] == "":
                        line2 = line2[0:-1]

                    #Define loci
                    if l == 2:
                        if len(line2) < 2:
                            raise MissingLine("The locus line causes an error.\nAt least two loci must be present.\nPlease ensure seperation by ',' tab or whitespace.\nLocus line:\n" + str(line2))

                        for lo in line2:
                            if lo.find(slp) > 0:
                                name = lo[0:lo.find(slp)]
                                pos = int(lo[lo.find(slp)+1:])
                            else:
                                name = lo
                                pos = 0
                            loc_list.append(Locus(name = name, position = pos, individuals={}))

                    #Define which population
                    elif len(line2) < 2:
                        pop +=  1

                    #Get individual data
                    else:
                        if pop == 0:
                            raise MissingLine("No Population line before the first Individual")

                        name = line2[0]
                        loci = line2[1:]
                        while loci[0] == "":
                            loci = loci[1:]

                        if len(loci) != len(loc_list):
                            raise WrongNumberLoci("Number of Loci of Locus Line: " + str(len(loc_list)) +
                            "\nNumber of Loci for " + name + ": " + str(len(loci)))

                        ind_list.append(Individual(name=name, population=pop))

                        #Fill locus data with individual allele data
                        for loc in range(len(loci)):
                            locus = loci[loc]
                            a1 = loci[loc][0:2]
                            a2 = loci[loc][2:]
                            if all(x in ["00","01","02","03","04"] for x in (a1,a2)):
                                loc_class = loc_list[loc]
                                ind_class = ind_list[i]
                                loc_class.individuals[ind_class] = (a1, a2)
                            else:
                                raise WrongTypeLoci("Wrong locus:\nIndividual = " + name +
                                "  Locus number = " + str(loc) + "  Locus data = " +  str(locus))
                        i += 1

    mess.append("Analysed loci: " + str(len(loc_list)))
    mess.append("Analysed individuals: " + str(len(ind_list)))
    mess.append("Analysed populations: " + str(pop))
    return(title, pop, loc_list, ind_list, mess)

#Assign new populations to the individuals
#Individuals without an assignment are excluded
def assign_new_pop(ind_l,file_sep,mess):
    new_pop = {}
    ind_in = []
    not_found = []
    l = 0
    with open(file_sep) as file1:
        for line in file1:
            l += 1
            line2 = line.replace(" ","\t").replace(",", "\t").replace("\t\t","\t").replace("\t\t","\t").rstrip("\n").split("\t")
            #Test for write format
            if len(line2) > 2:
                raise WrongPopSepFormat("Error in line " + str(l) + ": " + str(line))

            #Use data for assignment
            elif len(line2) == 2:
                if line2[0] != "" or line2[1] != "" or line2[1] != " ":
                    indiv = [ind for ind in ind_l if ind.name == line2[0]]
                    if len(indiv) == 1:
                        indi = indiv[0]
                        pop = line2[1]
                        if pop != "":
                            ind_in.append(indi)
                            assigned = new_pop.get(pop, [])
                            assigned.append(indi)
                            new_pop[pop] = assigned

                    #Write individual to a list if it is not in the original dataset
                    elif len(indiv) < 1:
                        not_found.append(("line " + str(l) + "; Individual: " + str(line2[0])))

            else:
                if line2[0] != "":
                    indiv = [ind for ind in ind_l if ind.name == line2[0]]
                    if len(indiv) < 1:
                        not_found.append(("line " + str(l) + "; Individual: " + str(line2[0])))

    #Write all individuals which were not part of the input file into the log file
    if len(not_found) > 0:
        mess.append("\nWARNING: Some of the Individuals in the Population assignment file were not part of input data:")
        mess.append([nf for nf in not_found])

    #Error if all individuals are excluded
    if len(ind_in) == 0:
        raise AllIndividualsLost("No individual was assigned to any population.\nPlease ensure that each individual is stored to a new line.\nPlease ensure that individual name and population are separated by the correct spacing symbole: ',' tab or whitespace")

    return(new_pop, ind_in, mess)

#Assign old populations to the individuals (for easier data handling)
def assign_old_pop(ind_l):
    old_pop = {}
    for ind in ind_l:
        pop = str(ind.population)
        assigned = old_pop.get(pop, [])
        assigned.append(ind)
        old_pop[pop] = assigned

    return(old_pop)

#Use user-specified Whitelist
#Include only SNPs which are in the Whitelist as well as in the given locus list
def set_ud_loc_list(loc_l, w_f, slp, mess):
    white = {}
    loc_inc = []
    not_found = []
    l = 0
    for line in open(w_f):
        l += 1
        line2 = line.replace(slp, "\t").replace(" ", "\t").rstrip("\n").split("\t")
        pos = white.get(line2[0], [])

        #Test data format
        if len(line2) > 2:
            raise WrongLocusListFormat("The Whitelist is in a wrong format\nSee line " + str(l) + ": " + str(line))

        elif len(line2) == 2:
            pos.append(line2[1])
        else:
            pos.append(0)
        white[line2[0]] = pos

    for loc, pos_l in white.items():
        for pos in pos_l:
            loci = [lo for lo in loc_l if (str(lo.name) == loc and str(lo.position) == pos)]

            #Test if the locus from the white list is also in the locus list
            if len(loci) == 1:
                loc_inc.append(loci[0])
            elif len(loci) < 1:
                not_found.append(str(loc) + slp + str(pos))

    #Write all loci which were not part of the input file into the log file
    if len(not_found) > 0:
        mess.append("\nWARNING: Some of the Whitelist SNPs were not part of input data:")
        mess.append([lo for lo in not_found])

    #Error if all loci are excluded
    if len(loc_inc) < 1:
        raise AllLociLost("Try different set of Whitelist SNPs")

    #sort locus positions by original list
    loc_possitions = [loc_l.index(locus) for locus in loc_inc]
    loc_inc_new = [x for (y,x) in sorted(zip(loc_possitions, loc_inc))]
    loc_inc = loc_inc_new
    #sorted(loc_inc, loc_l)
    return(loc_inc, mess)


#Use user-specified Blacklist
#Exclude SNPs which are part of the given locus list and which are found at the Blacklist
def set_ex_loc_list(loc_l, b_f, slp, mess):
    remove = []
    not_found = []
    black = {}
    l = 0
    for line in open(b_f):
        l += 1
        line2 = line.replace(slp, "\t").replace(" ", "\t").rstrip("\n").split("\t")

        #Test format of the black list
        if len(line2) > 2:
            raise WrongLocusListFormat("The Blacklist is in a wrong format\nSee line " + str(l) + ": " + str(line))

        pos = black.get(str(line2[0]), [])
        if len(line2) == 2:
            pos.append(str(line2[1]))
        else:
            pos.append(0)
        black[str(line2[0])] = pos

    for loc, pos_l in black.items():
        for pos in pos_l:
            loci = [lo for lo in loc_l if (str(lo.name) == loc and str(lo.position) == pos)]
            if len(loci) == 1:
                remove.append(loci[0])
            if len(loci) < 1:
                not_found.append(str(loc) + slp + str(pos))

    included = []
    for lo in loc_l:
        if lo not in remove:
            included.append(lo)

    #Write message to the log file if no locus was excluded by the Blacklist
    if len(included) == len(loc_l):
        mess.append("\nWARNING: No SNP was removed by the Blacklist")

    #Write all loci which were not part of the input file into the log file
    if len(not_found) > 0:
        mess.append("\nWARNING: Some of the Blacklist SNPs were not part of the input data:")
        mess.append([lo for lo in not_found])

    #Error if all loci are excluded
    if len(included) == 0:
        raise AllLociLost("Try different set of Blacklist SNPs")

    return(included,mess)

#Repeated update of data to ensure that thresholds remain after exclusion of Individuals/Loci
def repeated_update(ind_l, loc_l, maf, mds, mdi, f, mta, mess):
    #Run functions for the first time
    loc_l_in, maf_ex, mds_ex, both_ex, fixed, mta = filter_loci(ind_l, loc_l, maf, mds, f, mta)
    ind_l_in = filter_individuals(ind_l,loc_l_in,mdi)
    maf_l = []
    mds_l = []
    redo = False
    #Test if all loci met the filtering criteria
    for lo in loc_l_in:
        fr, gts, af, dl = lo.allele_frequency_statistics(ind_l_in)
        maf_l.append(af)
        mds_l.append(dl)
        if f == True:
            if len(gts) == 1 or len(list(fr.keys())) == 1:
                redo = True
        if mta == True:
            if len(list(fr.keys())) > 2:
                redo = True

    #Repeat until all values remain accepted
    while sorted(maf_l)[0] < maf or sorted(mds_l)[-1] > mds or redo == True:
        #Run functions
        loc_l_in2, maf_ex2, mds_ex2, both_ex2, fixed2, mta2 = filter_loci(ind_l_in, loc_l_in, maf, mds,f,mta)
        ind_l_in2 = filter_individuals(ind_l_in,loc_l_in2,mdi)
        maf_ex = maf_ex + maf_ex2
        mds_ex = mds_ex + mds_ex2
        both_ex = both_ex + both_ex2
        fixed = fixed + fixed2
        mta = mta + mta2
        loc_l_in = loc_l_in2
        ind_l_in = ind_l_in2
        maf_l = []
        mds_l = []
        redo = False
        #Test if all loci met the filtering criteria
        for lo in loc_l_in:
            fr, gts, af, dl = lo.allele_frequency_statistics(ind_l_in)
            maf_l.append(af)
            mds_l.append(dl)
            if f == True:
                if len(gts) == 1 or len(list(fr.keys())) == 1:
                    redo = True
            if mta == True:
                if len(list(fr.keys())) > 2:
                    redo = True

    #Write information to log file
    mess.append("\nSNPs excluded due to minor allele frequency: " + str(maf_ex))
    mess.append("SNPs excluded due to missing data frequency: " + str(mds_ex))
    mess.append("SNPs excluded due to minor allele and missing data frequency: " + str(both_ex))

    if f == True:
        mess.append("SNPs excluded as fixed: " + str(fixed))
    if mta == True:
        mess.append("SNPs excluded as to variable: " + str(mta))

    mess.append("Individuals excluded due to missing data frequency: " + str(len(ind_l) - len(ind_l_in)))

    return(ind_l_in, loc_l_in, mess)

#Filter multiple SNPss per locus to include only a single SNP
#Take the first position in the locus list for multiple SNPs per locus
def filter_m_snps_p_loci(loc_list, mess):
    dic_locs = {}
    loc_list_in = []
    l = 0
    for lo in loc_list:
        name = lo.name
        pos = lo.position
        found = dic_locs.get(name,[])
        found.append(pos)
        dic_locs[name] = found
        if len(found) == 1:
            loc_list_in.append(lo)
        else:
            l += 1

    #Write number of excluded loci to log file
    mess.append("\nPositions excluded due to multiple SNPs: " + str(l))

    return(loc_list_in, mess)

#Data Output

#Output data in the genepop format
#Use loci from individual list, locus_list and population assignment
#Use given locus separator and individual separator
def output_pop(ind_l, loc_l, pops,title,filename,slp,isp, mess):
    filename2 = str(filename)
    output_p = open(filename2, "w")

    #Write title line and loucs list
    output_p.write(title + "\n")
    for loc in loc_l[:-1]:
        output_p.write(str(loc.name) + slp +str(loc.position) + ",")
    loc = loc_l[-1]
    output_p.write(str(loc.name) + slp +str(loc.position))

    #Write Individual data per population
    for pop in sorted(list(pops.keys())):
        inds = pops.get(pop)
        output_p.write("\nPOP")
        for ind in inds:
            if ind in ind_l:
                output_p.write("\n" + str(ind.name) + isp)
                for loc in loc_l[:-1]:
                    output_p.write(ind.get_genotype(loc) + "\t")
                loc = loc_l[-1]
                output_p.write(ind.get_genotype(loc))

    output_p.close()

    #Write statistics to log file
    mess.append("\nGenepop output")
    mess.append("Loci exported: " + str(len(loc_l)))
    mess.append("Individuals exported: " + str(len(ind_l)))
    mess.append("Populations exported: " + str(len(list(pops.keys()))))

    return(mess)

#Output data in the fasta format
#Use loci from individual list and locus_list
#Output also the locus list with defined separator
def output_fasta(ind_l, loc_l,filename, slp, mess):
    output_fa = open(filename, "w")
    output_fa = open(filename, "a")

    #Write fasta data
    for ind in ind_l[:-1]:
        fasta = ind.fasta(loc_l)
        output_fa.write(">" + str(ind.name) + "\n" + fasta + "\n")
    fasta = ind_l[-1].fasta(loc_l)
    output_fa.write(">" + str(ind_l[-1].name) + "\n" + fasta)

    output_fa.close()

    #Write locus list information
    filename2 = str(filename) + "_loci"
    output_fa_l = open(filename2, "w")
    output_fa_l = open(filename2, "a")

    for lo in loc_l[:-1]:
        output_fa_l.write(str(lo.name) + slp + str(lo.position) + "\n")
    output_fa_l.write(str(loc_l[-1].name) + slp + str(loc_l[-1].position))

    output_fa_l.close()

    #Write fasta output statistics to log file
    mess.append("\nFasta output")
    mess.append("Loci exported: " + str(len(loc_l)))
    mess.append("Individuals exported: " + str(len(ind_l)))

    return(mess)

#Generate the lfmm format
def lfmm(ind_l, loc_l, maf, mds, mdi, f, mta, mess):
    lfmm_ref = {}
    lfmm_dic = {}
    if f != True or mta != True:
        #Filter data to ensure that they met the lfmm criteria
        mess.append("\nlfmm Filtering:")
        mess.append("Exclusion of fixed SNPs and SNPs with more than two alleles")
        ind_l_in, loc_l_in, mess = repeated_update(ind_l, loc_l, maf, mds, mdi, True, True, mess)
    else:
        #Filtering not necessary
        mess.append("\nlfmm Filtering omitted, as fixed SNPs and SNPs with more than two alleles were excluded previously")
        ind_l_in = ind_l
        loc_l_in = loc_l

    #For each locus set reference allele and assign allele state to each individual
    for lo in loc_l_in:
        ref = str(list(lo.allele_frequency_statistics(ind_l_in)[0].keys())[0])
        for ind in ind_l_in:
            a1 = str(lo.individuals.get(ind)[0])
            a2 = str(lo.individuals.get(ind)[1])
            lfmm_i = lfmm_dic.get(ind, [])
            if a1 == "00":
                lfmm_i.append(9)
            else:
                if a1 == a2:
                    if a1 == ref:
                        lfmm_i.append(0)
                    else:
                        lfmm_i.append(2)
                else:
                    lfmm_i.append(1)

            lfmm_dic[ind] = lfmm_i

    return(lfmm_dic,ind_l_in,loc_l_in,mess)

#Output data in the lfmm format
#Use loci from individual list and locus_list
##Output also the locus list with defined separator as well as a individual list
def output_lfmm(ind_l, loc_l,maf, mds, mdi, f, mta,filename,slp, mess):

    #Calculate lfmm data for given locus list and individual list
    lfmm_dic,ind_l_in,loc_l_in,mess = lfmm(ind_l, loc_l,maf, mds, mdi, f, mta,mess)

    output_lfmm_d = open(filename, "w")
    output_lfmm_d = open(filename, "a")
    output_lfmm_i = open(filename + "_ind", "w")
    output_lfmm_i = open(filename + "_ind", "a")
    output_lfmm_l = open(filename + "_loci", "w")
    output_lfmm_l = open(filename + "_loci", "a")

    #Use lfmm data to generate lfmm file
    #Output individuals included in the lfmm file
    for ind in ind_l_in[:-1]:
        output_lfmm_i.write(str(ind.name) + "\n")

        lfmm_d = lfmm_dic.get(ind)
        for d in lfmm_d[:-1]:
            output_lfmm_d.write(str(d) + " ")
        output_lfmm_d.write(str(lfmm_d[-1]) + "\n")

    output_lfmm_i.write(str(ind_l_in[-1].name))

    lfmm_d = lfmm_dic.get(ind_l_in[-1])
    for d in lfmm_d[:-1]:
        output_lfmm_d.write(str(d) + " ")
    output_lfmm_d.write(str(lfmm_d[-1]))

    #Output loci included in the lfmm file
    for lo in loc_l_in[:-1]:
        output_lfmm_l.write(str(lo.name) + slp + str(lo.position) + "\n")
    output_lfmm_l.write(str(loc_l_in[-1].name) + slp + str(loc_l_in[-1].position))

    output_lfmm_d.close()
    output_lfmm_i.close()
    output_lfmm_l.close()

    #Write lfmm output statistics to log file
    mess.append("lfmm output")
    mess.append("Loci exported: " + str(len(loc_l_in)))
    mess.append("Individuals exported: " + str(len(ind_l_in)))

    return(mess)


def output_fstat(ind_l, loc_l, pops, filename, slp, mess):
    output_fst_d = open(filename, "w")
    output_fst_d = open(filename, "a")
    output_fst_i = open(filename + "_ind", "w")
    output_fst_i = open(filename + "_ind", "a")

    output_fst_d.write(str(len(ind_l)) + "\t" + str(len(loc_l)) + "\t" + str(4) + "\t" + str(2))
    for loc in loc_l:
        output_fst_d.write("\n" + str(loc.name) + slp +str(loc.position))


    #Write Individual data per population
    i = 0
    for pop in sorted(list(pops.keys())):
        inds = pops.get(pop)
        for ind in inds:
            if ind in ind_l:
                output_fst_d.write("\n" + str(i+1))
                output_fst_i.write(ind.name + "\n")
                for loc in loc_l:
                    output_fst_d.write("\t" + ind.get_genotype(loc))
        i += 1

    output_fst_d.close()
    output_fst_i.close()

    mess.append("\nFstat output")
    mess.append("Loci exported: " + str(len(loc_l)))
    mess.append("Individuals exported: " + str(len(ind_l)))
    mess.append("Populations exported: " + str(i))

    return(mess)

def output_bayescan(ind_l, loc_l, pops, filename, slp, mess):
    output_bay_d = open(filename, "w")
    output_bay_l = open(filename + "_loci", "w")

    output_bay_d.write("[loci]=" + str(len(loc_l)))
    output_bay_d.write("\n\n[populations]=" + str(len(pops.keys())))

    for loc in loc_l:
        count_alleles, count_all = loc.alleles(ind_l)

    p = 1
    for pop in sorted(list(pops.keys())):
        inds = pops.get(pop)
        output_bay_d.write("\n\n[pop]=" + str(p))
        l = 1
        for loc in loc_l:
            count_alleles, counts_all = loc.alleles(inds)
            output_bay_d.write("\n" + str(l) + "\t" + str(counts_all) + "\t" + str(len(loc.order)) + "\t")
            for allele in loc.order:
                output_bay_d.write(str(count_alleles.get(allele,0)) + " ")
            l += 1
        p += 1

    output_bay_d.close()

    for loc in loc_l:
        output_bay_l.write(str(loc.name) + slp + str(loc.position) + "\n")

    mess.append("\nBayescan output")
    mess.append("Loci exported: " + str(len(loc_l)))
    mess.append("Individuals exported: " + str(len(ind_l)))
    mess.append("Populations exported: " + str(p - 1))

    return(mess)

#Test input name and output
def test_names(input1, gpop, lf, fas,bay, test):
    mess = []
    if input1 == "" and test == False:
        raise NoInputDefind
    elif input1 == "" and test == True:
        mess = ["The program is tested for its functionality"]
        mess = test_program(mess)
    elif input1 != "" and test == True:
        mess = ["The program is tested for its functionality"]
        mess = test_program(mess)
        if gpop == "" and lf == "" and fas =="" and fst == "" and bay == "":
            raise NoOutputDefind

    else:
        if gpop == "" and lf == "" and fas =="" and fst == "" and bay == "":
            raise NoOutputDefind
    return(mess)

#Validate program with defined dataset and parameters
def test_program(mess):
    mess2 = []
    title, pop, locs, inds, mess2 = load_data("validation/test.pop",".",mess2)

    #Test population assignment
    popn, inds2, mess2 = assign_new_pop(inds,"validation/pop_assign", mess2)
    assert_equal([ind.name for ind in inds2], ["Ind1","Ind2","Ind8"])
    assert_items_equal(list(popn.keys()), ["p1","12"])
    assert_items_equal([ind.name for ind in popn.get("p1")], ["Ind1","Ind8"])
    mess.append("Population assignment test successful")

    # Test Whitelist
    slp3 = "."
    locs2, mess2 = set_ud_loc_list(locs, "validation/white_black", ".", mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","1","2"])
    assert_equal([str(lo.position) for lo in locs2], ["1","3","1"])
    mess.append("Whitelist test successful")

    #Test Blacklist
    locs2, mess2 = set_ex_loc_list(locs, "validation/white_black", ".", mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","2","3","4","5"])
    assert_equal([str(lo.position) for lo in locs2], ["2","4","0","1","2"])
    mess.append("Blacklist test successful")

    #Test repeated update function with different settings
    inds2, locs2, mess2 = repeated_update(inds, locs, 0.2, 1, 1, False, False, mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","1","1","2","2","3"])
    assert_equal([str(lo.position) for lo in locs2], ["1","2","3","1","4","0"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0.1, 1, 1, False, False, mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","1","1","2","2","3","4","5"])
    assert_equal([str(lo.position) for lo in locs2], ["1","2","3","1","4","0","1","2"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0, 0.3, 1, False, False, mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","2","2","3","4","5"])
    assert_equal([str(lo.position) for lo in locs2], ["1","1","4","0","1","2"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0, 1, 0.3, False, False, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind1","Ind2","Ind3","Ind4","Ind6","Ind7"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0, 0.3, 0.3, False, False, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind2","Ind3","Ind4","Ind6","Ind7"])
    assert_equal([str(lo.name) for lo in locs2], ["1","2","2","3","4","5"])
    assert_equal([str(lo.position) for lo in locs2], ["1","1","4","0","1","2"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0.2, 0.3, 0.3, False, False, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind2","Ind3","Ind4","Ind5","Ind6","Ind7"])
    assert_equal([str(lo.name) for lo in locs2], ["1","2","2","3"])
    assert_equal([str(lo.position) for lo in locs2], ["1","1","4","0"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0.2, 0.5, 0.5, True, False, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind1","Ind2","Ind3","Ind4","Ind5","Ind6","Ind7"])
    assert_equal([str(lo.name) for lo in locs2], ["1","2","2"])
    assert_equal([str(lo.position) for lo in locs2], ["3","1","4"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0.2, 0.5, 0.5, False, True, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind1","Ind2","Ind3","Ind4","Ind5","Ind6","Ind7"])
    assert_equal([str(lo.name) for lo in locs2], ["1","1","3"])
    assert_equal([str(lo.position) for lo in locs2], ["1","3","0"])

    #Test filtering of only one SNP per locus
    locs2, mess2 = filter_m_snps_p_loci(locs, mess2)
    assert_equal([str(lo.name) for lo in locs2], ["1","2","3","4","5"])
    assert_equal([str(lo.position) for lo in locs2], ["1","1","0","1","2"])

    inds2, locs2, mess2 = repeated_update(inds, locs, 0.2, 0.3, 0.3, False, False, mess2)
    locs2, mess2 = filter_m_snps_p_loci(locs2, mess2)
    inds2, locs2, mess2 = repeated_update(inds2, locs2, 0.2, 0.3, 0.3, False, False, mess2)
    assert_equal([str(ind.name) for ind in inds2], ["Ind2","Ind3","Ind4","Ind6","Ind7"])
    assert_equal([str(lo.name) for lo in locs2], ["1","2","3"])
    assert_equal([str(lo.position) for lo in locs2], ["1","1","0"])
    mess.append("Filtering test successful")

    #Test fasta format
    Ind1_fasta = [ind.fasta(locs) for ind in inds if ind.name=="Ind1"][0]
    Ind8_fasta = [ind.fasta(locs) for ind in inds if ind.name=="Ind8"][0]
    assert_equal(Ind1_fasta,"NKTMNSGG")
    assert_equal(Ind8_fasta,"ANNNNNNA")
    mess.append("Fasta test successful")

    #Test lfmm format
    lfmm_dic,inds2,locs2,mess2 = lfmm(inds, locs, 0, 1, 1, False, False, mess2)
    refs = [str(list(lo.allele_frequency_statistics(inds2)[0].keys())[0]) for lo in locs2]
    ind1_l = []
    ind2_l = []
    ind3_l = []
    ind8_l = [9,9]
    if refs[0] == "04":
        ind1_l.append(0)
        ind2_l.append(0)
        ind3_l.append(2)
    else:
        ind1_l.append(2)
        ind2_l.append(2)
        ind3_l.append(0)
    if refs[1] == "03":
        ind1_l.append(0)
        ind2_l.append(2)
        ind3_l.append(0)
    else:
        ind1_l.append(2)
        ind2_l.append(0)
        ind3_l.append(2)

    if refs[2] == "03":
        ind1_l.append(0)
        ind2_l.append(0)
        ind3_l.append(0)
        ind8_l.append(2)
    else:
        ind1_l.append(2)
        ind2_l.append(2)
        ind3_l.append(2)
        ind8_l.append(0)

    Ind1 = [ind for ind in inds if ind.name=="Ind1"][0]
    Ind2 = [ind for ind in inds if ind.name=="Ind2"][0]
    Ind3 = [ind for ind in inds if ind.name=="Ind3"][0]
    Ind8 = [ind for ind in inds if ind.name=="Ind8"][0]
    assert_equal(lfmm_dic.get(Ind1), ind1_l)
    assert_equal(lfmm_dic.get(Ind2), ind2_l)
    assert_equal(lfmm_dic.get(Ind3), ind3_l)
    assert_equal(lfmm_dic.get(Ind8), ind8_l)
    mess.append("LFMM test successful")

    return(mess)

#Define command-line variables
parser = argparse.ArgumentParser(description = "Formating and filtering of genepop files", formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', '--input',
    nargs="?",
    default="",
    help="input file in genepop format")

parser.add_argument('-I', '--info',
    action="store_true",
    default=False,
    help="generates a file including input parameters and output overview\nDefault: information is written to standard output")

parser.add_argument('--maf',
    nargs="?",
    type=float,
    default=0,
    help="minor allele frequency per SNP\nDefault: 0")

parser.add_argument('--mds',
    nargs="?",
    type=float,
    default=1,
    help="maximum percentage of missing data per SNP\nDefault: 1")

parser.add_argument('--mdi',
    nargs="?",
    type=float,
    default=1,
    help="maximum percentage of missing data per Individual\nDefault: 1")

parser.add_argument('--fix',
    action="store_true",
    default=False,
    help="exclude fixed SNPs\nDefault: not excluded")

parser.add_argument('--mta',
    action="store_true",
    default=False,
    help="exclude SNPs with more than two alleles\nDefault: not excluded")

parser.add_argument('--single',
    action="store_true",
    default=False,
    help="if multiple SNPs are included per locus, use only first SNP\nDefault: all SNPs per locus used")

parser.add_argument('-p', '--pop',
    nargs="?",
    type=str,
    help="file to assign new populations\nIndividuals not assigned to any population or not included in the file are excluded from the further analysis\nFormat:\nOne line per individual\nIndividual name;tab or whitespace;population\nDefault: no assignment of new populations")

parser.add_argument('-w', '--whitelist',
    nargs="?",
    type=str,
    help="Only consider SNPs in this file\nFormat:\nOne line per SNP\nSNP name;locus position separator3;population\nDefault: no Whitelist")

parser.add_argument('-b', '--blacklist',
    nargs="?",
    type=str,
    help="Exclude all SNPs from this File from the dataset\nFormat:\nOne line per SNP\nSNP name;locus position separator3;population\nDefault: no Blacklist")

parser.add_argument('-s', '--slp1',
    nargs="?",
    type=str,
    default="_",
    help="separator among locus and position\nseparator for input file\nif slp2 is not defined also separator for output file\nDefault: '_'")

parser.add_argument('--slp2',
    nargs="?",
    type=str,
    help="separator among locus and position\nseparator for output file\nif slp3 is not defined also separator Whitelist/Blacklist")

parser.add_argument('--slp3',
    nargs="?",
    type=str,
    help="separator among locus and position\nseparator for Whitelist/Blacklist\n")

parser.add_argument('-G', '--genepop',
    nargs="?",
    default="",
    help="generate genepop output\nDefault: not generated")

parser.add_argument('-L', '--lfmm',
    nargs="?",
    default="",
    help="generate LFMM output\nDefault: not generated")

parser.add_argument('-F', '--fasta',
    nargs="?",
    default="",
    help="generate fasta output\nDefault: not generated")

parser.add_argument('-f', '--fstat',
    nargs="?",
    default="",
    help="generate fstat output\nDefault: not generated")

parser.add_argument('-B', '--bayescan',
    nargs="?",
    default="",
    help="generate bayescan output\nDefault: not generated")

parser.add_argument('--isp',
    nargs="?",
    type=str,
    default=" ,  ",
    help="separator among individual and loci in genepop format\nDefault: ' ,  '")

parser.add_argument('-c', '--code',
    nargs = 4,
    default=["01","02","03","04"],
    help="use alternative genepop code for a, c, g, t\nDefault: 01 02 03 04")

parser.add_argument('-V', '--validate',
    action="store_true",
    default=False,
    help="Validate the functionality of the program")

args = parser.parse_args()

#File names
input1 = args.input
sep_file1 = args.pop
whitelist_file1 = args.whitelist
blacklist_file1 = args.blacklist

#Locus position separators
slp1 = args.slp1
slp2 = args.slp2
slp3 = args.slp3
if slp2 == None:
    slp2 = slp1
if slp3 == None:
    slp3 = slp2
isp = args.isp

#Filter options
maf = args.maf
mds = args.mds
mdi = args.mdi
fix = args.fix
mta = args.mta
single = args.single

#Outputs
info = args.info
gpop = args.genepop
lf = args.lfmm
fas = args.fasta
fst = args.fstat
bay = args.bayescan
test = args.validate

#Fasta code
code = args.code
a = code[0]
c = code[1]
g = code[2]
t = code[3]



#Run program
try:
    #Test input/output
    message = test_names(input1, gpop, lf, fas, bay, test)

    if input1 != "":
        #Load data
        title, pop, locs, inds, message = load_data(input1,slp1,message)

        #Define population => exclude individuals
        if sep_file1 != None:
            popn, inds, message = assign_new_pop(inds,sep_file1, message)
        else:
            popn = assign_old_pop(inds)

        #Include Whitelist SNPs if option is set
        if whitelist_file1 != None:
            locs, message = set_ud_loc_list(locs, whitelist_file1, slp3, message)

        #Exclude Blacklist SNPs if option is set
        if blacklist_file1 != None:
            locs, message = set_ex_loc_list(locs, blacklist_file1, slp3, message)

        #Filter data if any filter option is set
        if maf != 0 or mds != 1 or mdi != 1 or fix == True or mta == True:
            inds, locs, message = repeated_update(inds, locs, maf, mds, mdi, fix, mta, message)
        else:
            message.append("No filter options set")

        #Filter multiple SNPs if single option is set
        if single == True:
            locs, message = filter_m_snps_p_loci(locs, message)
            if maf != 0 or mds != 1 or mdi != 1 or fix == True or mta == True:
                message.append("Repeat filtering to ensure thresholds")
                inds, locs, message = repeated_update(inds, locs, maf, mds, mdi, fix, mta, message)

        #Generate output
        #Output genepop
        if gpop != "":
            message = output_pop(inds, locs, popn, title, gpop, slp2,isp,message)
        #Output fasta
        if fas != "":
            message = output_fasta(inds, locs, fas, slp2, message)
        #Output lfmm
        if lf != "":
            message = output_lfmm(inds, locs, maf, mds, mdi, fix, mta ,lf, slp2, message)

        #Output lfmm
        if fst != "":
            message = output_fstat(inds, locs, popn, fst, slp2, message)

        if bay != "":
            message = output_bayescan(inds, locs, popn, bay, slp2, message)


    #Generate info in new file if info option is set
    if info == True:
        inf = open((output1 + ".info"), "w")
        inf = open((output1 + ".info"), "a")
        for m in message:
            inf.write(m)
            inf.write("\n")
        inf.close()
    else:
        for m in message:
            print(m)

#Print specific Errors
except IOError as ex:
    print("Could not open the file: " + ex.strerror)
    print("File: " + ex.filename)
except MissingLine as ex:
    print("Wrong input file format")
    print(ex.args[0])
except WrongNumberLoci as ex:
    print("Wrong number of loci")
    print("Please ensure separation by ',' tab or whitespace")
    print(ex.args[0])
except WrongTypeLoci as ex:
    print("Locus with wrong Character")
    print("Please ensure that only '01', '02', '03', '04' and '00' are used")
    print(ex.args[0])
except AllIndividualsLost as ex:
    print("All loci were excluded after filtering")
    print(ex.args[0])
except AllLociLost as ex:
    print("All individuals were excluded after filtering")
    print(ex.args[0])
    #for snakemake
    for outp in [lf,fas,fst,gpop]:
        out = open(outp, "w")
        out.write("No data")
        out.close

except WrongLocusListFormat as ex:
    print("Error in Whitelist/Blacklist")
    print("Please ensure that each SNP is stored to a new line")
    print("Please ensure that locus names and positions are separated by the correct spacing symbol: " + slp3 + "or tab or whitespace")
    print(ex.args[0])
except WrongPopSepFormat as ex:
    print("Error with Population assignment file")
    print("Please ensure that each individual is stored to a new line")
    print("Please ensure that individual name and population are separated by the correct spacing symbol: ',' tab or whitespace")
    print(ex.args[0])
except NoInputDefind:
    print("No name is given for the input file and the program is not set to test mode")
    print("Please give an input file or set the -V flag ")
    print("Use -h or --help to get detailed information on usage")
except NoOutputDefind:
    print("None of the output flags is activated")
    print("Please chose to either output in Genepop, Fasta, Fstat or LFMM format")
