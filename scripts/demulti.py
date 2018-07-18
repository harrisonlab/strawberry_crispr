#!/usr/bin/python

'''

'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import numpy as np
from sets import Set
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument('--FastqF',required=True,type=str,help='uncompressed fastq file of forward reads')
ap.add_argument('--FastqR',required=True,type=str,help='uncompressed fastq file of reverse reads')
ap.add_argument('--primer_loci',required=True,type=str,nargs='+',help='Locus names associated with primer pairs')
ap.add_argument('--primersF',required=True,type=str,nargs='+',help='Forward primer sequences')
ap.add_argument('--primersR',required=True,type=str,nargs='+',help='Reverse primer sequences')

conf = ap.parse_args()


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
# print conf.FastqF
# print conf.FastqR
recordsF = SeqIO.parse(conf.FastqF, "fastq")
recordsR = SeqIO.parse(conf.FastqR, "fastq")
loci = conf.primer_loci
primersF = conf.primersF
primersR = conf.primersR

locus_dictF = defaultdict(list)
locus_dictR = defaultdict(list)

for readF, readR in zip(recordsF, recordsR):
    hit_any_locus = False
    # print readF.seq
    # print readR.seq
    # for locus, primerF, primerR in zip(loci, primersF, primersR):
    #     # print locus
    #     # primerF = "r\"^" + primerF + "\""
    #     # primerR = "r\"^" + primerR + "\""
    #     # print primerF
    #     if any(re.search(x, str(readF.seq)) for x in [r"^" + primerF]):
    #         # print "\thitF"
    #     # print primerR
    #     if any(re.search(x, str(readR.seq)) for x in [r"^" + primerR]):
    #         # print "\thitR"
    #     if (re.search(primerF, str(readF.seq)) and re.search(primerR, str(readR.seq))):
    #         hit_this_locus = True
    #     else:
    #         hit_this_locus = False
    #     if hit_this_locus == True:
    #         locus_dictF[locus].append(readF)
    #         locus_dictR[locus].append(readR)
    #         hit_any_locus = True

    # # ---
    # # Force both reads to contain F and R primers
    # # ---
    # for locus, primerF, primerR in zip(loci, primersF, primersR):
    #     if (any(re.search(x, str(readF.seq)) for x in [r"^" + primerF]) and any(re.search(x, str(readR.seq)) for x in [r"^" + primerR])):
    #         hit_this_locus = True
    #     else:
    #         hit_this_locus = False
    #     if hit_this_locus == True:
    #         locus_dictF[locus].append(readF)
    #         locus_dictR[locus].append(readR)
    #         hit_any_locus = True
    # if hit_any_locus == False:
    #     locus_dictF['ambiguous'].append(readF)
    #     locus_dictR['ambiguous'].append(readR)

    # ---
    # Take reads if one contains a F or R primer
    # ---
    for locus, primerF, primerR in zip(loci, primersF, primersR):
        if any(re.search(x, str(readF.seq)) for x in [r"^" + primerF]):
            hit_this_locus = True
        elif any(re.search(x, str(readR.seq)) for x in [r"^" + primerR]):
            hit_this_locus = True
        else:
            hit_this_locus = False
        if hit_this_locus == True:
            locus_dictF[locus].append(readF)
            locus_dictR[locus].append(readR)
            hit_any_locus = True
    if hit_any_locus == False:
        locus_dictF['ambiguous'].append(readF)
        locus_dictR['ambiguous'].append(readR)


print loci
loci.append('ambiguous')
for locus in loci:
    print locus
    # F reads:
    out_prefix = conf.FastqF.split('.')[0]
    outfile = out_prefix + '_' + locus + '.fq'
    records = locus_dictF[locus]
    SeqIO.write(records, outfile, 'fastq')
    print(len(records))
    # R reads:
    out_prefix = conf.FastqR.split('.')[0]
    outfile = out_prefix + '_' + locus + '.fq'
    records = locus_dictR[locus]
    SeqIO.write(records, outfile, 'fastq')
    print(len(records))
