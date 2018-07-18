#!/usr/bin/python

'''
This program adds taxon information to a fasta file of OTUs
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
ap.add_argument('--fasta',required=True,type=str,help='fasta file of OTUs')
ap.add_argument('--taxa',required=True,type=str,help='Taxa IDs of OTUs')

conf = ap.parse_args()


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

locus_dictF = defaultdict(list)

with open(conf.fasta) as f:
    fasta_lines = f.readlines()
with open(conf.taxa) as f:
    taxa_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Assign taxa to dictionary
#-----------------------------------------------------

taxa_dict = defaultdict(str)
for line in taxa_lines:
    species = ''
    genus = ''
    line = line.rstrip()
    split_line = line.split(',')
    OTU = split_line[0]
    species = split_line[-2]
    species = species.replace('.07FU', '') # This corrects outputs of the ITS database
    genus = split_line[-4]
    taxa_dict[OTU] = ",".join([OTU, genus, species])

# print taxa_dict

for line in fasta_lines:
    line = line.rstrip()
    if line.startswith('>'):
        OTU = line.replace('>', '')
        # print OTU
        new_header = ">" + taxa_dict[OTU]
        print new_header
    else:
        print line
