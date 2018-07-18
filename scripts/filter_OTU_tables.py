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


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--table',required=True,type=str,help='')
ap.add_argument('--threshold',required=True,type=int,help='')
conf = ap.parse_args()

threshold = conf.threshold
with open(conf.table) as f:
    table_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# Defin classes and functions
#-----------------------------------------------------

class Run_obj(object):
    # """Data associated with a single set of barcodes / primer set.
    # Attributes:
    #     OTU_obj: A counts object for OTU data
    #     species_obj: A counts object for OTU data summarised by species
    # """
    def __init__(self):
        """Return a Run_obj"""
        self.test = "hello"
        self.OTU_obj = counts_obj()
        self.species_obj = counts_obj()

    def OTUs_to_species(self):
        """"""
        OTUs = self.OTU_obj.taxa_dict.keys()
        for OTU in OTUs:
            count = self.OTU_obj.taxa_dict[OTU]
            species = OTU.split(',')[2]
            self.species_obj.add_count(species, count)

class counts_obj(object):
    # """Data associated with a single set of barcodes / primer set.
    # Attributes:
    #     OTU_obj:
    #     species_obj:
    # """
    def __init__(self):
        """Return a counts_obj"""
        self.taxa_dict = defaultdict(lambda:0)
    def add_count(self, OTU, count):
        """Add count information"""
        self.taxa_dict[OTU] += count


#-----------------------------------------------------
# Step 3
#   Import read information into an object for each run
#-----------------------------------------------------

# species_dict = defaultdict(str)

header_line = table_lines[0]
header_line = header_line.rstrip()
split_line = header_line.split("\t")
run_conditions = split_line[1:]
results_dict = defaultdict(list)
# results_dict = {}
for condition in run_conditions:
    results_dict[condition] = Run_obj()

# for condition in run_conditions:
#     print "\t".join([condition, results_dict[condition].test])

for line in table_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    OTU = split_line[0]
    for i,count in enumerate(split_line[1:]):
        condition = run_conditions[i]
        results_dict[condition].OTU_obj.add_count(OTU, int(count))

# for condition in run_conditions:
#     print "\t".join([condition, str(results_dict[condition].OTU_obj.taxa_dict)])

#-----------------------------------------------------
# Step 4
# Summarise run information by species rather than OTU
#-----------------------------------------------------

species_set = set()
for condition in run_conditions:
    results_dict[condition].OTUs_to_species()
    species_present = results_dict[condition].species_obj.taxa_dict.keys()
    for x in species_present:
        species_set.add(x)

# print species_set
out_line = ['Species']
out_line.extend(run_conditions)
print "\t".join(out_line)
for species in species_set:
    species_count = 0
    out_line = [species]
    for condition in run_conditions:
        count = results_dict[condition].species_obj.taxa_dict[species]
        count = np.subtract(count, threshold)
        count = max(0, count)
        out_line.append(str(count))
        species_count += count
    if species_count > 0:
        print "\t".join(out_line)

# for condition in run_conditions:
#     print "\t".join([condition, str(results_dict[condition].species_obj.taxa_dict)])
