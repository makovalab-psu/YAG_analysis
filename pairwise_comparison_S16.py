

from __future__ import print_function
import os,sys
import argparse
import re
# import math
# import numpy as np
from collections import defaultdict
from itertools import combinations
# import string
# import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

# from collections import Counter
# import errno
# import shutil


def plot(filename):
	palette = {
	'human1': 'tab:blue',
	'human2': 'tab:green',
	'chimp1': 'tab:orange',
	'chimp2': 'tab:red',
	'bonobo' : 'tab:purple',
	'gorilla' : 'tab:grey',
	'B.orang' : 'tab:pink',
	"S.orang" : 'black'
	}
	indata = pd.read_csv(filename)
	# #for fam in ['BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 'RBMY', 'TSPY', 'VCY', 'XKRY']:
	for species in ['human1','human2', 'chimp1', 'chimp2','bonobo','gorilla','B.orang','S.orang']:
		sns.countplot(data=indata, x="gene_family", hue="species", palette = palette)
		plt.xlabel("Gene family")
		plt.ylabel("Count")
		plt.title("Shared with {0}".format(species))
		plt.savefig(filename+"{0}.pdf".format(species))
		plt.clf()


def get_shared_between_species(gene_fam_shared_table, s1, s2):
	identical = 0
	for h in gene_fam_shared_table:
		if s1 in gene_fam_shared_table[h] and s2 in gene_fam_shared_table[h]:
			identical += gene_fam_shared_table[h][s1]
			identical += gene_fam_shared_table[h][s2]
	return identical

def get_nr_shared(infile, typ ='hash'):
	# initialize counts table
	
	shared_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

	for line in open(infile, 'r'):
		hash_,gene_family,species,transcript_id,seq = line.split(',')
		if typ == 'hash':
			shared_table[gene_family][hash_][species] += 1
		elif typ == 'sequence':
			shared_table[gene_family][seq][species] += 1

	for s1,s2 in combinations(['human1','human2', 'chimp1', 'chimp2','bonobo','gorilla','B.orang','S.orang'], 2):
		print(s1,s2)
		l = []
		print('BPY2,CDY,DAZ,HSFY,PRY,RBMY,TSPY,VCY,XKRY')
		for fam in ['BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 'RBMY', 'TSPY', 'VCY', 'XKRY']:
			identical = get_shared_between_species(shared_table[fam], s1,s2)
			l.append(identical)
		print(','.join([str(n) for n in l]))
	return shared_table


def main(args):
	outfile = open(args.outfile, 'w')
	outfile.write('hash,gene_family,species,transcript_id,seq\n')

	for line in open(args.infile, 'r'):
		#print(line.split())
		hash_, acc, seq = line.split()
		gene_family = acc.split('_')[0][1:]
		species = acc.split('_')[1]
		transcript_id = '_'.join(acc.split('_')[2:])
		if species == 'h' or species == 'cb':
			species = acc.split('_')[2]
			transcript_id = '_'.join(acc.split('_')[3:])

		outfile.write('{0},{1},{2},{3},{4}\n'.format(hash_, gene_family,species, transcript_id,seq.strip()))
	outfile.close()

	print('NUMBER IDENTICAL (HASH)')
	get_nr_shared(outfile.name, typ ='hash')
	print()
	print('NUMBER IDENTICAL (SEQUENCE)')
	get_nr_shared(outfile.name, typ ='sequence')

	# plot(outfile.name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('infile', type=str, help='All transcripts in space separated format from excel')
    parser.add_argument('outfile', type=str, help='Path to processed csv')

    args = parser.parse_args()

    main(args)
