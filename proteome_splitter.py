#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import pickle
import argparse

from Bio import SeqIO


def main(input_files, output_dir):
	
	if os.path.isdir(output_dir) is False:
		os.mkdir(output_dir)

	# import proteomes
	proteomes = [os.path.join(input_files, f) for f in os.listdir(input_files) if f.endswith('.fasta')]

	# divide into Swiss-Prot and TrEMBL records
	sp = {}
	trembl = {}
	descriptions = {}
	for file in proteomes:
	    for rec in SeqIO.parse(file, 'fasta'):
	        recid = rec.id
	        prot = str(rec.seq)
	        desc = rec.description
	        if recid.startswith('tr'):
	            trembl[recid] = prot
	        elif recid.startswith('sp'):
	            sp[recid] = prot
	        descriptions[recid] = desc

	# save to file
	tr_file = os.path.join(output_dir, 'trembl_prots.fasta')
	with open(tr_file, 'w') as tout:
	    records = ['>{0}\n{1}'.format(k, v) for k, v in trembl.items()]
	    rectext = '\n'.join(records)
	    tout.write(rectext)

	sp_file = os.path.join(output_dir, 'sp_prots.fasta')
	with open(sp_file, 'w') as tout:
	    records = ['>{0}\n{1}'.format(k, v) for k, v in sp.items()]
	    rectext = '\n'.join(records)
	    tout.write(rectext)

	# save descriptions
	descriptions_file = os.path.join(output_dir, 'descriptions')
	with open(descriptions_file, 'wb') as dout:
		pickle.dump(descriptions, dout)


def parse_arguments():

	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)


	parser.add_argument('-i', type=str, required=True,
						dest='input_files',
						help='Path to directory with reference proteomes.')

	parser.add_argument('-o', type=str, required=True,
						dest='output_dir',
						help='Path to output directory.')

	args = parser.parse_args()

	return args


if __name__ == '__main__':

	args = parse_arguments()

	main(**vars(args))
