#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script splits records in a set of UniProt proteomes into
two Fasta files, one with records from Swiss-Prot and another
with records from TrEMBL. It saves the description in the header
of each record into a binary file create with the Pickle module.

Code documentation
------------------
"""


import os
import pickle
import argparse

from Bio import SeqIO


def main(input_directory, output_directory):
	
	if os.path.isdir(output_directory) is False:
		os.mkdir(output_directory)

	# list proteomes in input directory
	proteomes = [os.path.join(input_directory, f)
                 for f in os.listdir(input_directory)
                 if f.endswith('.fasta')]

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
	tr_file = os.path.join(output_directory, 'trembl_prots.fasta')
	with open(tr_file, 'w') as tout:
	    records = ['>{0}\n{1}'.format(k, v) for k, v in trembl.items()]
	    rectext = '\n'.join(records)
	    tout.write(rectext+'\n')

	sp_file = os.path.join(output_directory, 'sp_prots.fasta')
	with open(sp_file, 'w') as tout:
	    records = ['>{0}\n{1}'.format(k, v) for k, v in sp.items()]
	    rectext = '\n'.join(records)
	    tout.write(rectext+'\n')

	# save descriptions
	descriptions_file = os.path.join(output_directory, 'descriptions')
	with open(descriptions_file, 'wb') as dout:
		pickle.dump(descriptions, dout)


def parse_arguments():

	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-i', type=str, required=True,
						dest='input_directory',
						help='Path to directory with UniProt '
                             'proteomes in Fasta format.')

	parser.add_argument('-o', type=str, required=True,
						dest='output_directory',
						help='Path to the output directory.')

	args = parser.parse_args()

	return args


if __name__ == '__main__':

	args = parse_arguments()
	main(**vars(args))
