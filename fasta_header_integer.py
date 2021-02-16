#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:21:39 2021

@author: pcerqueira
"""
import os
import csv

schema_dir = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/schema_seed2"

output_dir = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/integer_headers"

reps_dir = os.path.join(schema_dir, 'short')
schema_fastas = [os.path.join(reps_dir, f) for f in os.listdir(reps_dir) if f.endswith('.fasta') is True]

# schema_fastas = [os.path.join(schema_dir, f) for f in os.listdir(schema_dir) if f.endswith('.fasta') is True]

schema_fastas_sorted = sorted(schema_fastas)

dict_list = []

fasta_dict = {}
for sf in schema_fastas_sorted:
    # read lines from FASTA file
    with open(sf, 'r') as in_fasta:
        lines = in_fasta.readlines()
    
    # change sequence headers
    integer_id = 1
    # fasta_dict = {os.path.basename(sf): {}}
    for l in range(len(lines)):
        if lines[l].startswith('>'):
            # print(lines[l])
            fasta_dict['>{0}_{1}'.format(os.path.basename(sf).split("_")[0], int(integer_id))] = lines[l]
            lines[l] = '>{0}_{1}\n'.format(os.path.basename(sf).split("_")[0], int(integer_id))
            integer_id += 1
            
    dict_list.append(fasta_dict)
    
    # write lines with changed headers to output FASTA
    with open(os.path.join(output_dir, os.path.basename(sf)), 'w') as out_fasta:
        out_fasta.writelines(lines)
    
        sequence_count = sum([1 for line in lines if '>' in line])
        print('Wrote {0} sequences to {1}'.format(sequence_count,
                                                  os.path.join(output_dir, os.path.basename(sf))))

matcher_path = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/proteomes/matcher_results_allele_call"

with open(os.path.join(matcher_path, "sp_blastout.tsv"), "r") as spb:
    reader_sp = csv.reader(spb, delimiter='\t')
    lines = [line for line in reader_sp]
    

new_lines = []
for l in lines:
    header = '>{0}'.format(l[0])
    # print(header)
    # break
    if header in fasta_dict:
        # print("yes")
        new_lines.append([fasta_dict[header].rstrip().replace(">", ""), l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], l[11], l[12]])

with open(os.path.join(matcher_path, "sp_blastout_corrections_auto.tsv"), "w") as spbo:
    writer = csv.writer(spbo, delimiter='\t')
    writer.writerows(new_lines)
    

### TREMbl

with open(os.path.join(matcher_path, "tr_blastout.tsv"), "r") as trb:
    reader_tr = csv.reader(trb, delimiter='\t')
    lines2 = [line for line in reader_tr]
    

new_lines2 = []
for l2 in lines2:
    header2 = '>{0}'.format(l2[0])
    # print(header)
    # break
    if header2 in fasta_dict:
        # print("yes")
        new_lines2.append([fasta_dict[header2].rstrip().replace(">", ""), l2[1], l2[2], l2[3], l2[4], l2[5], l2[6], l2[7], l2[8], l2[9], l2[10], l2[11], l2[12]])

with open(os.path.join(matcher_path, "tr_blastout_corrections_auto.tsv"), "w") as trbo:
    writer2 = csv.writer(trbo, delimiter='\t')
    writer2.writerows(new_lines2)
        
# for i in dict_list:
#     if "GCF-000006885-protein330_short.fasta" in i:
#         ola = i["GCF-000006885-protein330_short.fasta"]
