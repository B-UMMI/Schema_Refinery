#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script merge the annotations from different sources into one single file,
when we are creating a new schema.
"""
import os
import csv


def mergeDict(dict1, dict2):
   """ Merge dictionaries and keep values of common keys in list."""
   dict3 = {**dict1, **dict2}
   for key, value in dict3.items():
       if key in dict1 and key in dict2:
               dict3[key] = [value , dict1[key]]
   return dict3


def main(proteome_matches_files, genbank_path, uniprot_finder_path, output_dir):
    """

    Parameters
    ----------
    proteome_matches_files : str
        Path to the directory of the output files of proteome_matcher.py.
    genbank_path : str
        Path to the directory of the genbank files of the schema input.
    uniprot_finder_path : str
        Path to the directory of the output UniprotFinder files of 
        the schema input.
    output_dir : str
        Path to output directory.

    Returns
    -------
    None.

    """

    # open tr blast results
    with open(os.path.join(proteome_matches_files, "tr_blastout.tsv"), "r") as trb:
        # tr_blast_results = list(csv.reader(trb, delimiter='\t'))
        tr_reader = csv.reader(trb, delimiter='\t')
        tr_dict = {rows[0]:rows[1] for rows in tr_reader}


    # switch keys with values
    tr_blast = {v: k for k, v in tr_dict.items()}

    # open tr proteome annotations
    with open(os.path.join(proteome_matches_files, "tr_annotations.tsv"), "r") as tra:
        tr_reader2 = csv.reader(tra, delimiter='\t')
        next(tr_reader2, None)  # skip the headers
        proteome_annotations = {rows[0]: [rows[1], rows[2], rows[3]] for rows in tr_reader2}
    
    
    # merge proteome annotations with TrEMBL BLAST results
    merge = mergeDict(tr_blast, proteome_annotations)

    # remove keys
    keys_to_remove = []
    
    for k, v in merge.items():
        if len(v[0]) != 3:
            keys_to_remove.append(k)
    
    for key in keys_to_remove:
      del merge[key]

    # locus as the key
    merge_tr = {v[1]: [k, v[0][0], v[0][1], v[0][2]] for k, v in merge.items()}

    # open sp blast results
    with open(os.path.join(proteome_matches_files, "sp_blastout.tsv"), "r") as spb:
        reader_sp = csv.reader(spb, delimiter='\t')
        mydict_sp = {rows[0]:rows[1] for rows in reader_sp}

    # switch keys with values
    sp_blast = {v: k for k, v in mydict_sp.items()}

    # open sp proteome annotations
    with open(os.path.join(proteome_matches_files, "sp_annotations.tsv"), "r") as spa:
        sp_reader3 = csv.reader(spa, delimiter='\t')
        next(sp_reader3, None)  # skip the headers
        proteome_sp_annotations = {rows[0]: [rows[1], rows[2], rows[3]] for rows in sp_reader3}

    # merge proteome annotations with SwissProt BLAST results
    merge_sp = mergeDict(sp_blast, proteome_sp_annotations)

    # locus as key
    merge_sp_switch = {v[1]: [k, v[0][0], v[0][1], v[0][2]] for k, v in merge_sp.items()}

    # merge all annotations into one dict
    merge_all = mergeDict(merge_tr, merge_sp_switch)

    # optionally write a file with proteome info only
    # rows = []
    # # out_header = 'Locus\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME'
    # for k, v in merge_all.items():
    #     if len(v) == 2:
    #         rows.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(k, v[1][0], v[1][1], v[1][2], v[1][3], v[0][0], v[0][1], v[0][2], v[0][3]))
    #     else:
    #         rows.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(k, v[0], v[1], v[2], v[3], '-', '-', '-', '-'))


    # open genbank annotations
    with open(os.path.join(genbank_path, "genbank_annotations.tsv"), "r") as gbka:
        gbk_reader = csv.reader(gbka, delimiter='\t')
        next(gbk_reader, None)
        genbank_annotations_dict = {rows[0]: [rows[1], rows[2], rows[3], rows[4]] for rows in gbk_reader}

    # merge genbank annotations with proteome annotations
    merge_all_gbk = mergeDict(genbank_annotations_dict, merge_all)

    # get a list of the schema loci
    merge_all_keys = list(merge_all_gbk.keys())

    # add .fasta to loci names
    merge_all_keys_modified = [k.replace("_1", ".fasta") for k in merge_all_keys]

    # create a new dict with new loci names
    fasta_dict = dict(zip(merge_all_keys_modified, list(merge_all_gbk.values()))) 

    # open uniprot finder results
    with open(os.path.join(uniprot_finder_path,"new_protids.tsv"), "r") as np:
        uniprot_finder_reader = csv.reader(np, delimiter='\t')
        next(uniprot_finder_reader, None)
        new_protids_annotations = {rows[0]: [rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]] for rows in uniprot_finder_reader}
    
    # merge uniprot finder annotations with the proteome+genbank annotations
    final_merge_dict = mergeDict(fasta_dict, new_protids_annotations)

    # create rows for output tsv
    # oh boy is this ugly
    # really need to fix this in the future
    rows_final = []
    for k, v in final_merge_dict.items():
        if len(v) == 7:
            rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(k, v[0], v[1], v[2], v[3], v[4], v[5], v[6], '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'))
        else:
            if len(v[1]) == 2:
                if len(v[1][0]) == 2:
                    # has everything
                    rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(k, v[0][0], v[0][1], v[0][2], v[0][3], v[0][4], v[0][5], v[0][6], v[1][0][1][0], v[1][0][1][1], v[1][0][1][2], v[1][0][1][3], v[1][0][0][0], v[1][0][0][1], v[1][0][0][2], v[1][0][0][3], v[1][1][0], v[1][1][1], v[1][1][2], v[1][1][3]))
                else:
                    # does not have sp
                    rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(k, v[0][0], v[0][1], v[0][2], v[0][3], v[0][4], v[0][5], v[0][6], v[1][0][0], v[1][0][1], v[1][0][2], v[1][0][3], '-', '-', '-', '-', v[1][1][0], v[1][1][1], v[1][1][2], v[1][1][3]))
            else:
                if len(v[1]) == 4:
                    # only trembl
                    rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(k, v[0][0], v[0][1], v[0][2], v[0][3], v[0][4], v[0][5], v[0][6], v[1][0], v[1][1], v[1][2], v[1][3], '-', '-', '-', '-', '-', '-', '-', '-'))
                else:
                    # only genbank
                    rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(k, v[0][0], v[0][1], v[0][2], v[0][3], v[0][4], v[0][5], v[0][6], '-', '-', '-', '-', '-', '-', '-', '-', v[1][0], v[1][1], v[1][2], v[1][3]))
             
    # sort keys to ensure the same order
    rows_final_sort = sorted(rows_final)

    # write the output file
    out_header = 'Locus\tGenome\tcontig\tStart\tStop\tprotID\tname\turl\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME\torigin_id\torigin_product\torigin_name\torigin_bsr'
    merged_annotations = os.path.join(output_dir, 'merged_annotations.tsv')
    with open(merged_annotations, 'w') as trout:
        outlines = [out_header] + rows_final_sort
        outtext = '\n'.join(outlines)
        trout.write(outtext)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', type=str, required=True,
                        dest='proteome_matches_files',
                        help='')

    parser.add_argument('-g', type=str, required=True,
                        dest='genbank_path',
                        help='')
    
    parser.add_argument('-u', type=str, required=True,
                        dest='uniprot_finder_path',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))