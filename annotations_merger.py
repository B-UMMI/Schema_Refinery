#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script merge the annotations from different sources into one single file,
when we are creating a new schema.
"""
import os
import csv
import argparse
import itertools

from collections import defaultdict

def mergeDict(dict1, dict2):
   """ Merge dictionaries and keep values of common keys in list."""
   dict3 = {**dict1, **dict2}
   for key, value in dict3.items():
       if key in dict1 and key in dict2:
               dict3[key] = [value , dict1[key]]
   return dict3


def flatten_list(list_to_flatten):
    """Flattens one level of a nested list

        Args:
            list_to_flatten (list)

        Returns:
            flattened list

        Example:

            >>> flatten_list([[[1,2],[3,4]]])
            [[1, 2], [3, 4]]

    """

    return list(itertools.chain(*list_to_flatten))


proteome_matches_files = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/proteomes/matcher_results_allele_call"

genbank_path = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/genbank_files_complete_genomes/annotations_allele_call_nuc"

uniprot_finder_path = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA"

output_dir = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA"

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
    with open(os.path.join(proteome_matches_files, "tr_blastout_corrections_auto.tsv"), "r") as trb:
        # tr_blast_results = list(csv.reader(trb, delimiter='\t'))
        tr_reader = csv.reader(trb, delimiter='\t')
        tr_dict = {rows[0]:rows[1] for rows in tr_reader}

    # switch keys with values
    # tr_blast = {v: k for k, v in tr_dict.items()}
    
    tr_blast = {} 
    for key, value in tr_dict.items(): 
        if value in tr_blast: 
            tr_blast[value].append(key) 
        else: 
            tr_blast[value]=[key] 

    # open tr proteome annotations
    with open(os.path.join(proteome_matches_files, "tr_annotations.tsv"), "r") as tra:
        tr_reader2 = csv.reader(tra, delimiter='\t')
        next(tr_reader2, None)  # skip the headers
        proteome_annotations = {rows[0]: [rows[1], rows[2], rows[3], "tr"] for rows in tr_reader2}
    
    tr_annots = {}
    
    for locus, annot in tr_dict.items():
        if annot in proteome_annotations:
            tr_annots[locus] = [
                annot,
                proteome_annotations[annot][0],
                proteome_annotations[annot][1],
                proteome_annotations[annot][2],
                proteome_annotations[annot][3]
                ]

    
    # merge proteome annotations with TrEMBL BLAST results
    # merge = mergeDict(tr_blast, proteome_annotations)


    # remove keys
    # keys_to_remove = []
    
    # for k, v in merge.items():
    #     if len(v[0]) != 3:
    #         keys_to_remove.append(k)
    
    # if len(keys_to_remove) > 0:
    #     for key in keys_to_remove:
    #       del merge[key]

    # # locus as the key
    # merge_tr = {v[1][0]: [k, v[0][0], v[0][1], v[0][2]] for k, v in merge.items()}
    
    # merge_tr = {}
    # for k, v in merge.items():
    #     # print(k)
    #     # print(v)
    #     # print(v[0])
    #     # print(v[1])
    #     merge_tr[v[1][0]] = [k, v[0][0], v[0][1], v[0][2]]


    # open sp blast results
    with open(os.path.join(proteome_matches_files, "sp_blastout_corrections_auto.tsv"), "r") as spb:
        reader_sp = csv.reader(spb, delimiter='\t')
        mydict_sp = {rows[0]:rows[1] for rows in reader_sp}

    # switch keys with values
    # sp_blast = {v: k for k, v in mydict_sp.items()}
    
    # sp_blast = dict(zip(list(mydict_sp.values()), list(mydict_sp.keys()))) 
    
    sp_blast = {} 
    for key, value in mydict_sp.items(): 
       if value in sp_blast: 
           sp_blast[value].append(key) 
       else: 
           sp_blast[value]=[key] 
           
    # open sp proteome annotations
    with open(os.path.join(proteome_matches_files, "sp_annotations.tsv"), "r") as spa:
        sp_reader3 = csv.reader(spa, delimiter='\t')
        next(sp_reader3, None)  # skip the headers
        proteome_sp_annotations = {rows[0]: [rows[1], rows[2], rows[3], "sp"] for rows in sp_reader3}
        
    sp_annots = {}
    
    for locus, annot in mydict_sp.items():
        if annot in proteome_sp_annotations:
            sp_annots[locus] = [
                annot,
                proteome_sp_annotations[annot][0],
                proteome_sp_annotations[annot][1],
                proteome_sp_annotations[annot][2],
                proteome_sp_annotations[annot][3]
            ]
            
    
            
    merge_tr_sp = mergeDict(sp_annots, tr_annots)
    
    
    rows_final = []
    for kp, vp in merge_tr_sp.items():
        if len(vp) == 5:
            if "tr" in vp:
                rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t-\t-\t-\t-'.format(kp, vp[0], vp[1], vp[2], vp[3]))
            elif "sp" in vp:
                rows_final.append('{0}\t-\t-\t-\t-\t{1}\t{2}\t{3}\t{4}'.format(kp, vp[0], vp[1], vp[2], vp[3]))
        else:
            rows_final.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(kp, vp[0][0], vp[0][1], vp[0][2], vp[0][3], vp[1][0], vp[1][1], vp[1][2], vp[1][3]))
             
   
    # sort keys to ensure the same order
    rows_final_sort = sorted(rows_final)
    
    
    
    
    # open genbank annotations
    # with open(os.path.join(genbank_path, "genbank_annotations.tsv"), "r") as gbka:
    with open("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/genbank_files_complete_genomes/genbank_annotations_all.tsv", "r") as gbka:
        gbk_reader = csv.reader(gbka, delimiter='\t')
        next(gbk_reader, None)
        genbank_annotations_dict = {rows[0]: [rows[1], rows[2], rows[3], rows[4]] for rows in gbk_reader}

    
    merge_tr_sp_gbk = mergeDict(merge_tr_sp, genbank_annotations_dict)
    
    # rows_final = []
    # for kf, vf in merge_tr_sp_gbk:
    #     if len(vf) == 4:
    #         # only tr
    #         if "tr" in vf:
    #             pass
    #         elif "sp" in vf:
    #             pass
    #     else:
    #         if
    
    rows_final_with_gbk = []
    for r in rows_final_sort:
        # print(r)
        locus_name = r.split("\t")[0]
        # print(locus_name)
        if locus_name in genbank_annotations_dict:
            rows_final_with_gbk.append(
                "{0}\t{1}\t{2}\t{3}\t{4}".format(
                    r,
                    genbank_annotations_dict[locus_name][0],
                    genbank_annotations_dict[locus_name][1],
                    genbank_annotations_dict[locus_name][2],
                    genbank_annotations_dict[locus_name][3]
                )
            )
            # r.extend(genbank_annotations_dict[locus_name])
        else:
            rows_final_with_gbk.append(
                "{0}\t{1}\t{2}\t{3}\t{4}".format(
                    r,
                    "-",
                    "-",
                    "-",
                    "-"
                )
            )
            # r.extend(["-", "-", "-", "-"])
            
    rows_final_with_gbk_sort = sorted(rows_final_with_gbk)
    
    rows_final_with_gbk_sort_dict = {r.split("\t")[0]: [r.split("\t")[1], r.split("\t")[2], r.split("\t")[3], r.split("\t")[4], r.split("\t")[5], r.split("\t")[6], r.split("\t")[7], r.split("\t")[8], r.split("\t")[9], r.split("\t")[10], r.split("\t")[11], r.split("\t")[12]] for r in rows_final_with_gbk_sort}
            
    
    #     sp_blast = {} 
    # for key, value in mydict_sp.items(): 
    #    if value in sp_blast: 
    #        sp_blast[value].append(key) 
    #    else: 
    #        sp_blast[value]=[key] 
           
           
    test = {}
    for t in rows_final_with_gbk_sort_dict.keys():
        new_key = "{0}.fasta".format(t.split("_")[0])
        # print(new_key)
        if new_key in test:
            test[new_key].append(rows_final_with_gbk_sort_dict[t])
        else:
            test[new_key] = [rows_final_with_gbk_sort_dict[t]]
        # break
    
    
    out_header = 'Locus\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME\torigin_id\torigin_product\torigin_name\torigin_bsr'
    merged_annotations = os.path.join(output_dir, 'merged_annotations_new.tsv')
    with open(merged_annotations, 'w') as trout:
        outlines = [out_header] + rows_final_with_gbk_sort
        outtext = '\n'.join(outlines)
        trout.write(outtext)

            
    # open uniprot finder results
    with open(os.path.join(uniprot_finder_path,"new_protids.tsv"), "r") as np:
        uniprot_finder_reader = csv.reader(np, delimiter='\t')
        next(uniprot_finder_reader, None)
        # new_protids_annotations_list = [[rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]] for rows in uniprot_finder_reader]
        new_protids_annotations = {rows[0]: [rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]] for rows in uniprot_finder_reader}
        
    
    final_all = []
    for npk in new_protids_annotations.keys():
        aux = []
        if npk in test:
            # aux.extend([npk, *new_protids_annotations[npk]])
            # print(aux)
            for r in test[npk]:
                aux.append([npk, *new_protids_annotations[npk], *r])
        # else:
        #     aux.extend([])
        # print(aux)
        final_all.append(aux)
        # final_all.append("\t".join(aux))
        # break
    
    # count = 0
    # final_2 = []
    # for i in final_all:
    #     if not i == []:
    #         final_2.append(i)
            
            # count += 1
    
    final_flat = sorted(flatten_list(final_all))
    
    rows_mega_final = []
    for ff in final_flat:
        rows_mega_final.append("\t".join(ff))
        
    out_header = 'Locus\tGenome\tcontig\tStart\tStop\tprotID\tname\turl\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME\torigin_id\torigin_product\torigin_name\torigin_bsr'
    merged_annotations = os.path.join(output_dir, 'merged_annotations_final.tsv')
    with open(merged_annotations, 'w') as trout:
        outlines = [out_header] + rows_mega_final
        outtext = '\n'.join(outlines)
        trout.write(outtext)
    
    
    ref_gbk_path = "/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/genbank_files_complete_genomes"
    with open(os.path.join(ref_gbk_path, "genbank_annotations_GCF_003967135.1_ASM396713v1.tsv"), "r") as gbkb:
        gbkb_reader = csv.reader(gbkb, delimiter='\t')
        next(gbkb_reader, None)
        genbank_annotations_ref_dict = {rows[0]: [rows[1], rows[2], rows[3], rows[4]] for rows in gbkb_reader}
        
    
    # get a list of the schema loci
    genbank_annotations_ref_dict_keys = list(genbank_annotations_ref_dict.keys())

    # add .fasta to loci names
    genbank_annotations_ref_dict_keys_modified = [k.replace("_1", ".fasta") for k in genbank_annotations_ref_dict_keys]

    # create a new dict with new loci names
    fasta_dict = dict(zip(genbank_annotations_ref_dict_keys_modified, list(genbank_annotations_ref_dict.values())))
    
    
    with open("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/merged_annotations_final_copy.tsv", "r") as m:
        m_reader = csv.reader(m, delimiter="\t")
        next(m_reader, None)
        merged_final_dict = {rows[0]: [rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7], rows[8], rows[9], rows[10], rows[11], rows[12], rows[13], rows[14], rows[15], rows[16], rows[17], rows[18], rows[19]] for rows in m_reader}

    
    ref_genome_all = []
    for mk in merged_final_dict.keys():
        # print(mk)
        if mk in fasta_dict:
            fd_data = fasta_dict[mk]
            # print(merged_final_dict[mk])
            # print(fd_data)
            ref_genome_all.append([mk, *merged_final_dict[mk], *fd_data])
        else:
            ref_genome_all.append([mk, *merged_final_dict[mk], "-", "-", "-", "-"])
        # break
    
    
    rows_ref_genome_final = []
    for rref in ref_genome_all:
        rows_ref_genome_final.append("\t".join(rref))
        
    out_header_ref = 'Locus\tGenome\tcontig\tStart\tStop\tprotID\tname\turl\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME\torigin_id\torigin_product\torigin_name\torigin_bsr\tReference_origin_id\tReference_origin_product\tReference_origin_name\tReference_origin_bsr'
    merged_annotations_ref = os.path.join(output_dir, 'merged_annotations_final_reference.tsv')
    with open(merged_annotations_ref, 'w') as trout:
        outlines = [out_header_ref] + rows_ref_genome_final
        outtext = '\n'.join(outlines)
        trout.write(outtext)
        
    
    with open("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/GAS_schema_annotations/sdse_gas_matches_allele_call/matches.tsv", "r") as mat:
        matches_reader = csv.reader(mat, delimiter="\t")
        gas_matches = {rows[0]: [rows[1].split("_")[0], rows[2]] for rows in matches_reader}
        
    with open("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/GAS_schema_annotations/loci_annotations_20210118_GAS.tsv", "r") as t:
        reader11 = csv.reader(t, delimiter='\t')
        next(reader11, None)
        gas_annotations = {rows[0]: [rows[6], rows[7], rows[8], rows[9], rows[10], rows[11], rows[12], rows[13], rows[14], rows[15], rows[16], rows[17], rows[18], rows[19], rows[20], rows[21], rows[22], rows[23], rows[24], rows[25], rows[28], rows[29], rows[30], rows[31], rows[32], rows[33], rows[34], rows[35], rows[36], rows[37], rows[38]] for rows in reader11}
        
    
    matches_rows = {}
    for g in gas_matches:
        try:
            ga = gas_annotations[gas_matches[g][0]]
            matches_rows["{0}.fasta".format(g)] = [gas_matches[g][0], gas_matches[g][1], ga[0], ga[1], ga[2], ga[3], ga[4], ga[5], ga[6], ga[7], ga[8], ga[9], ga[10], ga[11], ga[12], ga[13], ga[14], ga[15], ga[16], ga[17], ga[18], ga[19], ga[20], ga[21], ga[22], ga[23], ga[24], ga[25], ga[26], ga[27], ga[28], ga[29], ga[30]]
        except KeyError:
            pass
        
    with open("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/loci_annotations_20210215.tsv", "r") as loci:
        reader12 = csv.reader(loci, delimiter="\t")
        next(reader12, None)
        loci_annots = {rows[0]: [rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7], rows[8], rows[9], rows[10], rows[11], rows[12], rows[13], rows[14], rows[15], rows[16], rows[17], rows[18], rows[19], rows[20], rows[21], rows[22], rows[23], rows[24], rows[25], rows[28], rows[29], rows[30], rows[31], rows[32], rows[33], rows[34], rows[35]] for rows in reader12}

    
    final_rows = []
    
    final_keys_for_gas = sorted(list(loci_annots.keys()))
    
    for fk in final_keys_for_gas:
        try:
            matches_rows[fk]
            final_rows.append('{0}\t{1}\t{2}\t{3}\t{4}'.format(fk, matches_rows[fk][0], matches_rows[fk][1], matches_rows[fk][4], matches_rows[fk][5]))
        except KeyError:
             final_rows.append('{0}\t{1}\t{2}\t{3}\t{4}'.format(fk, '-', '-', '-', '-'))
             
    final_matches_rows_sort = sorted(final_rows)

    out_gas_matches_header = 'Locus\tGAS_locus_match\tGAS_locus_match_BSR\tGAS_User_locus_name\tGAS_Custom_Annotation'
    sdse_gas_merged_annotations = os.path.join("/home/pcerqueira/DATA/DYSGALACTIAE_SCHEMA/GAS_schema_annotations/sdse_gas_matches_allele_call", 'sdse_gas_merged_annotations.tsv')
    with open(sdse_gas_merged_annotations, 'w') as trout:
        outlines = [out_gas_matches_header] + final_matches_rows_sort
        outtext = '\n'.join(outlines)
        trout.write(outtext)


    
    # for k in matches_rows.keys():
    #     # if k in loci_annots:
    #     #     print(k)
    #     #     # print(loci_annots[k])
    #     #     print(matches_rows[k])
    #     # break
    #     try:
    #         loci_annots[k]
    #         to_add = [k, *loci_annots[k], matches_rows[k][0], matches_rows[k][1], matches_rows[k][4], matches_rows[k][5]]
    #         # to_add.extend([matches_rows[k][0], matches_rows[k][1], matches_rows[k][4], matches_rows[k][5]])
    #         # to_add = loci_annots[k].update([matches_rows[k][0], matches_rows[k][1], matches_rows[k][4], matches_rows[k][5]])
    #         final_rows.append(to_add)
    #     except KeyError:
    #         continue
            

    # final_flat_2 = flatten_list(final_2)
    
    # final_all_sort = sorted(final_all)
    
    # for rf in rows_final_with_gbk_sort_dict.keys():
    #     locus_rows_final = rf.split("_")[0]
    #     if locus_row_final
    
    
    # for nk in new_protids_annotations.keys():
    #     locus_protids = os.path.splitext(nk)[0]
    #     # print(locus_protids)
    #     if locus_protids in rows_final_with_gbk_sort_dict:
    #         final_all = [
    #             "{0}\t{1}\t{2}\t{3\t{4}\t{5}\t{6}\{7}".format(
    #                 nk,
    #                 new_protids_annotations[nk][0],
    #                 new_protids_annotations[nk][1],
    #                 new_protids_annotations[nk][2],
    #                 new_protids_annotations[nk][3],
    #                 new_protids_annotations[nk][4],
    #                 new_protids_annotations[nk][5],
    #                 new_protids_annotations[nk][6],
    #             )
    #         ] 
        

    new_protids_annotations_list2 = []
    for n in new_protids_annotations_list:
        new_protids_annotations_list2.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]))

    
    sorted_new_protids = sorted(new_protids_annotations_list2)
    
    out_protids_header = 'Fasta\tGenome\tcontig\tStart\tStop\tprotID\tname\turl'
    protids_sorted = os.path.join(output_dir, 'new_protids_sorted.tsv')
    with open(protids_sorted, 'w') as trout:
        outlines = [out_protids_header] + sorted_new_protids
        outtext = '\n'.join(outlines)
        trout.write(outtext)


    # merge proteome annotations with SwissProt BLAST results
    merge_sp = mergeDict(sp_blast, proteome_sp_annotations)

    # locus as key
    # merge_sp_switch = {v[1]: [k, v[0][0], v[0][1], v[0][2]] for k, v in merge_sp.items()}
    
    merge_sp_switch = {}
    
    for k, v in merge_sp.items():
        # print(k)
        # print(v)
        # print(v[0])
        # print(v[0][1])
        if type(v) == str:
            merge_sp_switch[v] = k
        else:
            merge_sp_switch[v[1]] = [k, v[0][0], v[0][1], v[0][2]]
            

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