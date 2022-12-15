"""
This Script merges annotations from different sources into one single comprehensible.

Inputs:

    --prot_species and --prot_genus, are outputs received from ProtFinder chewBBACA module

    --genbank_file , output from genebank annotatios script in Schema_refinery

    --proteome_trembl and --proteome_swiss , outputs from proteome_matcher script located in Schema_refinery

    --match_to_add , annotation of another schema, needed with --matched_schemas

    --matched_schemas , schemas matched loci, obtained from match_schemas in Schema_refinery, needed with match_to_add

    -o , path to output file

Outputs:
    Annotation TSV file with selected merged items.
"""

import os
import pandas as pd
import argparse

def merger(table_1, table_2):

    """
    merges two annotation os Locus_ID:

    Arguments:
        table_1 : pandas table
            main table

        table_2 : pandas table
            table that will be added

    Returns :
        locus_mean : pandas table
            merged table
    """

    merged = pd.merge(table_1,
                        table_2,
                        on =['Locus_ID'],
                        how ='left')

    return merged

def main(species, genus, genebank, proteome_trembl, proteome_swiss, match_to_add, matched, output_path):
    
    """
    Imports and convertions to right types, with some changes to columns names
    """

    table_species = pd.read_csv(species, delimiter="\t")

    table_species.convert_dtypes()

    if genus != '':
        table_genus = pd.read_csv(genus, delimiter="\t")

        table_genus.convert_dtypes()

    if genebank != '':
        table_genebank = pd.read_csv(genebank, delimiter="\t")

        table_genebank.convert_dtypes()

        table_genebank.rename({"origin_id":"genebank_origin_id","origin_product":"genebank_origin_product",
                            "origin_name":"genebank_origin_name","origin_bsr":"genebank_origin_bsr"},axis='columns',inplace=True)

    if proteome_trembl != '':
        table_proteome_t = pd.read_csv(proteome_trembl, delimiter="\t")

        table_proteome_t.convert_dtypes()


    if proteome_swiss != '':
        table_proteome_s = pd.read_csv(proteome_swiss, delimiter="\t")

        table_proteome_s.convert_dtypes()

    if match_to_add != '' and matched != '':
        match_add = pd.read_csv(match_to_add, delimiter="\t")
        matched_table = pd.read_csv(matched, delimiter="\t",names=['Locus_ID','Locus','BSR_schema_match'])

        match_add.convert_dtypes()
        matched_table.convert_dtypes()        

        match_add = match_add[['Locus','User_locus_name','Custom_annotatiom']]

    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    
    merged_table = table_species

    """
    Following lines merges the inputs mostly by using Locus_ID similiarities,
    choosing only the essential columns from each table.
    """

    if genus != '':

        merged_table = pd.merge(merged_table,
                                table_genus[['Locus_ID','Proteome_ID','Proteome_Product',
                                'Proteome_Gene_Name','Proteome_Species','Proteome_BSR']],
                                suffixes=("_species","_genus"),
                                on =['Locus_ID'],
                                how ='left')

    if genebank != '':

        merged_table = merger(merged_table,table_genebank)


    if proteome_trembl != '':

        merged_table = merger(merged_table,table_proteome_t)

    if proteome_swiss != '':

        merged_table = merger(merged_table,table_proteome_s)  

            
    if match_to_add != '' and matched != '':

        # merge columns so that both table to add and reference have locus_ID
        merged_match = pd.merge(match_add,matched_table, on = 'Locus', 
                                how = 'left')
        # change columns names
        merged_match.columns = ['Locus_GAS','User_locus_name_GAS',
                                'Custom_annotatiom_GAS',
                                'Locus_ID','BSR_schema_match']


        merged_table = pd.merge(merged_table,
                                merged_match,
                                on = ['Locus_ID'],
                                how = 'left')
        
    return merged_table.to_csv(os.path.join(output_path,"merged_file.tsv"),sep='\t',index=False)

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--prot_species', type=str, required=True,
                        dest='species',
                        help='protfinder output for species')

    parser.add_argument('--prot_genus', type=str, required=False,
                        dest='genus',
                        default = '',
                        help='protfinder output for genus')

    parser.add_argument('--genbank_file', type=str, required=False,
                        dest='genebank',
                        default = '',
                        help='genebank file')

    parser.add_argument('--proteome_trembl', type=str, required=False,
                        dest='proteome_trembl',
                        default = '',
                        help='proteome output trembl')

    parser.add_argument('--proteome_swiss', type=str, required=False,
                        dest='proteome_swiss',
                        default = '',
                        help='proteome output swiss')

    parser.add_argument('--match_to_add', type=str, required=False,
                        dest='match_to_add',
                        default = '',
                        help='proteome matcher file output')

    parser.add_argument('--matched_schemas', type=str, required=False,
                        dest='matched',
                        default = '',
                        help='proteome matcher file output')
    
    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='output dir')
    

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))