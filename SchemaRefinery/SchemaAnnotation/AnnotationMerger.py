"""
This Script merges annotations from different sources into one single comprehensible.

Inputs:

    --uniprot_species and --uniprot_genus, are outputs received from ProtFinder chewBBACA module

    --genbank_file , output from genbank annotations script in Schema_refinery

    --proteome_trembl and --proteome_swiss , outputs from proteome_matcher script located in Schema_refinery

    --match_to_add , annotation of another schema, needed with --matched_schemas

    --matched_schemas , schemas matched loci, obtained from match_schemas in Schema_refinery, needed with match_to_add

    -o , path to output file

Outputs:
    Annotation TSV file with selected merged items.
"""

import os
import pandas as pd

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

def annotation_merger(uniprot_species, uniprot_genus, genebank, proteome_trembl, proteome_swiss, match_to_add, old_schema_columns, matched_schemas, output_path):
    
    """
    Imports and convertions to right types, with some changes to columns names
    """

    if uniprot_species:
        table_species = pd.read_csv(uniprot_species, delimiter="\t")

        table_species.convert_dtypes()

    if uniprot_genus:
        table_genus = pd.read_csv(uniprot_genus, delimiter="\t")

        table_genus.convert_dtypes()

    if genebank:
        table_genebank = pd.read_csv(genebank, delimiter="\t")

        table_genebank.convert_dtypes()

        table_genebank.rename({"origin_id":"genebank_origin_id","origin_product":"genebank_origin_product",
                            "origin_name":"genebank_origin_name","origin_bsr":"genebank_origin_bsr"},axis='columns',inplace=True)

    if proteome_trembl:
        table_proteome_t = pd.read_csv(proteome_trembl, delimiter="\t")

        table_proteome_t.convert_dtypes()


    if proteome_swiss:
        table_proteome_s = pd.read_csv(proteome_swiss, delimiter="\t")

        table_proteome_s.convert_dtypes()

    if match_to_add and matched_schemas:
        match_add = pd.read_csv(match_to_add, delimiter="\t")
        matched_table = pd.read_csv(matched_schemas, delimiter="\t",names=['Locus_ID','Locus','BSR_schema_match'])

        match_add.convert_dtypes()
        matched_table.convert_dtypes()        

        match_add = match_add[['Locus','User_locus_name'] + old_schema_columns]

    merged_table = pd.DataFrame({'Locus_ID': []})

    """
    Following lines merges the inputs mostly by using Locus_ID similiarities,
    choosing only the essential columns from each table.
    """

    if uniprot_species:

        merged_table = merger(table_species, merged_table)

    if uniprot_genus:

        merged_table = pd.merge(table_genus[['Locus_ID','Proteome_ID','Proteome_Product',
                                'Proteome_Gene_Name','Proteome_Species','Proteome_BSR']],
                                merged_table,
                                suffixes=("_species","_genus"),
                                on =['Locus_ID'],
                                how ='left')

    if genebank:

        merged_table = merger(table_genebank, merged_table)


    if proteome_trembl:

        merged_table = merger(table_proteome_t, merged_table)

    if proteome_swiss:

        merged_table = merger(table_proteome_s, merged_table)  

            
    if match_to_add and matched_schemas:

        # merge columns so that both table to add and reference have locus_ID
        merged_match = pd.merge(match_add, matched_table, on = 'Locus', 
                                how = 'left')
        # change columns names
        merged_match.columns = ['Locus_GAS','User_locus_name_GAS',
                                'Custom_annotatiom_GAS',
                                'Locus_ID','BSR_schema_match']


        merged_table = pd.merge(merged_match,
                                merged_table,
                                on = ['Locus_ID'],
                                how = 'left')
        
    return merged_table.to_csv(os.path.join(output_path,"merged_file.tsv"),sep='\t',index=False)
