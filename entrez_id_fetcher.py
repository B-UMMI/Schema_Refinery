#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""



"""


import csv
import time
import argparse

from Bio import Entrez


def main(input_file, id_type, email):

    # define email to make requests
    Entrez.email = email

    # read identifiers
    with open(input_file, 'r') as infile:
        identifiers = infile.readlines()

    # decide what to do based on type of identifier

Entrez.email = 'rmamede@medicina.ulisboa.pt'

# path to file with metadata
metadata_file = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/datasets/Davies_et_al/metadata.tsv'
with open(metadata_file, 'r') as infile:
    metadata = list(csv.reader(infile, delimiter='\t'))

# get Assembly identifiers
metadata_ids = [(l[4], l)
                for l in metadata[1:]
                if l[4] != 'NA']

# get RefSeq and Genbank identifiers
assembly_ids = {}
for m in metadata_ids:
    # Assembly identifier in metadata
    identifier = m[0]
    # get UID
    handle = Entrez.esearch(db='assembly', term=identifier)
    record = Entrez.read(handle)
    # get full report for UIDs
    for uid in record['IdList']:
        # Get Assembly Summary
        esummary_handle = Entrez.esummary(db='assembly', id=uid, report='full')
        esummary_record = Entrez.read(esummary_handle, validate=False)

        # get RefSeq identifier
        refseq_accession = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('RefSeq', '')
        assembly_ids[identifier] = [refseq_accession]

        genbank_accession = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('Genbank', '')
        assembly_ids[identifier].append(genbank_accession)

        # get Biosample accession number
        biosample_id = esummary_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleId', '')
        biosample_accession = esummary_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleAccn', '')
        assembly_ids[identifier].append(biosample_accession)

        if biosample_id != '':
            # Get Biosample UID
            biosample_handle = Entrez.esearch(db='biosample', term=biosample_id)
            biosample_record = Entrez.read(biosample_handle)

            # Use Biosample UID to get SRA accession number
            biosample_uid = biosample_record['IdList'][0]
            biosample_handle2 = Entrez.esummary(db='biosample', id=biosample_uid, report='full')
            biosample_record2 = Entrez.read(biosample_handle2, validate=False)

            # Get SRA identifier
            sra_id = biosample_record2['DocumentSummarySet']['DocumentSummary'][0]['Identifiers']
            if 'SRA:' in sra_id:
                sra_id = sra_id.split('SRA: ')[-1]
            else:
                sra_id = ''
            assembly_ids[identifier].append(sra_id)

        print(identifier, assembly_ids[identifier])
        time.sleep(0.5)

# alter file with metadata
for m in metadata[1:]:
    if m[4] in assembly_ids:
        m[0] = assembly_ids[m[4]][0]
        m[1] = assembly_ids[m[4]][1]
        m[2] = assembly_ids[m[4]][2]
        m[3] = assembly_ids[m[4]][3]

new_lines = ['\t'.join(m) for m in metadata]
new_text = '\n'.join(new_lines)
metadata_out = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/datasets/Davies_et_al/metadata_out.tsv'
with open(metadata_out, 'w') as outfile:
    outfile.write(new_text+'\n')

# find other identifiers for samples that only have SRA id
# read new file
with open(metadata_out, 'r') as infile:
    latest_metadata = list(csv.reader(infile, delimiter='\t'))

# get SRA identifiers
metadata_ids = [(l[5], l)
                for l in latest_metadata[1:]
                if l[4] == 'NA' and l[5] != '-']

# get other ids
sra_ids = {}
for m in metadata_ids:
    # Assembly identifier in metadata
    identifier = m[0]
    found = False
    while found is False:
        try:
            # get UID
            handle = Entrez.esearch(db='sra', term=identifier)
            record = Entrez.read(handle)
            # get full report for UIDs
            for uid in record['IdList']:
                # Get Assembly Summary
                esummary_handle = Entrez.esummary(db='sra', id=uid, report='full')
                esummary_record = Entrez.read(esummary_handle, validate=False)
        
                # get Biosample accession number
                biosample_accession = esummary_record[0]['ExpXml'].split('<Biosample>')[-1].split('</Biosample>')[0]
        
                # Get Biosample UID
                biosample_handle = Entrez.esearch(db='biosample', term=biosample_accession)
                biosample_record = Entrez.read(biosample_handle)
                biosample_uid = biosample_record['IdList'][0]
        
                # Use Biosample UID to get Assembly accession number
                biosample_handle2 = Entrez.esummary(db='biosample', id=biosample_uid, report='full')
                biosample_record2 = Entrez.read(biosample_handle2, validate=False)
        
                # find link to Assembly database
                elink_handle = Entrez.elink(dbfrom='biosample', db='assembly', id=biosample_uid)
                elink_record = Entrez.read(elink_handle)
        
                # get Assembly UID
                assembly_uid = elink_record[0]['LinkSetDb']
                if len(assembly_uid) > 0:
                    assembly_uid = assembly_uid[0]['Link'][0]['Id']
                else:
                    assembly_uid = ''
        
                if assembly_uid != '':
                    # Get Assembly Summary
                    assembly_handle = Entrez.esummary(db='assembly', id=assembly_uid, report='full')
                    assembly_record = Entrez.read(assembly_handle, validate=False)
            
                    # get RefSeq identifier
                    refseq_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('RefSeq', '')
            
                    genbank_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('Genbank', '')
                else:
                    refseq_accession = ''
                    genbank_accession = ''
        
                sra_ids[identifier] = [refseq_accession, genbank_accession, biosample_accession, identifier]
        
                print(identifier, sra_ids[identifier])
                time.sleep(0.1)
    
            found = True
        except Exception as e:
            print(identifier, e)

# alter file with metadata
for m in latest_metadata[1:]:
    if m[5] in sra_ids:
        m[0] = sra_ids[m[5]][0]
        m[1] = sra_ids[m[5]][1]
        m[2] = sra_ids[m[5]][2]
        m[3] = sra_ids[m[5]][3]

new_lines = ['\t'.join(m) for m in latest_metadata]
new_text = '\n'.join(new_lines)
metadata_final = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/datasets/Davies_et_al/metadata_final.tsv'
with open(metadata_final, 'w') as outfile:
    outfile.write(new_text+'\n')
