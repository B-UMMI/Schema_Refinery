#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""



"""


import re
import csv
import time
import argparse

from Bio import Entrez


ids_types = ['assembly', 'biosample', 'bioproject', 'sra']

database_patterns = {'biosample': 'SAM[E|D|N][A-Z]?[0-9]+',
                     'bioproject': 'PRJ[E|D|N][A-Z][0-9]+',
                     'sra': '[E|D|S]RR[0-9]{6,}',
                     'refseq': 'GCF_[0-9]{9}.[0-9]+',
                     'genbank': 'GCA_[0-9]{9}.[0-9]+'}


def determine_id_type(identifier):
    """
    """

    match = None
    for db, pat in database_patterns.items():
        db_match = re.findall(pat, identifier)
        if db_match != []:
            return db

    return match


def get_esearch_record(identifier, database):
    """
    """

    handle = Entrez.esearch(db=database, term=identifier)
    record = Entrez.read(handle)

    return record


def get_esummary_record(identifier, database):
    """
    """

    esummary_handle = Entrez.esummary(db=database, id=identifier, report='full')
    esummary_record = Entrez.read(esummary_handle, validate=False)

    return esummary_record


def get_elink_record(identifier, fromdb, todb):
    """
    """

    elink_handle = Entrez.elink(dbfrom=fromdb, db=todb, id=identifier)
    elink_record = Entrez.read(elink_handle)

    return elink_record


def get_elink_id(elink_record):
    """
    """

    elink_id = elink_record[0]['LinkSetDb']
    if len(elink_id) > 0:
        elink_ids = elink_id[0]['Link']
        elink_ids = [i['Id'] for i in elink_ids]
    else:
        elink_ids = ''

    return elink_ids


def fetch_sra_accessions(identifiers):
    """
    """

    sra_accessions = []
    for i in identifiers:
        # Get SRA summary
        sra_record = get_esummary_record(i, 'sra')

        # get SRA identifier
        sra_accession = re.findall(database_patterns['sra'],
                                   sra_record[0]['Runs'])
        if len(sra_accession) > 0:
            sra_accessions.append(sra_accession[0])

    return sra_accessions


def fetch_assembly_accessions(identifiers):
    """
    """

    refseq_accessions = []
    genbank_accessions = []
    biosample_ids = []
    biosample_accessions = []
    for i in identifiers:
        # Get Assembly Summary
        assembly_record = get_esummary_record(i, 'assembly')
        # get RefSeq identifier
        refseq_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('RefSeq', '')
        refseq_accessions.append(refseq_accession)
        # get GenBank identifier
        genbank_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('Genbank', '')
        genbank_accessions.append(genbank_accession)

        # get Biosample accession number
        biosample_id = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleId', '')
        biosample_ids.append(biosample_id)
        biosample_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleAccn', '')
        biosample_accessions.append(biosample_accession)

    return [refseq_accessions,
            genbank_accessions,
            list(set(biosample_ids)),
            list(set(biosample_accessions))]


input_file = '/home/rfm/Desktop/test_schema_refinery/ids_ncbi_test_1_100.txt'
output_file = '/home/rfm/Desktop/test_schema_refinery/converted_ids.tsv'
email = 'rmamede@medicina.ulisboa.pt'
def main(input_file, output_file, email):

    # read identifiers
    with open(input_file, 'r') as infile:
        identifiers = infile.read().splitlines()

    # define email to make requests
    Entrez.email = email

    database_identifiers = {}
    # detect identifier type
    for i in identifiers:
        match = determine_id_type(i)
        if match is not None:
            print('{0:<} : {1}'.format(i, match.upper()))
            if match in ['refseq', 'genbank']:
                match_db = 'assembly'
            else:
                match_db = match
        else:
            pass
            print('{0:<} : {1}'.format('Could not determine database type.'))

        # get record data for identifier
        record = get_esearch_record(i, match_db)
        record_ids = record['IdList']

        # one Assembly identifier should only match one BioSample record?
        if match_db == 'assembly':
            refseq_accessions, genbank_accessions,\
            biosample_ids, biosample_accessions = fetch_assembly_accessions(record_ids)

            if len(biosample_ids) > 0:
                # link to and get SRA accessions
                elink_record = get_elink_record(biosample_ids[0], 'biosample', 'sra')
                sra_ids = get_elink_id(elink_record)
                sra_accessions = fetch_sra_accessions(sra_ids)

            database_identifiers[i] = [refseq_accessions[0],
                                       genbank_accessions[0],
                                       biosample_accessions[0],
                                       ','.join(sra_accessions)]

        elif match_db == 'biosample':
            # link and get Assembly accessions
            elink_record = get_elink_record(record_ids[0], 'biosample', 'assembly')
            assembly_ids = get_elink_id(elink_record)
            refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

            # link to and get SRA accessions
            elink_record = get_elink_record(record_ids[0], 'biosample', 'sra')
            sra_ids = get_elink_id(elink_record)
            sra_accessions = fetch_sra_accessions(sra_ids)

            database_identifiers[i] = [','.join(refseq_accessions),
                                       ','.join(genbank_accessions),
                                       i,
                                       ','.join(sra_accessions)]

        elif match_db == 'sra':
            # Get SRA Summary
            esummary_record = get_esummary_record(record_ids[0], 'sra')

            # get Biosample accession number (possible for one SRA Run accession to match multiple BioSample ids?)
            biosample_accession = esummary_record[0]['ExpXml'].split('<Biosample>')[-1].split('</Biosample>')[0]
            biosample_record = get_esearch_record(biosample_accession, 'biosample')
            biosample_id = biosample_record['IdList'][0]

            # find links to Assembly database
            # possible to get multiple RefSeq and GenBank ids
            elink_record = get_elink_record(biosample_id, 'biosample', 'assembly')
            assembly_ids = get_elink_id(elink_record)
            refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

            database_identifiers[i] = [','.join(refseq_accessions),
                                       ','.join(genbank_accessions),
                                       biosample_accession,
                                       i]

        print('RefSeq: {0:<}\n'
              'GenBank: {1:<}\n'
              'BioSample: {2:<}\n'
              'SRA: {3:<}\n'.format(*database_identifiers[i]))

    # write output table
    output_header = 'RefSeq\tGenBank\tBioSample\tSRA'
    output_lines = [output_header] + ['\t'.join(v) for k, v in database_identifiers.items()]
    output_text = '\n'.join(output_lines)
    with open(output_file, 'w') as outfile:
        outfile.write(output_text+'\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='Path to text file with NCBI database '
                             'identifiers (supports identifiers from '
                             'the Assembly, BioSample and SRA '
                             'databases).')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to output TSV file with linked '
                             'identifiers between supported databases.')

    parser.add_argument('--email', type=str,
                        required=True, dest='email',
                        help='Email to perform requests with.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
