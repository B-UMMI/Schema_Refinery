#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a list of accession numbers from the NCBI
databases and searches for linked identifiers to create a TSV
file with linked identifiers between the NCBI\'s databases.

Code documentation
------------------
"""


import os
import re
import time
import argparse
import concurrent.futures
from itertools import repeat
from typing import Any, Dict, List, Optional, Tuple

from tqdm import tqdm
from Bio import Entrez


# regex expressions to identify identifier type
database_patterns: Dict[str, str] = {
    'biosample': 'SAM[E|D|N][A-Z]?[0-9]+',
    'bioproject': 'PRJ[E|D|N][A-Z][0-9]+',
    'sra': '[E|D|S]RR[0-9]{6,}',
    'refseq': 'GCF_[0-9]{9}.[0-9]+',
    'genbank': 'GCA_[0-9]{9}.[0-9]+'
}

def determine_id_type(identifier: str) -> Optional[str]:
    """
    Determine the origin database for an accession number.

    Parameters
    ----------
    identifier : str
        Accession number from one of NCBI's databases.

    Returns
    -------
    Optional[str]
        The name of the database the identifier belongs to. Returns None if it is not possible 
        to identify the identifier type.
    """
    for db, pat in database_patterns.items():
        db_match: List[str] = re.findall(pat, identifier)
        if db_match:
            return db
    return None

def get_esearch_record(identifier: str, database: str) -> Dict[str, Any]:
    """
    Run an Entrez search and return the parsed results.

    Parameters
    ----------
    identifier : str
        Accession number from one of NCBI's databases.
    database : str
        NCBI database to query.

    Returns
    -------
    Dict[str, Any]
        Multilevel data structure of Python lists and dictionaries with the query results, 
        including the primary IDs linked to the accession number.
    """
    handle = Entrez.esearch(db=database, term=identifier)
    record: Dict[str, Any] = Entrez.read(handle)
    return record

def get_esummary_record(identifier: str, database: str) -> Dict[str, Any]:
    """
    Retrieve document summaries and return parsed results.

    Parameters
    ----------
    identifier : str
        Primary ID for a NCBI record.
    database : str
        NCBI database to query.

    Returns
    -------
    Dict[str, Any]
        Multilevel data structure of Python lists and dictionaries with summary data about the record.
    """
    esummary_handle = Entrez.esummary(db=database, id=identifier, report='full')
    esummary_record: Dict[str, Any] = Entrez.read(esummary_handle, validate=False)
    return esummary_record

def get_elink_record(identifier: str, fromdb: str, todb: str) -> List[Dict[str, Any]]:
    """
    Retrieve primary IDs for links to NCBI databases.

    Parameters
    ----------
    identifier : str
        Primary ID for a NCBI record.
    fromdb : str
        Database the input identifier belongs to.
    todb : str
        Database to search for linked identifiers.

    Returns
    -------
    List[Dict[str, Any]]
        Multilevel data structure of Python lists and dictionaries with the query results, 
        including primary IDs from `todb` linked to the input identifier.
    """
    elink_handle = Entrez.elink(dbfrom=fromdb, db=todb, id=identifier)
    elink_record: List[Dict[str, Any]] = Entrez.read(elink_handle)
    return elink_record

def get_elink_id(elink_record: List[Dict[str, Any]]) -> List[str]:
    """
    Extract linked identifiers from parsed results from Entrez.elink.

    Parameters
    ----------
    elink_record : List[Dict[str, Any]]
        Multilevel data structure of Python lists and dictionaries with the query results,
        including primary IDs from `todb` linked to the input identifier.

    Returns
    -------
    List[str]
        List of primary IDs for one of NCBI's databases.
    """
    elink_id = elink_record[0]['LinkSetDb']
    if len(elink_id) > 0:
        elink_ids: List[str] = [i['Id'] for i in elink_id[0]['Link']]
    else:
        elink_ids = []

    return elink_ids

def fetch_sra_accessions(identifiers: List[str]) -> Tuple[List[str], List[str]]:
    """
    Retrieve SRA accession numbers.

    Parameters
    ----------
    identifiers : List[str]
        List of primary IDs for the SRA records.

    Returns
    -------
    Tuple[List[str], List[str]]
        sra_accessions : List[str]
            List with SRA accession numbers.
        sequencing_platforms : List[str]
            List with the sequencing platform for each SRA record.
    """
    sra_accessions: List[str] = []
    sequencing_platforms: List[str] = []
    for i in identifiers:
        # Get SRA summary
        sra_record: Dict[str, Any] = get_esummary_record(i, 'sra')

        # Get SRA identifier
        sra_accession: List[str] = re.findall(database_patterns['sra'], sra_record['Runs'])
        if sra_accession:
            sra_accessions.append(sra_accession[0])

        # Get sequencing platform
        sequencing_platform: str = sra_record['ExpXml'].split('</Platform>')[0].split('>')[-1]
        sequencing_platforms.append(sequencing_platform)

    return sra_accessions, sequencing_platforms

def fetch_assembly_accessions(identifiers: List[str]) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Retrieve Assembly accession numbers.

    Parameters
    ----------
    identifiers : List[str]
        List with primary IDs for Assembly records.

    Returns
    -------
    Tuple[List[str], List[str], List[str], List[str]]
        refseq_accessions : List[str]
            List with RefSeq accession numbers.
        genbank_accessions : List[str]
            List with GenBank accession numbers.
        biosample_ids : List[str]
            List with primary IDs for BioSample records linked to the Assembly records.
        biosample_accessions : List[str]
            List with accession numbers for BioSample records linked to the Assembly records.
    """
    refseq_accessions: List[str] = []
    genbank_accessions: List[str] = []
    biosample_ids: List[str] = []
    biosample_accessions: List[str] = []
    for i in identifiers:
        # Get Assembly Summary
        assembly_record: Dict[str, Any] = get_esummary_record(i, 'assembly')
        # Get RefSeq identifier
        refseq_accession: str = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('RefSeq', '')
        if refseq_accession:
            refseq_accessions.append(refseq_accession)
        # Get GenBank identifier
        genbank_accession: str = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('Genbank', '')
        if genbank_accession:
            genbank_accessions.append(genbank_accession)

        # Get BioSample accession number
        biosample_id: str = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleId', '')
        if biosample_id:
            biosample_ids.append(biosample_id)
        biosample_accession: str = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleAccn', '')
        if biosample_accession:
            biosample_accessions.append(biosample_accession)

    return refseq_accessions, genbank_accessions, list(set(biosample_ids)), list(set(biosample_accessions))

def multi_threading(i: str, retry: int) -> Tuple[str, List[str]]:
    """
    This function is called by concurrent.futures in order to increase the
    number of downloads per second.

    Parameters
    ----------
    i : str
        Assembly id.
    retry: int
        Number of retries per download.

    Returns
    -------
    Tuple[str, List[str]]
        i : str
            Id of iteration of multithreading.
        identifiers : List[str]
            List that contains all found identifiers for that biosample.
    """
    rtry: int = 0
    identifiers: List[str] = []
    while rtry < retry:
        try:
            match: Optional[str] = determine_id_type(i)

            if match is not None:
                if match in ['refseq', 'genbank']:
                    match_db: str = 'assembly'
                else:
                    match_db = match
            else:
                print(f'{i} : Could not determine database type.')
                match_db = "None"

            # Get record data for identifier
            record: Dict[str, Any] = get_esearch_record(i, match_db)
            record_ids: List[str] = record['IdList']

            # One Assembly identifier should only match one BioSample record?
            if match_db == 'assembly':
                refseq_accessions, genbank_accessions, biosample_ids, biosample_accessions = fetch_assembly_accessions(record_ids)

                if biosample_ids:
                    # Link to and get SRA accessions
                    elink_record: List[Dict[str, Any]] = get_elink_record(biosample_ids[0], 'biosample', 'sra')
                    sra_ids: List[str] = get_elink_id(elink_record)
                    sra_accessions, sequencing_platforms = fetch_sra_accessions(sra_ids)

                identifiers = [
                    refseq_accessions[0],
                    genbank_accessions[0],
                    biosample_accessions[0],
                    ','.join(sra_accessions),
                    ','.join(sequencing_platforms)
                ]

            elif match_db == 'biosample':
                # Link and get Assembly accessions
                elink_record = get_elink_record(record_ids[0], 'biosample', 'assembly')
                assembly_ids = get_elink_id(elink_record)
                refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

                # Link to and get SRA accessions
                elink_record = get_elink_record(record_ids[0], 'biosample', 'sra')
                sra_ids = get_elink_id(elink_record)
                sra_accessions, sequencing_platforms = fetch_sra_accessions(sra_ids)

                identifiers = [
                    ','.join(refseq_accessions),
                    ','.join(genbank_accessions),
                    i,
                    ','.join(sra_accessions),
                    ','.join(sequencing_platforms)
                ]

            elif match_db == 'sra':
                # Get SRA Summary
                esummary_record = get_esummary_record(record_ids[0], 'sra')

                # Get BioSample accession number (possible for one SRA Run accession to match multiple BioSample ids?)
                biosample_accession: str = esummary_record['ExpXml'].split('<Biosample>')[-1].split('</Biosample>')[0]
                biosample_record: Dict[str, Any] = get_esearch_record(biosample_accession, 'biosample')
                biosample_id: str = biosample_record['IdList'][0]

                # Find links to Assembly database
                # Possible to get multiple RefSeq and GenBank ids
                elink_record = get_elink_record(biosample_id, 'biosample', 'assembly')
                assembly_ids = get_elink_id(elink_record)
                refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

                identifiers = [
                    ','.join(refseq_accessions),
                    ','.join(genbank_accessions),
                    biosample_accession,
                    i
                ]
        except Exception:
            rtry += 1
            time.sleep(1)
        else:
            break

    return i, identifiers


def main(input_file: str, output_file: str, email: str, threads: int, retry: int, api_key: Optional[str]) -> None:
    """
    Main function of ncbi_linked_ids.py, this function links the input ids to
    other ids present in the databases.

    Parameters
    ----------
    input_file : str
        File containing assemblies ids.
    output_file : str
        String containing output file path.
    email : str
        Email for NCBI usage.
    threads : int
        Number of threads.
    retry : int
        Number of retries to perform in case of failure.
    api_key : Optional[str]
        String key for API of the NCBI.

    Returns
    -------
    None
        Creates a file with linked ids.
    """
    # Read identifiers
    with open(input_file, 'r', encoding='utf-8') as infile:
        identifiers: List[str] = infile.read().splitlines()

    # Define email to make requests
    Entrez.email = email

    if api_key is not None:
        Entrez.api_key = api_key

    database_identifiers: Dict[str, List[str]] = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for res in list(tqdm(executor.map(multi_threading, identifiers, repeat(retry)), total=len(identifiers))):
            database_identifiers[res[0]] = res[1]

    # Write output table
    output_header: str = 'InputID\tRefSeq\tGenBank\tBioSample\tSRA\tSequencingPlatform'
    output_lines: List[str] = [output_header]
    output_lines.extend(['\t'.join([k] + v) for k, v in database_identifiers.items()])
    output_text: str = '\n'.join(output_lines)
    with open(output_file, 'w', encoding='utf-8') as outfile:
        outfile.write(output_text + '\n')

    # Delete .dtd files
    dtd_files: List[str] = [file for file in os.listdir(os.getcwd()) if file.endswith('.dtd')]
    for file in dtd_files:
        os.remove(file)


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
