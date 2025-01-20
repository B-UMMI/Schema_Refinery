import re
import time
import concurrent.futures
from itertools import repeat
from typing import Any, Dict, List, Optional, Tuple

from Bio import Entrez

try:
    from utils import (print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (print_functions as pf)

# regex expressions to identify identifier type
database_patterns: Dict[str, str] = {
    'biosample': 'SAM[E|D|N][A-Z]?[0-9]+',
    'bioproject': 'PRJ[E|D|N][A-Z][0-9]+',
    'sra': '[E|D|S]RR[0-9]{6,}',
    'refseq': 'GCF_[0-9]{9}.[0-9]+',
    'genbank': 'GCA_[0-9]{9}.[0-9]+'
}

def chunk_list(lst: List[str], max_chars: int = 8000) -> List[List[str]]:
    """
    Split a list of identifiers into chunks, ensuring that the total length of identifiers
    in each chunk does not exceed a specified maximum number of characters.

    Parameters
    ----------
    lst : list of str
        The list of identifiers to be chunked.
    max_chars : int, optional
        The maximum number of characters allowed in each chunk. Default is 8000.

    Returns
    -------
    list of list of str
        A list of chunks, where each chunk is a list of identifiers.
    """
    chunks: List[List[str]] = []  # List to store the chunks
    current_chunk: List[str] = []  # List to store the current chunk
    current_length: int = 0  # Current length of the identifiers in the current chunk

    for identifier in lst:
        item_length: int = len(identifier) + 4  # Include " OR " in the length calculation
        # Check if adding the current identifier exceeds the max_chars limit
        if current_length + item_length > max_chars:
            chunks.append(current_chunk)  # Add the current chunk to the list of chunks
            current_chunk = []  # Reset the current chunk
            current_length = 0  # Reset the current length
        current_chunk.append(identifier)  # Add the identifier to the current chunk
        current_length += item_length  # Update the current length

    if current_chunk:
        chunks.append(current_chunk)  # Add the last chunk to the list of chunks

    return chunks

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


def categorize_identifiers(identifiers: List[str]) -> Dict[str, List[str]]:
    """
    Categorize identifiers based on their type.

    Parameters
    ----------
    identifiers : list of str
        List of identifiers to categorize.

    Returns
    -------
    dict
        Dictionary with keys as identifier types and values as lists of identifiers.
    """
    categorized_identifiers: Dict[str, List[str]] = {}
    for identifier in identifiers:
        id_type = determine_id_type(identifier)
        if id_type in ['refseq', 'genbank']:
            categorized_identifiers.setdefault('assembly', []).append(identifier)
        elif id_type:
            categorized_identifiers.setdefault(id_type, []).append(identifier)
        else:
            categorized_identifiers.setdefault('unknown', []).append(identifier)

    return categorized_identifiers

def fetch_assembly_accessions(identifiers: List[str]) -> Dict[str, Dict[str, str]]:
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
    identified_ids: Dict[str, Dict[str, str]] = {}
    # Get Assembly Summary
    assemblies_records: Dict[str, Any] = get_esummary_record(identifiers, 'assembly')
    for assembly_record in assemblies_records['DocumentSummarySet']['DocumentSummary']:
        if assembly_record['RsUid'] in identifiers:
            internal_id = assembly_record['RsUid']
        else:
            internal_id = assembly_record['GbUid']
        identified_ids.setdefault(internal_id, {})

        # Get RefSeq identifier
        refseq_accession: str = assembly_record['Synonym'].get('RefSeq', '')
        identified_ids[internal_id].setdefault('RefSeq', refseq_accession)
        # Get GenBank identifier
        genbank_accession: str = assembly_record['Synonym'].get('Genbank', '')
        identified_ids[internal_id].setdefault('Genbank', genbank_accession)

        # Get BioSample accession number
        biosample_id: str = assembly_record.get('BioSampleId', '')
        identified_ids[internal_id].setdefault('BioSampleId', biosample_id)
        biosample_accession: str = assembly_record.get('BioSampleAccn', '')
        identified_ids[internal_id].setdefault('BioSampleAccn', biosample_accession)

    return identified_ids


def fetch_linked_ids(identifiers: List[str], retry: int, database: str) -> Dict[str, Dict[str, str]]:
    """
    Fetch linked identifiers from the NCBI databases.

    Parameters
    ----------
    identifiers : List[str]
        List of primary IDs for the NCBI records.
    retry : int
        Number of retries to perform in case of failure.
    database : str
        NCBI database to query.

    Returns
    -------
    Dict[str, List[str]]
        Dictionary with keys as the input identifiers and values as lists of linked identifiers.
    """
    i = 0
    while i < retry:
        try:
            if database == 'assembly':
                linked_ids = fetch_assembly_accessions(identifiers)
                # if biosample_ids:
                #     elink_records = get_elink_record(biosample_ids, 'biosample', 'sra')
                #     for elink_record in elink_records:
                #         sra_ids = get_elink_id(elink_record)
                #         sra_accession, sequencing_platform = fetch_sra_accessions(sra_ids)


            elif database == 'biosample':
                elink_records = get_elink_record(identifiers, 'biosample', 'assembly')
                assembly_ids = get_elink_id(elink_records)
                linked_ids = fetch_assembly_accessions(assembly_ids)

                # elink_records = get_elink_record(identifiers, 'biosample', 'sra')
                # sra_ids = get_elink_id(elink_records)
                # sra_accessions, sequencing_platforms = fetch_sra_accessions(sra_ids)

        except Exception:
            i += 1
            time.sleep(1)
        else:
            break
    
    return linked_ids

def get_esearch_record(identifiers: List[str], database: str, retry: int) -> Dict[str, Any]:
    """
    Run an Entrez search and return the parsed results.

    Parameters
    ----------
    identifiers : List[str]
        List of primary IDs for the NCBI records.
    database : str
        NCBI database to query.

    Returns
    -------
    Dict[str, Any]
        Multilevel data structure of Python lists and dictionaries with the query results, 
        including the primary IDs linked to the accession number.
    """
    i = 0
    record: Dict[str, Any] = {}
    while i < retry:
        try:
            search_query = ' OR '.join(identifiers)
            handle = Entrez.esearch(db=database, term=search_query, retmax=len(identifiers))
            record = Entrez.read(handle)
        except Exception:
            i += 1
            time.sleep(1)
        else:
            break
    return record

def get_esummary_record(identifiers: List[str], database: str) -> Dict[str, Any]:
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
    esummary_handle = Entrez.esummary(db=database, id=",".join(identifiers), report='full')
    esummary_record: Dict[str, Any] = Entrez.read(esummary_handle, validate=False)
    return esummary_record

def get_elink_record(identifiers: List[str], fromdb: str, todb: str) -> List[Dict[str, Any]]:
    """
    Retrieve primary IDs for links to NCBI databases.

    Parameters
    ----------
    identifier : List[str]
        Primary IDs for a NCBI record.
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
    elink_handle = Entrez.elink(dbfrom=fromdb, db=todb, id=','.join(identifiers))
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

def main(input_file: str, linked_ids_file: str, email: str, threads: int, retry: int, api_key: Optional[str]) -> None:
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

    dict_identifiers = categorize_identifiers(identifiers)

    # Define email to make requests
    Entrez.email = email

    if api_key is not None:
        Entrez.api_key = api_key

    database_identifiers: Dict[str, List[str]] = {}

    pf.print_message('Fetching internal IDs identifiers...', 'info')
    internal_dict_identifiers: Dict[str, List[List[str]]] = {}
    ids_matches: Dict[str, Dict[str, Dict[str, str]]] = {database: {} for database in dict_identifiers.keys()}
    for database, ids in dict_identifiers.items():
        split_ids: List[List[str]] = chunk_list(ids)
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for record in list(executor.map(get_esearch_record, split_ids, repeat(database), repeat(retry))):
                if record:
                    internal_dict_identifiers.setdefault(database, []).append(list(record['IdList']))

    pf.print_message('Fetching linked IDs...', 'info')
    for database, ids in internal_dict_identifiers.items():
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for linked_ids in list(executor.map(fetch_linked_ids, ids, repeat(retry), repeat(database))):
                ids_matches[database].update(linked_ids)

    pf.print_message('Writing output file...', 'info')
    # Write output table
    output_header: str = 'InputID\tRefSeq\tGenBank\tBioSample\n'
    processed_ids: List[str] = []
    with open(linked_ids_file, 'w', encoding='utf-8') as outfile:
        outfile.write(output_header)
        for database, matches in ids_matches.items():
            for internal_id, values in matches.items():
                if database == 'assembly':
                    input_id = values['RefSeq'] or values['Genbank']
                elif database == 'biosample':
                    input_id = values['BioSampleAccn']
                outfile.write(f"{input_id}\t{values.get('RefSeq', '')}\t{values.get('Genbank', '')}\t{values.get('BioSampleAccn', '')}\n")
                processed_ids.append(input_id)
            not_matched_ids = set(dict_identifiers[database]) - set(processed_ids)
            for id in not_matched_ids:
                if database == 'assembly':
                    if determine_id_type(id) == 'refseq':
                        outfile.write(f"{id}{id}\t\t\n")
                    else:
                        outfile.write(f"{id}\t{id}\t\n")
                else:
                    outfile.write(f"{id}\t\t\t{id}\n")
