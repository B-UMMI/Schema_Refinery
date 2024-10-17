import os
from typing import Dict, List

try:
    from utils import (
                        sequence_functions as sf,
                        blast_functions as bf,
                        linux_functions as lf,
                        file_functions as ff,
                        alignments_functions as af,
                        statistics as stats,
                        clustering_functions as cf,
    )
except ModuleNotFoundError:
    from SchemaRefinery.utils import (
                                    sequence_functions as sf,
                                    blast_functions as bf,
                                    linux_functions as lf,
                                    file_functions as ff,
                                    alignments_functions as af,
                                    statistics as stats,
                                    clustering_functions as cf,
    )

def get_schema_files(schema_directory: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Identify all of the FASTA files in the schema directory and its 'short' subdirectory.

    Parameters
    ----------
    schema_directory : str
        Path to the directory containing schema FASTA files.

    Returns
    -------
    Tuple[Dict[str, str], Dict[str, str]]
        A tuple containing two dictionaries:
        - The first dictionary maps loci names to their file paths in the schema directory.
        - The second dictionary maps loci names to their file paths in the 'short' subdirectory.
    """
    # Identify all of the FASTA files in the schema directory
    fasta_files_dict: Dict[str, str] = {
        loci.split('.')[0]: os.path.join(schema_directory, loci)
        for loci in os.listdir(schema_directory)
        if os.path.isfile(os.path.join(schema_directory, loci)) and loci.endswith('.fasta')
    }
    
    # Identify all of the FASTA files in the 'short' subdirectory
    short_folder: str = os.path.join(schema_directory, 'short')
    fasta_files_short_dict: Dict[str, str] = {
        loci.split('.')[0].split('_')[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }
    
    return fasta_files_dict, fasta_files_short_dict

def match_schemas(query_schema, subject_schema, output_directory, bsr, translation_table, cpu, run_mode):

    # Query schema files
    query_files, query_files_short = get_schema_files(query_schema)
    # Subject schema files
    subject_files, subject_files_short = get_schema_files(subject_schema)
    # Choose what files to use for the BLAST search
    query_fastas = query_files if run_mode.split('_')[0] == 'alleles' else query_files_short
    subject_fastas = subject_files_short if run_mode.split('_')[-1] == 'rep' else subject_files
    # Create the output directory
    blast_folder = os.path.join(output_directory, 'Blast')
    ff.create_directory(blast_folder)
    query_translation_folder = os.path.join(output_directory, 'Query_Translation')
    ff.create_directory(query_translation_folder)
    
    len_query_fastas, len_subject_fasta = len(query_fastas), len(subject_fastas)
    # Process query FASTA files
    query_translation_dict: Dict[str, str] = {}
    query_ids: Dict[str, List[str]] = {}
    query_translations_paths: Dict[str, str] = {}
    for query_loci, path in query_fastas.items():
        print(f"\rTranslated query loci FASTA: {i}/{len_query_fastas}", end='', flush=True)
        i += 1
        # Get the fasta sequences for the query
        fasta_dict = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        query_ids.setdefault(query_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        query_fasta_translation = os.path.join(query_translation_folder, f"{query_loci}-translation.fasta")
        # Translate sequences and update translation dictionary
        query_translations_paths[query_loci] = query_fasta_translation
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        query_fasta_translation,
                                                        None,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        False)
        # Update the query translation dictionary
        query_translation_dict.update(trans_dict)
    
    # Process subject FASTA files
    subject_translation_dict: Dict[str, str] = {}
    subject_ids: Dict[str, List[str]] = {}
    subject_translations_paths: Dict[str, str] = {}
    master_file_path = os.path.join(blast_folder, 'master_file.fasta')
    for subject_loci, path in subject_fastas.items():
        print(f"\rTranslated subject loci FASTA: {i}/{len_subject_fasta}", end='', flush=True)
        i += 1
        # Get the fasta sequences for the query
        fasta_dict = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        subject_ids.setdefault(subject_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        subject_fasta_translation = os.path.join(query_translation_folder, f"{subject_loci}-translation.fasta")
        # Translate sequences and update translation dictionary
        subject_translations_paths[subject_loci] = subject_fasta_translation
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        subject_fasta_translation,
                                                        None,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        False)
        # Write the sequences to the master file
        with open(master_file_path) as master:
            for id_, sequence in trans_dict.items():
                master.write(f">{id_}\n{sequence}\n")

        # Update the subject translation dictionary
        subject_translation_dict.update(trans_dict)
    