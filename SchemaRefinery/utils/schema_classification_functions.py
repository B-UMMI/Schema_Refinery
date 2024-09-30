import os

try:
    from utils import (file_functions as ff,
                       iterable_functions as itf,
                       sequence_functions as sf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      iterable_functions as itf,
                                      sequence_functions as sf)

def dropped_loci_to_file(schema_loci, dropped, results_output):
    """
    Write dropped loci information to a file.

    Parameters
    ----------
    schema_loci : set
        A set of loci IDs that are part of the schema.
    dropped : list
        A list of loci IDs that have been dropped.
    results_output : str
        The directory where the dropped loci file will be saved.

    Returns
    -------
    None

    Notes
    -----
    - The function creates a file named "dropped.tsv" in the specified results_output directory.
    - The file contains three columns: ID, Reason, and Where from.
    - The Reason column is always 'Dropped_due_to_cluster_frequency_filtering'.
    - The Where from column indicates whether the loci was dropped from the schema or from possible new loci.

    Examples
    --------
    >>> schema_loci = {'locus1', 'locus2', 'locus3'}
    >>> dropped = ['locus1', 'locus4']
    >>> results_output = '/path/to/results'
    >>> dropped_loci_to_file(schema_loci, dropped, results_output)
    """
    dropped_file = os.path.join(results_output, "dropped.tsv")
    
    with open(dropped_file, 'w') as d:
        d.write("ID\tReason\tWhere from\n")
        for drop in dropped:
            if drop in schema_loci:
                dropped_from = 'from_schema'
            else:
                dropped_from = 'from_possible_new_loci'
            d.write(f"{drop}\t{'Dropped_due_to_cluster_frequency_filtering'}\t{dropped_from}\n")

def process_new_loci(schema_folder, allelecall_directory, constants, processing_mode, results_output):
    """
    Process new loci by translating sequences, counting frequencies, and preparing files for BLAST.

    Parameters
    ----------
    schema_folder : str
        Path to the folder containing schema FASTA files.
    allelecall_directory : str
        Path to the directory containing allele call results.
    constants : list
        A list of constants used for processing.
    processing_mode : str
        Mode of processing, which determines how sequences are handled, four types, reps_vs_reps
        reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.
    results_output : str
        Path to the directory where results will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - alleles (dict): Dictionary of alleles with loci IDs as keys.
        - master_file_path (str): Path to the master FASTA file.
        - translation_dict (dict): Dictionary of translated sequences.
        - frequency_in_genomes (dict): Dictionary of loci frequencies in genomes.
        - to_blast_paths (dict): Dictionary of paths to sequences to be used for BLAST.
        - all_alleles (dict): Dictionary of all alleles with loci IDs as keys.

    Notes
    -----
    - The function creates a master FASTA file containing all sequences to be used for BLAST.
    - It also translates sequences and counts their frequencies in genomes.
    - The function handles different processing modes to determine which sequences to use for BLAST.

    Examples
    --------
    >>> schema_folder = '/path/to/schema'
    >>> allelecall_directory = '/path/to/allelecall'
    >>> constants = [0, 1, 2, 3, 4, 5, 6]
    >>> processing_mode = 'alleles_rep'
    >>> results_output = '/path/to/results'
    >>> process_new_loci(schema_folder, allelecall_directory, constants, processing_mode, results_output)
    """
    # Create a dictionary of schema FASTA files
    schema = {fastafile: os.path.join(schema_folder, fastafile) for fastafile in os.listdir(schema_folder) if fastafile.endswith('.fasta')}
    
    # Create a dictionary of short schema FASTA files
    schema_short_dir = os.path.join(schema_folder, 'short')
    schema_short = {fastafile: os.path.join(schema_short_dir, fastafile) for fastafile in os.listdir(schema_short_dir) if fastafile.endswith('.fasta')}
    
    # Path to the master FASTA file
    master_file_path = os.path.join(results_output, 'master.fasta')
    
    # Create a directory for possible new loci translations
    possible_new_loci_translation_folder = os.path.join(results_output, 'schema_translation_folder')
    ff.create_directory(possible_new_loci_translation_folder)

    # Determine paths to sequences to be used for BLAST
    to_blast_paths = schema if processing_mode.split('_')[0] == 'alleles' else schema_short
    to_run_against = schema_short if processing_mode.split('_')[-1] == 'rep' else schema

    # Initialize dictionaries for alleles, translations, and frequencies
    all_alleles = {}
    alleles = {}
    translation_dict = {}
    frequency_in_genomes = {}
    temp_frequency_in_genomes = {}
    
    # Path to the CDS presence file
    cds_present = os.path.join(allelecall_directory, "temp", "2_cds_preprocess/cds_deduplication/distinct.hashtable")
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)
    
    # Process alleles to run
    for loci in to_blast_paths.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        alleles.setdefault(loci_id, {})
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():
            alleles.setdefault(loci_id, {}).update({allele_id: str(sequence)})

    # Write master file to run against
    for loci in to_run_against.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        alleles.setdefault(loci_id, {})
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():
            alleles.setdefault(loci_id, {}).update({allele_id: str(sequence)})
            # Write to master file
            write_type = 'a' if os.path.exists(master_file_path) else 'w'
            with open(master_file_path, write_type) as master_file:
                master_file.write(f">{allele_id}\n{str(sequence)}\n")

    # Count loci presence and translate all of the alleles
    for loci in schema.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        all_alleles.setdefault(loci_id, [])
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():  
            all_alleles[loci_id].append(allele_id)
            hashed_seq = sf.seq_to_hash(str(sequence))
            # If CDS sequence is present in the schema, count the number of genomes it is found in
            if hashed_seq in decoded_sequences_ids:
                # Count frequency of only presence, do not include the total CDS in the genomes
                temp_frequency_in_genomes.setdefault(loci_id, []).append(len(set(decoded_sequences_ids[hashed_seq][1:])))

        frequency_in_genomes.setdefault(loci_id, sum(temp_frequency_in_genomes[loci_id]))

        # Translate sequences and update translation dictionary
        trans_path_file = os.path.join(possible_new_loci_translation_folder, f"{loci_id}.fasta")
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        trans_path_file,
                                                        None,
                                                        constants[5],
                                                        False,
                                                        constants[6],
                                                        False)
        
        translation_dict.update(trans_dict)
                
    return alleles, master_file_path, translation_dict, frequency_in_genomes, to_blast_paths, all_alleles