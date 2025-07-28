#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import pandas as pd
import shutil
from typing import Any, List, Dict, Optional, Set, Tuple


try:
    from SchemaAnnotation import (consolidate as cs)
    from utils import (core_functions as cof,
                                        file_functions as ff,
                                        sequence_functions as sf,
                                        iterable_functions as itf,
                                        blast_functions as bf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf,
                                        constants as ct,
                                        Types as tp,
                                        print_functions as pf,
                                        logger_functions as logf,
                                        time_functions as ti,
                                        globals as gb)
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation import (consolidate as cs)
    from SchemaRefinery.utils import (core_functions as cof,
                                        file_functions as ff,
                                        sequence_functions as sf,
                                        iterable_functions as itf,
                                        blast_functions as bf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf,
                                        constants as ct,
                                        Types as tp,
                                        print_functions as pf,
                                        logger_functions as logf,
                                        time_functions as ti,
                                        globals as gb)

def create_directories(output_directory: str, run_mode: str) -> Tuple[str, str, str, str, str, Optional[str], Optional[str], str, str]:
    """
    Create necessary directories for the processing pipeline.

    Parameters
    ----------
    output_directory : str
        Path to the base output directory.
    run_mode : str
        Mode of running the process, either 'unclassified_cds' or 'schema'.

    Returns
    -------
    List[Optional[str]]
        A list of paths to the created directories.
    """
    # Create base output directory
    ff.create_directory(output_directory)
    output_d= os.path.abspath(output_directory)
    
    # Create initial processing output directory based on run mode
    initial_processing_output: str
    schema_folder: Optional[str]
    if run_mode == 'unclassified_cds':
        initial_processing_output = os.path.join(output_d, '1_CDS_processing')
        if os.path.isdir(initial_processing_output):
            shutil.rmtree(initial_processing_output)
        ff.create_directory(initial_processing_output)
        schema_folder = None
    else:
        initial_processing_output = os.path.join(output_d, '1_schema_processing')
        if os.path.isdir(initial_processing_output):
            shutil.rmtree(initial_processing_output)
        schema_folder = os.path.join(initial_processing_output, 'schema')
        ff.create_directory(schema_folder)

    # Create BLAST processing directories
    blast_output: str = os.path.join(output_d, '2_BLAST_processing')
    if os.path.isdir(blast_output):
            shutil.rmtree(blast_output)
    blastp_output: str = os.path.join(blast_output, '1_BLASTp_processing')
    blast_db: str = os.path.join(blastp_output, 'blast_db_prot')
    ff.create_directory(blast_db)

    # Create representatives BLASTn folder if run mode is 'unclassified_cds'
    representatives_blastn_folder: Optional[str]
    if run_mode == 'unclassified_cds':
        blastn_output: str = os.path.join(blast_output, '2_BLASTn_processing')
        if os.path.isdir(blastn_output):
            shutil.rmtree(blastn_output)
        representatives_blastn_folder = os.path.join(blastn_output, 'cluster_representatives_fastas_dna')
        ff.create_directory(representatives_blastn_folder)
        representatives_blastn_folder_against = os.path.join(blastn_output, 'cluster_representatives_fastas_against')
        ff.create_directory(representatives_blastn_folder_against)
    else:
        representatives_blastn_folder = None
        representatives_blastn_folder_against = None
        
    # Create results output directories
    results_output: str = os.path.join(output_d, '3_processing_results')
    if os.path.isdir(results_output):
            shutil.rmtree(results_output)
    blast_results: str = os.path.join(results_output, 'blast_results')
    ff.create_directory(blast_results)
    
    return initial_processing_output, schema_folder, blast_output, blastp_output, blast_db, representatives_blastn_folder, representatives_blastn_folder_against, results_output, blast_results


def calculate_frequency(schema_folder: str, 
                        allelecall_directory: str,
                        second_schema_folder: str, 
                        second_allelecall_directory: str,
                        temp_paths: List[str],
                        constants: List[Any], 
                        initial_processing_output: str,
                        run_mode: str, results_output: str) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, List[str]], Dict[str, str], Dict[str, str], Dict[str, int]]:
    """
    Calculates the frequency in the genomes in 3 different was depending in the run mode.

    Parameters
    ----------
    schema_folder : str
        Path to the folder containing schema FASTA files.
    allelecall_directory : str
        Path to the directory containing allele call results.
    second_schema_folder : str
        Path to the folder containing schema FASTA files of the second schema.
    second_allelecall_directory : str
        Path to the directory containing allele call results for the second schema.
    temp_paths : List[str]
        List of temporary paths.
    constants : List[Any]
        List of constants used in the process.
    initial_processing_output: str
        initial processing output directory based on run mode.
    run_mode : str
        Mode for running the module.
        schema, unclassified_cds, schema_vs_schema
    results_output : str
        Path to the directory where results will be saved.

    Returns
    -------
    Tuple
        - frequency_in_genomes (Dict[str, int]): Dictionary of loci frequencies in genomes.
        - frequency_in_genomes_second_schema (Dict[str, int]): Dictionary of loci frequencies in genomes from the second given schema.
        - frequency_cds (Dict[str, int]): Dictionary of cds frequencies in genomes.
        - cds_presence_in_genomes (Dict[str, List[str]]): Dictionary of the presence of the cds in the genomes.
        - all_nucleotide_sequences (Dict[str, str]): Dictionary with CDS IDs as keys and sequences as values.
        - dropped_alleles (Dict[str, str]): Dropped CDS dictionary.
        - cds_size (Dict[str, int]): Dictionary of the size of the cds.
    """

    # Initialize the output dictionary
    frequency_in_genomes: Dict[str, int] = {}
    frequency_in_genomes_second_schema: Dict[str, int] = {}
    cds_present: str = ""
    cds_presence_in_genomes: Dict[str, List[str]] = {}
    all_nucleotide_sequences: Dict[str, str] = {}
    dropped_alleles: Dict[str, str] = {}
    frequency_cds: Dict[str, int] = {}
    cds_size: Dict[str, int] = {}

    if run_mode == 'schema':

        pf.print_message("Calculate frequencies in genomes for the schema run_mode.", 'info')
        schema = {fastafile: os.path.join(schema_folder, fastafile) for fastafile in os.listdir(schema_folder) if fastafile.endswith('.fasta')}
        cds_present = os.path.join(allelecall_directory, "results_alleles.tsv")
        df = pd.read_csv(cds_present, sep = '\t', dtype = object)
        allele_columns = df.columns[1:]

        pf.print_message('Calculating frequencies of each locus in each genome...', 'info')
        for loci, loci_path in schema.items():
            loci_id = ff.file_basename(loci).split('.')[0]
            # For each locus count the frequency (don't count LNF, ASM or ALM)
            matching_cols = [col for col in allele_columns if loci_id in col]
            # Change the values in the dataframe into 0 (LNF, ASM, ALM) or 1 (the locus has found seen in the genome)
            presence_mask = df[matching_cols].applymap(lambda x: 0 if str(x) == 'LNF' or str(x) == 'ASM' or str(x) == 'ALM' else 1)
            # Count the frequency of each locus in all genome
            genome_presence = presence_mask.any(axis=1)
            frequency_in_genomes[loci_id] = genome_presence.sum()

    if run_mode == 'schema_vs_schema':

        pf.print_message("Calculate frequencies in genomes of the first schema for the schema_vs_schema run_mode.")
        schema = {fastafile: os.path.join(schema_folder, fastafile) for fastafile in os.listdir(schema_folder) if fastafile.endswith('.fasta')}
        cds_present = os.path.join(allelecall_directory, "results_alleles.tsv")
        df = pd.read_csv(cds_present, sep = '\t', dtype = object)
        allele_columns = df.columns[1:]

        pf.print_message('Calculating frequencies of each locus in each genome...', 'info')
        for loci, loci_path in schema.items():
            loci_id = ff.file_basename(loci).split('.')[0]
            # For each locus count the frequency (don't count LNF, ASM or ALM)
            matching_cols = [col for col in allele_columns if loci_id in col]
            # Change the values in the dataframe into 0 (LNF, ASM, ALM) or 1 (the locus has found seen in the genome)
            presence_mask = df[matching_cols].applymap(lambda x: 0 if str(x) == 'LNF' or str(x) == 'ASM' or str(x) == 'ALM' else 1)
            # Count the frequency of each locus in all genome
            genome_presence = presence_mask.any(axis=1)
            frequency_in_genomes[loci_id] = genome_presence.sum()
        
        pf.print_message("Calculate frequencies in genomes of the second schema for the schema_vs_schema run_mode.")
        second_schema = {fastafile: os.path.join(second_schema_folder, fastafile) for fastafile in os.listdir(second_schema_folder) if fastafile.endswith('.fasta')}
        cds_present_ss = os.path.join(second_allelecall_directory, "results_alleles.tsv")
        df = pd.read_csv(cds_present_ss, sep = '\t', dtype = object)
        allele_columns_ss = df.columns[1:]

        pf.print_message('Calculating frequencies of each locus in each genome...', 'info')
        for loci, loci_path in second_schema.items():
            loci_id = ff.file_basename(loci).split('.')[0]
            # For each locus count the frequency (don't count LNF, ASM or ALM)
            matching_cols = [col for col in allele_columns_ss if loci_id in col]
            # Change the values in the dataframe into 0 (LNF, ASM, ALM) or 1 (the locus has found seen in the genome)
            presence_mask = df[matching_cols].applymap(lambda x: 0 if str(x) == 'LNF' or str(x) == 'ASM' or str(x) == 'ALM' else 1)
            # Count the frequency of each locus in all genome
            genome_presence_ss = presence_mask.any(axis=1)
            frequency_in_genomes_second_schema[loci_id] = genome_presence_ss.sum()
        
    if run_mode == 'unclassified_cds':

        pf.print_message("Calculate frequencies in genomes for the unclassified_cds run_mode.")
        temp_folder: str = temp_paths[0]
        file_path_cds: str = temp_paths[1]

        # Verify if the dataset is small, if it is, keep minimum genomes in which
        # specific CDS cluster is present to 5 if not to 1% of the dataset size.
        if not constants[2]:
            ccf.set_minimum_genomes_threshold(temp_folder, constants)

        pf.print_message("Identifying CDS not present in the schema...", "info")
        # Get dict with CDS ids as key and sequence as values.
        all_nucleotide_sequences = sf.fetch_fasta_dict(file_path_cds, True)
        
        # Make IDS universally usable
        for key, value in list(all_nucleotide_sequences.items()):
            all_nucleotide_sequences[itf.replace_by_regex(key, '_', '-')] = all_nucleotide_sequences.pop(key)

        pf.print_message("Filtering missing CDS in the schema...", "info")
        cds_size, all_nucleotide_sequences, dropped_alleles = ccf.filter_cds_by_size(all_nucleotide_sequences, constants[5])

        # Count the number of CDS not present in the schema and write CDS sequence
        # into a FASTA file.

        pf.print_message("Identifying CDS present in the schema and counting frequency of missing CDSs in the genomes...", "info")
        cds_present, frequency_cds, cds_presence_in_genomes = ccf.process_cds_not_present(initial_processing_output,
                                                                                          temp_folder,
                                                                                          all_nucleotide_sequences)

    
    return frequency_in_genomes, frequency_in_genomes_second_schema, frequency_cds, cds_presence_in_genomes, all_nucleotide_sequences, dropped_alleles, cds_size




def identify_spurious_genes(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
                            annotation_paths: List[str], constants: List[Any], temp_paths: List[str],
                            run_mode: str, processing_mode: str, cpu: int, bsr: float, translation_table: int, no_cleanup: bool) -> None:
    """
    Identify spurious genes in the given schema.

    Parameters
    ----------
    schema_directory : List[str]
        Path to the schema directory.
    output_directory : str
        Path to the output directory.
    allelecall_directory : List[str]
        Path to the allele call directory.
    annotation_paths : List[str]
        Paths for the files with the annotations.
    constants : List[Any]
        List of constants used in the process.
    temp_paths : List[str]
        List of temporary paths.
    run_mode : str
        Mode of running the process.
    processing_mode : str
        Mode of processing.
    cpu : int
        Number of CPUs to use.
    bsr : int
        Value of the threshold for the BSR.
    translation_table : int
        Number of the translation table to be used.
    no_cleanup : bool
        Flag to indicate whether to clean up temporary files.

    Returns
    -------
    None
    """
    
    output_d= os.path.abspath(output_directory)

    # Create directories structure.
    (initial_processing_output,
     schema_folder,
     blast_output,
     blastn_output,
     blast_db,
     representatives_blastn_folder,
     representatives_blastn_folder_against,
     results_output,
     blast_results) = create_directories(output_d, run_mode)
    
    # Calculate the frequencies depending of the run mode
    if run_mode == 'unclassified_cds':
        
        # Calculate the frequency of the cds in genomes, cds presence, nucleotide sequences and dropped alleles
        _, _, frequency_cds, cds_presence_in_genomes, all_nucleotide_sequences, dropped_alleles, cds_size = calculate_frequency(None, None, None, None, 
                                                                                                                        temp_paths, 
                                                                                                                        constants, 
                                                                                                                        initial_processing_output, 
                                                                                                                        run_mode, results_output)

        pf.print_message("Translating and deduplicating CDS...", "info")
        all_translation_dict: Dict[str, str]
        protein_hashes: Dict[str, List[str]]
        cds_translation_size: Dict[str, int]
        # Deduplicate protein sequences to cluster sequences without problems of similar
        # having same proteins and make it quicker
        all_translation_dict, protein_hashes, cds_translation_size = ccf.translate_and_deduplicate_cds(
                                                                                                    all_nucleotide_sequences,
                                                                                                    initial_processing_output,
                                                                                                    constants
                                                                                                )

        pf.print_message("Extracting minimizers for the translated sequences and clustering...", "info")
        # Remove CDS that did not pass filtering criteria
        all_translation_dict = ccf.remove_dropped_cds(all_translation_dict, dropped_alleles, protein_hashes)
        # Sort proteins by size
        all_translation_dict = ccf.sort_by_protein_size(all_translation_dict)
        # Set types
        all_alleles: Dict[str, List[str]]
        reps_sequences: Dict[str, str]
        reps_groups: Dict[str, List[str]]
        prot_len_dict: Dict[str, int]
        # Run clustering by minimizers
        all_alleles, reps_sequences, reps_groups, prot_len_dict = ccf.cluster_by_minimizers(all_translation_dict, constants)
        # Reformat all_alleles
        all_alleles = ccf.reformat_clusters(all_alleles, protein_hashes)

        # Add again the deduplicated sequences (Since we run dna alleles, further down the pipeline
        # we would require the protein for each DNA even though they were deduplicated)
        for hash, elements in protein_hashes.items():
            deduplicated_protein = all_translation_dict[elements[0]]
            prot_len = len(deduplicated_protein)
            for element in elements[1:]:
                all_translation_dict.setdefault(element, deduplicated_protein)
                prot_len_dict.setdefault(element, prot_len)

        # Calculate the total number of clusters
        total_number_clusters: int = len(all_alleles)
        pf.print_message(f"{len(all_translation_dict)} unique proteins have been clustered into {total_number_clusters} clusters.", "info")
        
        # Calculate the number of singleton clusters
        singleton_clusters: int = len([cluster for cluster in all_alleles.values() if len(cluster) == 1])
        pf.print_message(f"Out of those clusters, {singleton_clusters} are singletons.", "info")

        # Calculate the number of clusters with more than one CDS
        multi_cds_clusters: int = total_number_clusters - singleton_clusters
        pf.print_message(f"Out of those clusters, {multi_cds_clusters} have more than one CDS.", "info")
        
        pf.print_message("\nFiltering clusters...", "info")
        # Calculate the frequency of each cluster in the genomes.
        frequency_in_genomes: Dict[str, int] = {
            rep: sum(frequency_cds[entry] for entry in cluster_members)
            for rep, cluster_members in all_alleles.items()
        }
        # Add reason for filtering out CDS.
        dropped_alleles.update({
            cds_id: 'Dropped_due_to_cluster_frequency_filtering'
            for cds_id in itf.flatten_list([
                all_alleles[rep] for rep in all_alleles if frequency_in_genomes[rep] < constants[2]
            ])
        })
        # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
        filtered_alleles: Dict[str, List[str]] = {
            rep: cluster_members for rep, cluster_members in all_alleles.items()
            if frequency_in_genomes[rep] >= constants[2]
        }

        pf.print_message(f"After filtering by CDS frequency in the genomes (>= {constants[2]}),"
                    f" out of {total_number_clusters} clusters, {len(filtered_alleles)} remained.", "info")

        # Update all_alleles with the filtered results
        all_alleles = filtered_alleles
        # Filter also the all_translation_dict dict
        all_translation_dict = {key: value for key, value in all_translation_dict.items() if key in itf.flatten_list(all_alleles.values())}
        pf.print_message("Retrieving kmers similarity and coverage between representatives...", "info")
        reps_translation_dict: Dict[str, str] = ccf.get_representative_translation_dict(all_translation_dict, all_alleles)
        # Choose which translation dict to use as representative
        if processing_mode.split('_')[0] == 'alleles':
            trans_dict = all_translation_dict
        else:
            trans_dict = reps_translation_dict
        # Calculate kmers similarity
        reps_kmers_sim: Dict[str, float] = ccf.calculate_kmers_similarity(trans_dict, reps_groups, prot_len_dict)

        # Remove filtered out elements from reps_kmers_sim
        reps_kmers_sim = {key: value for key, value in reps_kmers_sim.items() if key in itf.flatten_list(all_alleles.values())}
        pf.print_message("Replacing CDSs IDs with the cluster representative ID...", "info")
        cds_original_ids: Dict[str, List[str]] = ccf.replace_ids_in_clusters(all_alleles,
                                                    frequency_cds,
                                                    dropped_alleles,
                                                    all_nucleotide_sequences,
                                                    prot_len_dict,
                                                    all_translation_dict,
                                                    protein_hashes,
                                                    cds_presence_in_genomes,
                                                    reps_kmers_sim)

        to_blast_paths: Dict[str, str]
        master_file_path: str
        trans_paths: Dict[str, str]
        translation_files_paths = os.path.join(initial_processing_output, "cluster_rep_translations")
        ff.create_directory(translation_files_paths)
        seqid_files_paths = os.path.join(initial_processing_output, "cluster_rep_seqid")
        ff.create_directory(seqid_files_paths)
        # Prepare the schema to run blasts (similar to prepare_loci)
        pf.print_message("Prepare the CDS for runing blasts.", "info")
        (to_blast_paths,
        to_run_against, 
        master_file_path, 
        trans_paths, 
        new_max_hits, 
        seqid_file_dict) = ccf.prepare_files_to_blast(translation_files_paths,
                                                        seqid_files_paths,
                                                        representatives_blastn_folder,
                                                        representatives_blastn_folder_against,
                                                        all_alleles,
                                                        all_nucleotide_sequences,
                                                        processing_mode,
                                                        constants)


    if run_mode == 'schema':

        ff.copy_folder(schema_directory[0], schema_folder)

        # Calculate the frequency of the cds in genomes
        frequency_in_genomes, _, _, _, _, dropped_alleles, _ = calculate_frequency(schema_folder, allelecall_directory[0], None, None, 
                                                None, 
                                                None, 
                                                None, 
                                                run_mode, None)
        
        pf.print_message('Prepare loci files for Blast and count frequencies.', 'info')
        (all_nucleotide_sequences,
        master_file_path,
        trans_paths,
        to_blast_paths,
        all_alleles,
        group_reps_ids,
        group_alleles_ids,
        to_run_against,
        new_max_hits,
        seqid_file_dict) = cof.prepare_loci(schema_folder,
                                            constants,
                                            processing_mode,
                                            initial_processing_output)
        
    if run_mode == 'schema_vs_schema':

        ff.copy_folder(schema_directory[0], schema_folder)
        second_schema_folder = os.path.join(initial_processing_output, 'second_schema')
        ff.create_directory(second_schema_folder)
        ff.copy_folder(schema_directory[1], second_schema_folder)

        
        pf.print_message(f'First schema: {schema_directory[0]}', 'info')
        pf.print_message(f'Second schema: {schema_directory[1]}', 'info')

        # Calculate the frequency of the cds in genomes for the first schema and second schema
        frequency_in_genomes, frequency_in_genomes_second_schema, _, _, _, dropped_alleles, _ = calculate_frequency(schema_folder, allelecall_directory[0], second_schema_folder, allelecall_directory[1], 
                                                None, 
                                                None, 
                                                None, 
                                                run_mode, None)
        
        ff.merge_folders(schema_directory[0], schema_directory[1], schema_folder, cpu, bsr, translation_table)

        pf.print_message('Prepare loci files for Blast and count frequencies.', 'info')
        (all_nucleotide_sequences,
        master_file_path,
        trans_paths,
        to_blast_paths,
        all_alleles,
        group_reps_ids,
        group_alleles_ids,
        to_run_against,
        new_max_hits,
        seqid_file_dict) = cof.prepare_loci(schema_folder,
                                            constants,
                                            processing_mode,
                                            initial_processing_output)


    # Main part that is the same for all run modes
    # =============================================

    # Create BLAST db for the schema DNA sequences.
    pf.print_message("Creating BLASTp database...", "info")
    # Get the path to the makeblastdb executable.
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    blast_db_prot: str = os.path.join(blast_db, 'Blast_db_proteins')
    bf.make_blast_db(makeblastdb_exec, master_file_path, blast_db_prot, 'prot')


    # Run the BLASTn and BLASTp
    representative_blast_results: tp.BlastDict
    representative_blastn_results: tp.BlastDict
    representative_blast_results, representative_blastn_results, loci_too_big = cof.run_blasts(blast_db_prot,
                                                                                                all_alleles,
                                                                                                trans_paths,
                                                                                                to_blast_paths,
                                                                                                to_run_against,
                                                                                                blast_output,
                                                                                                new_max_hits,
                                                                                                seqid_file_dict,
                                                                                                constants,
                                                                                                reps_kmers_sim if run_mode == 'unclassified_cds' else None,
                                                                                                frequency_in_genomes,
                                                                                                frequency_in_genomes_second_schema if run_mode == 'schema_vs_schema' else None,
                                                                                                cpu,
                                                                                                output_d)


    pf.print_message("Filtering BLAST results into classes...", "info")
    # Separate results into different classes.
    classes_outcome: Tuple[str] = cof.separate_blast_results_into_classes(representative_blast_results, representative_blastn_results,
                                                           constants, ct.CLASSES_OUTCOMES)
    
    
    # Sort each entry based on their assigned classes
    sorted_blast_dict: tp.BlastDict = cof.sort_blast_results_by_classes(representative_blast_results,
                                                          classes_outcome)
    # Process the results_outcome dict and write individual classes to TSV file.
    processed_results: tp.ProcessedResults
    count_results_by_class: tp.CountResultsByClass
    count_results_by_class_with_inverse: tp.CountResultsByClassWithInverse
    reps_and_alleles_ids: tp.RepsAndAllelesIds
    drop_mark: List[str]
    all_relationships: tp.AllRelationships
    # Process and extract relevant information from the blast results
    pf.print_message("Processing classes...", "info")
    (processed_results,
     count_results_by_class,
     count_results_by_class_with_inverse,
     reps_and_alleles_ids,
     drop_mark,
     all_relationships) = cof.process_classes(sorted_blast_dict, representative_blastn_results,
                                classes_outcome,
                                all_alleles)

    count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)
    # Extract which clusters are to maintain and to display to user.
    clusters_to_keep: tp.MergedAllClasses
    dropped_loci_ids: Set[str]
    clusters_to_keep_1a, clusters_to_keep, dropped_loci_ids = cof.extract_clusters_to_keep(classes_outcome, count_results_by_class, drop_mark)
 
    
    # Add the loci/new_loci IDs of the 1a joined clusters to the clusters_to_keep
    clusters_to_keep_1a_renamed: Dict[str, List[str]] = {values[0]: values for key, values in clusters_to_keep_1a.items()}

    # Merge classes
    merged_all_classes: tp.MergedAllClasses = {'1a': clusters_to_keep_1a_renamed.copy()}
    merged_all_classes.update(clusters_to_keep)
    if run_mode == 'unclassified_cds':
        updated_frequency_in_genomes: Dict[str, int] = ccf.update_frequencies_in_genomes(clusters_to_keep_1a_renamed,  frequency_in_genomes)
    
        # Open dict to store IDs of the reps and alleles
        group_reps_ids = {}
        group_alleles_ids = {}
        # Count the number of reps and alleles again because clusters were joined
        group_reps_ids, group_alleles_ids = cof.count_number_of_reps_and_alleles(merged_all_classes,
                                                                                processing_mode,
                                                                                all_alleles,
                                                                                dropped_loci_ids,
                                                                                group_reps_ids,
                                                                                group_alleles_ids)

        pf.print_message("Adding remaining clusters that didn't match by BLASTn...", "info")
        # Add cluster not matched by BLASTn
        all_matched_clusters: List[str] = itf.flatten_list([v for v in {key: value for key, value in clusters_to_keep.items() if key != '1a'}.values()]) + itf.flatten_list([values for values in clusters_to_keep_1a_renamed.values()])
        clusters_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in all_alleles.keys() if cluster not in all_matched_clusters])

    processed_drop: List[str] = []
    # Add Ids of the dropped cases due to frequency during classification
    cof.add_cds_to_dropped_cds(dropped_loci_ids,
                            dropped_alleles,
                            clusters_to_keep,
                            clusters_to_keep_1a_renamed,
                            all_alleles,
                            'Dropped_due_to_smaller_genome_presence_than_matched_cluster',
                            processed_drop)

    pf.print_message("Extracting results...", "info")
    related_clusters: tp.RelatedClusters
    recommendations: tp.Recomendations
    # Extract the results from the processed results
    related_clusters, recommendations = cof.extract_results(processed_results,
                                                            count_results_by_class,
                                                            frequency_in_genomes,
                                                            frequency_in_genomes_second_schema if run_mode == 'schema_vs_schema' else None, 
                                                            merged_all_classes,
                                                            dropped_loci_ids,
                                                            classes_outcome)


    pf.print_message("Writing count_results_by_cluster.tsv, related_matches.tsv files and recommendations.tsv...", "info")
    # Write the results to files and return the paths to the files.
    reverse_matches: bool = True
    (related_matches_path,
     count_results_by_cluster_path,
     recommendations_file_path,
     count_classes_final) = cof.write_recommendations_summary_results(to_blast_paths,
                                                                            related_clusters,
                                                                            count_results_by_class_with_inverse,
                                                                            group_reps_ids,
                                                                            group_alleles_ids,
                                                                            frequency_in_genomes,
                                                                            frequency_in_genomes_second_schema if run_mode == 'schema_vs_schema' else None,
                                                                            recommendations,
                                                                            reverse_matches,
                                                                            classes_outcome,
                                                                            output_d)

    # Get all of the CDS that matched with loci or loci matched with loci
    is_matched: Dict[str, Any]
    is_matched_alleles: Dict[str, Any]
    is_matched, is_matched_alleles = cof.get_matches(all_relationships,
                                                    merged_all_classes,
                                                    sorted_blast_dict)

    pf.print_message("Writing classes and cluster results to files...", "info")
    report_file_path: str = os.path.join(blast_results, 'blast_all_matches.tsv')
    # Write all of the alignments results to a file.
    cof.alignment_dict_to_file(representative_blast_results,
                               report_file_path,
                               'w')
    cof.alignment_dict_to_file(representative_blastn_results,
                               report_file_path,
                               'w')
    # Write the processed results to a file alignments by clusters and classes.
    cof.write_processed_results_to_file(merged_all_classes,
                                    representative_blast_results,
                                    classes_outcome,
                                    all_alleles,
                                    is_matched,
                                    is_matched_alleles,
                                    blast_results)
    merged_classes_6 = {k: merged_all_classes[k] for k in {'6'}}
    cof.write_processed_results_to_file(merged_classes_6,
                                    representative_blastn_results,
                                    ['6'],
                                    all_alleles,
                                    is_matched,
                                    is_matched_alleles,
                                    blast_results)

    if run_mode == 'unclassified_cds':
        pf.print_message("Updating IDs and saving changes in cds_id_changes.tsv...", "info")
        # Update the IDs and save the changes in a file.
        ccf.update_ids_and_save_changes(merged_all_classes,
                                    all_alleles,
                                    cds_original_ids,
                                    dropped_alleles,
                                    all_nucleotide_sequences,
                                    results_output)
        pf.print_message("Writing dropped CDSs to file...", "info")
        # Write the dropped CDS to a file.
        ccf.write_dropped_cds_to_file(dropped_alleles, results_output)

    pf.print_message("Writing dropped possible new loci to file...", "info")
    # Write the dropped possible new loci to a file.
    drop_possible_loci_output = cof.write_dropped_possible_new_loci_to_file(dropped_loci_ids,
                                                                        dropped_alleles,
                                                                        output_d)
    # Write alot of alleles file
    alot_of_alleles_file = os.path.join(output_d, 'alot_of_alleles.txt')
    with open(alot_of_alleles_file, 'w') as alot_file:
        for loci in loci_too_big:
            alot_file.write(f"{loci}\n")

    # Print the classification results
    cof.print_classifications_results(merged_all_classes,
                                        dropped_loci_ids,
                                        to_blast_paths,
                                        all_alleles,
                                        count_classes_final)
    # Graphs are only created for unclassified CDS (see if needed for schema)
    if run_mode == 'unclassified_cds':
        pf.print_message("Writing temporary fastas to file...", "info")
        temp_fastas_folder: str = ccf.write_fastas_to_files(all_alleles, all_nucleotide_sequences, output_d)

        pf.print_message("Creating graphs for the BLAST results...", "info")
        cds_size_dicts: Dict[str, Any] = {'IDs': cds_size.keys(),
                        'Size': cds_size.values()}
        cds_translation_size_dicts: Dict[str, Any] = {'IDs': cds_size.keys(),
                                    'Size': [int(cds/3) for cds in cds_size.values()]}
        cof.create_graphs(report_file_path,
                    results_output,
                    'All_of_CDS_graphs',
                    [[cds_size_dicts, 'histogram', "Nucleotide Size", 'Size', 'CDS'],
                    [cds_translation_size_dicts, 'histogram','Protein Size' , 'Size', 'CDS']])
        
        for file in ff.get_paths_in_directory(os.path.join(blast_results, 'blast_results_by_class'), 'files'):
            cof.create_graphs(file,
                        results_output,
                        f"graphs_class_{os.path.basename(file).split('_')[-1].replace('.tsv', '')}")


    # Annotate the recommendation outputs with the given annotation files using consolidate
    pf.print_message("")
    consolidated_annotations = os.path.join(output_d, "recommendations_annotations.tsv")
    if annotation_paths:
        files: List[str]
        files = [recommendations_file_path] + annotation_paths
        pf.print_message("Consolidating annoations...", "info")
        consolidated_annotations: str = cs.consolidate_annotations(files,
                                    False,
                                    consolidated_annotations)
        pf.print_message('Annotation consolidation successfully completed.', 'info')
        pf.print_message('')

    # Clean up temporary files
    if not no_cleanup:
        pf.print_message("Cleaning up temporary files...", "info")
        # Remove temporary files
        ff.cleanup(output_d, [related_matches_path,
                                alot_of_alleles_file,
                                count_results_by_cluster_path,
                                recommendations_file_path,
                                consolidated_annotations if annotation_paths else None,
                                drop_possible_loci_output,
                                temp_fastas_folder if run_mode == 'unclassified_cds' else None,
                                logf.get_log_file_path(gb.LOGGER)])


def main(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
        annotation_paths: List[str],
        alignment_ratio_threshold: float, 
        pident_threshold: float, clustering_sim_threshold: float, clustering_cov_threshold:float,
        genome_presence: int, absolute_size: int, translation_table: int,
        bsr: float, size_ratio: float, run_mode: str, processing_mode: str, cpu: int,
        no_cleanup: bool) -> None:
    """
    Main function to identify spurious genes in a schema.

    Parameters
    ----------
    schema_directory : List[str]
        Path to the schema directory.
    output_directory : str
        Path to the output directory.
    allelecall_directory : List[str]
        Path to the allele call directory.
    annotation_paths : List[str]
        Paths for the files with the annotations.
    alignment_ratio_threshold : float
        Threshold for alignment ratio.
    pident_threshold : float
        Threshold for percentage identity.
    clustering_sim_threshold : float
        Similarity threshold for clustering.
    clustering_cov_threshold : float
        Coverage threshold for clustering.
    genome_presence : int
        Minimum genome presence required.
    absolute_size : int
        Absolute size threshold.
    translation_table : int
        Genetic code used for translation.
    bsr : float
        BLAST Score Ratio value.
    size_ratio : float
        Size ratio threshold.
    run_mode : str
        Mode of running the process.
    processing_mode : str
        Mode of processing.
    cpu : int
        Number of CPU cores to use.
    no_cleanup : bool
        Flag to indicate whether to clean up temporary files.

    Returns
    -------
    None
    """
    if run_mode == 'unclassified_cds':
        temp_paths: List[str] = [os.path.join(allelecall_directory[0], "temp"), 
                                os.path.join(allelecall_directory[0], "unclassified_sequences.fasta")]
        if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
            sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                        "was run using --no-cleanup and --output-unclassified flag.")
    
    # Put all constants in one dict in order to decrease number of variables
    # used around.
    constants: List[Any] = [alignment_ratio_threshold, 
                pident_threshold,
                genome_presence,
                clustering_sim_threshold,
                clustering_cov_threshold,
                absolute_size,
                translation_table,
                bsr,
                size_ratio]


    identify_spurious_genes(schema_directory,
                output_directory,
                allelecall_directory,
                annotation_paths,
                constants,
                temp_paths if run_mode == 'unclassified_cds' else [],
                run_mode,
                processing_mode,
                cpu,
                bsr,
                translation_table,
                no_cleanup)
