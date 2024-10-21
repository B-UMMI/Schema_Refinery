#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from typing import Any, List, Dict, Optional

try:
    from utils import (core_functions as cof,
                       file_functions as ff,
                       sequence_functions as sf,
                       iterable_functions as itf,
                       blast_functions as bf,
                       linux_functions as lf,
                       classify_cds_functions as ccf,
                       schema_classification_functions as scf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (core_functions as cof,
                                        file_functions as ff,
                                        sequence_functions as sf,
                                        iterable_functions as itf,
                                        blast_functions as bf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf,
                                        schema_classification_functions as scf)

def create_directories(output_directory: str, run_mode: str) -> List[Optional[str]]:
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
    
    # Create initial processing output directory based on run mode
    if run_mode == 'unclassified_cds':
        initial_processing_output = os.path.join(output_directory, '1_CDS_processing')
        ff.create_directory(initial_processing_output)
        schema_folder = None
    else:
        initial_processing_output = os.path.join(output_directory, '1_schema_processing')
        schema_folder = os.path.join(initial_processing_output, 'new_schema')
        ff.create_directory(schema_folder)

    # Create BLAST processing directories
    blast_output = os.path.join(output_directory, '2_BLAST_processing')
    blastn_output = os.path.join(blast_output, '1_BLASTn_processing')
    blast_db = os.path.join(blastn_output, 'blast_db_nucl')
    ff.create_directory(blast_db)

    # Create representatives BLASTn folder if run mode is 'unclassified_cds'
    if run_mode == 'unclassified_cds':
        representatives_blastn_folder = os.path.join(blastn_output, 'cluster_representatives_fastas')
        ff.create_directory(representatives_blastn_folder)
    else:
        representatives_blastn_folder = None
        
    # Create results output directories
    results_output = os.path.join(output_directory, '3_processing_results')
    blast_results = os.path.join(results_output, 'blast_results')
    ff.create_directory(blast_results)
    
    return [initial_processing_output, schema_folder, blast_output, blastn_output, blast_db, representatives_blastn_folder, results_output, blast_results]


def identify_spurious_genes(schema_directory: str, output_directory: str, allelecall_directory: str, possible_new_loci: str, 
                            constants: List[Any], temp_paths: List[str], run_mode: str, processing_mode: str, cpu: int) -> None:
    """
    Identify spurious genes in the given schema.

    Parameters
    ----------
    schema_directory : str
        Path to the schema directory.
    output_directory : str
        Path to the output directory.
    allelecall_directory : str
        Path to the allele call directory.
    possible_new_loci : str
        Path to possible new loci.
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

    Returns
    -------
    None
    """
    
    # Create directories structure.
    (initial_processing_output,
     schema_folder,
     blast_output,
     blastn_output,
     blast_db,
     representatives_blastn_folder,
     results_output,
     blast_results) = create_directories(output_directory, run_mode)
    
    # Process unclassified CDS, retrieving and clustering
    if run_mode == 'unclassified_cds':
        temp_folder: str = temp_paths[0]
        file_path_cds: str = temp_paths[1]

        # Verify if the dataset is small, if it is, keep minimum genomes in which
        # specific CDS cluster is present to 5 if not to 1% of the dataset size.
        if not constants[2]:
            ccf.set_minimum_genomes_threshold(temp_folder, constants)

        print("Identifying CDS not present in the schema...")
        # Get dict with CDS ids as key and sequence as values.
        all_nucleotide_sequences: Dict[str, str] = sf.fetch_fasta_dict(file_path_cds, True)
        
        # Make IDS universally usable
        for key, value in list(all_nucleotide_sequences.items()):
            all_nucleotide_sequences[itf.replace_by_regex(key, '_', '-')] = all_nucleotide_sequences.pop(key)

        print("\nFiltering missing CDS in the schema...")
        cds_size: Dict[str, int]
        dropped_alleles: Dict[str, str]
        cds_size, all_nucleotide_sequences, dropped_alleles = ccf.filter_cds_by_size(all_nucleotide_sequences, constants[5])

        # Count the number of CDS not present in the schema and write CDS sequence
        # into a FASTA file.
        frequency_cds: Dict[str, int] = {}
        cds_presence_in_genomes: Dict[str, int] = {}

        print("Identifying CDS present in the schema and counting frequency of missing CDSs in the genomes...")
        cds_present: Dict[str, str]
        cds_present, frequency_cds, cds_presence_in_genomes = ccf.process_cds_not_present(initial_processing_output,
                                                                                          temp_folder,
                                                                                          all_nucleotide_sequences)

        print("\nTranslate and deduplicate CDS...")
        all_translation_dict: Dict[str, str]
        protein_hashes: Dict[str, str]
        cds_translation_size: Dict[str, int]
        all_translation_dict, protein_hashes, cds_translation_size = ccf.translate_and_deduplicate_cds(
                                                                                                    all_nucleotide_sequences,
                                                                                                    initial_processing_output,
                                                                                                    constants
                                                                                                )

        print("\nExtracting minimizers for the translated sequences and clustering...")
        all_translation_dict = ccf.remove_dropped_cds(all_translation_dict, dropped_alleles, protein_hashes)
        all_translation_dict = ccf.sort_by_protein_size(all_translation_dict)
        all_alleles: Dict[str, List[str]]
        reps_sequences: Dict[str, str]
        reps_groups: Dict[str, List[str]]
        prot_len_dict: Dict[str, int]
        all_alleles, reps_sequences, reps_groups, prot_len_dict = ccf.cluster_by_minimizers(all_translation_dict, constants)
        all_alleles = ccf.reformat_clusters(all_alleles, protein_hashes)

        # Calculate the total number of clusters
        total_number_clusters: int = len(all_alleles)
        print(f"{len(all_translation_dict)} unique proteins have been clustered into {total_number_clusters} clusters.")

        # Calculate the number of singleton clusters
        singleton_clusters: int = len([cluster for cluster in all_alleles.values() if len(cluster) == 1])
        print(f"\tOut of those clusters, {singleton_clusters} are singletons")

        # Calculate the number of clusters with more than one CDS
        multi_cds_clusters: int = total_number_clusters - singleton_clusters
        print(f"\tOut of those clusters, {multi_cds_clusters} have more than one CDS.")
        
        print("\nFiltering clusters...")
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

        print(f"After filtering by CDS frequency in the genomes (>= {constants[2]}),"
            f" out of {total_number_clusters} clusters, {len(filtered_alleles)} remained.")

        # Update all_alleles with the filtered results
        all_alleles = filtered_alleles

        print("\nRetrieving kmers similarity and coverage between representatives...")
        reps_translation_dict: Dict[str, str] = ccf.get_representative_translation_dict(all_translation_dict, all_alleles)
        reps_kmers_sim: Dict[str, float] = ccf.calculate_kmers_similarity(reps_translation_dict, reps_groups, prot_len_dict)

        print("\nReplacing CDSs IDs with the cluster representative ID...")
        cds_original_ids: Dict[str, str] = ccf.replace_ids_in_clusters(all_alleles,
                                                    frequency_cds,
                                                    dropped_alleles,
                                                    all_nucleotide_sequences,
                                                    prot_len_dict,
                                                    all_translation_dict,
                                                    protein_hashes,
                                                    cds_presence_in_genomes,
                                                    reps_kmers_sim)

        to_blast_paths: List[str]
        master_file_path: str
        to_blast_paths, master_file_path = ccf.create_blast_files(representatives_blastn_folder, all_alleles, all_nucleotide_sequences, processing_mode)

    # Process loci
    else:
        if possible_new_loci:
            ff.merge_folders(schema_directory, possible_new_loci, schema_folder)
        else:
            ff.copy_folder(schema_directory, schema_folder)

        dropped_alleles = {}

        (all_nucleotide_sequences,
        master_file_path,
        all_translation_dict,
        frequency_in_genomes,
        to_blast_paths,
        all_alleles,
        cds_present,
        group_reps_ids,
        group_alleles_ids) = scf.process_new_loci(schema_folder, allelecall_directory, constants, processing_mode, initial_processing_output)

    # Create BLAST db for the schema DNA sequences.
    print("\nCreating BLASTn database...")
    # Get the path to the makeblastdb executable.
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    blast_db_nuc: str = os.path.join(blast_db, 'Blast_db_nucleotide')
    bf.make_blast_db(makeblastdb_exec, master_file_path, blast_db_nuc, 'nucl')
    
    print("\nRunning BLASTn...")
    # Run the BLASTn and BLASTp
    (representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _) = cof.run_blasts(blast_db_nuc,
                        all_alleles,
                        all_translation_dict,
                        to_blast_paths,
                        blast_output,
                        constants,
                        cpu,
                        all_alleles)
    
    # Add various results to the dict
    cof.add_items_to_results(representative_blast_results,
                         reps_kmers_sim if run_mode == 'unclassified_cds' else None,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_in_genomes,
                         [True, True])

    print("\nFiltering BLAST results into classes...")
    # Separate results into different classes.
    classes_outcome = cof.separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    
    print("\nProcessing classes...")
    sorted_blast_dict = cof.sort_blast_results_by_classes(representative_blast_results,
                                                          classes_outcome)
    # Process the results_outcome dict and write individual classes to TSV file.
    (processed_results,
     count_results_by_class,
     count_results_by_class_with_inverse,
     reps_and_alleles_ids,
     drop_mark,
     all_relationships) = cof.process_classes(sorted_blast_dict,
                                classes_outcome,
                                all_alleles)

    count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)
    # Extract which clusters are to maintain and to display to user.
    clusters_to_keep, dropped_loci_ids = cof.extract_clusters_to_keep(classes_outcome, count_results_by_class, drop_mark)
    
    # Add the loci/new_loci IDs of the 1a joined clusters to the clusters_to_keep
    clusters_to_keep['1a'] = {values[0]: values for key, values in clusters_to_keep['1a'].items()}
    if run_mode == 'unclassified_cds':
        updated_frequency_in_genomes = ccf.update_frequencies_in_genomes(clusters_to_keep, frequency_in_genomes)
    
        # Open dict to store IDs of the reps and alleles
        group_reps_ids = {}
        group_alleles_ids = {}

        cof.count_number_of_reps_and_alleles(clusters_to_keep,
                                            all_alleles,
                                            dropped_loci_ids,
                                            group_reps_ids,
                                            group_alleles_ids)

        print("\nAdd remaining cluster that didn't match by BLASTn...")
        # Add cluster not matched by BLASTn
        all_matched_clusters = itf.flatten_list([v for v in {key: value for key, value in clusters_to_keep.items() if key != '1a'}.values()]) + itf.flatten_list([values for values in clusters_to_keep['1a'].values()])
        clusters_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in all_alleles.keys() if cluster not in all_matched_clusters])

    processed_drop = []
    # Add Ids of the dropped cases due to frequency during classification
    cof.add_cds_to_dropped_cds(dropped_loci_ids,
                            dropped_alleles,
                            clusters_to_keep,
                            all_alleles,
                            'Dropped_due_to_smaller_genome_presence_than_matched_cluster',
                            processed_drop)

    print("\nFiltering problematic possible new loci based on NIPHS and NIPHEMS presence and"
        " writing results to niphems_and_niphs_groups.tsv...")
    # Identify problematic possible new loci based on NIPHs and NIPHEMs presence
    dropped_problematic = cof.identify_problematic_new_loci(clusters_to_keep,
                                                         all_alleles,
                                                         cds_present,
                                                         all_nucleotide_sequences,
                                                         constants,
                                                         results_output)
    # For unclassified cds remove NIPHS and NIPHEMS possible new loci
    # For Schema keep them because the user should make the choice to remove or keep
    # them, especially if there are joined clusters.
    if run_mode == 'unclassified_cds':
        # Add Ids of the dropped cases due to during NIPH and NIPHEMs
        dropped_loci_ids.update(dropped_problematic)
        # Add Ids of the dropped cases due to frequency during NIPH and NIPHEMs
        # classification
        cof.add_cds_to_dropped_cds(dropped_loci_ids,
                                dropped_alleles,
                                clusters_to_keep,
                                all_alleles,
                                'Dropped_due_high_presence_of_NIPHs_and_NIPHEMs_in_genomes',
                                processed_drop)

        # Remove from all relevant dicts
        ccf.remove_dropped_cds_from_analysis(dropped_alleles,
                                            all_nucleotide_sequences,
                                            all_translation_dict,
                                            protein_hashes)

    print("\nExtracting results...")
    related_clusters, recommendations = cof.extract_results(processed_results,
                                                            count_results_by_class,
                                                            frequency_in_genomes,
                                                            clusters_to_keep,
                                                            dropped_loci_ids,
                                                            classes_outcome)

    print("\nWriting count_results_by_cluster.tsv, related_matches.tsv files"
          " and recommendations.tsv...")
    reverse_matches = True
    cof.write_blast_summary_results(related_clusters,
                                count_results_by_class_with_inverse,
                                group_reps_ids,
                                group_alleles_ids,
                                frequency_in_genomes,
                                recommendations,
                                reverse_matches,
                                classes_outcome,
                                results_output)
    
    if run_mode != 'unclassified_cds':
        add_group_column = True
    else:
        add_group_column = False

    # Get all of the CDS that matched with loci
    is_matched, is_matched_alleles = cof.get_matches(all_relationships,
                                                    clusters_to_keep,
                                                    sorted_blast_dict)
    print("\nWriting classes and cluster results to files...")
    report_file_path = os.path.join(blast_results, 'blast_all_matches.tsv')
    # Write all of the BLASTn results to a file.
    cof.alignment_dict_to_file(representative_blast_results,
                               report_file_path,
                               'w',
                               add_group_column)

    cof.write_processed_results_to_file(clusters_to_keep,
                                    representative_blast_results,
                                    classes_outcome,
                                    all_alleles,
                                    is_matched,
                                    is_matched_alleles,
                                    blast_results)

    if run_mode == 'unclassified_cds':
        print("\nUpdating IDs and saving changes in cds_id_changes.tsv...")
        ccf.update_ids_and_save_changes(clusters_to_keep,
                                    all_alleles,
                                    cds_original_ids,
                                    dropped_alleles,
                                    all_nucleotide_sequences,
                                    results_output)
        print("\nWriting dropped CDSs to file...")
        ccf.write_dropped_cds_to_file(dropped_alleles, results_output)

    print("\nWriting dropped possible new loci to file...")
    cof.write_dropped_possible_new_loci_to_file(dropped_loci_ids,
                                                dropped_alleles,
                                                results_output)

    cof.print_classifications_results(clusters_to_keep,
                                        dropped_loci_ids,
                                        to_blast_paths,
                                        all_alleles)

    if run_mode == 'unclassified_cds':
        print("\nWriting temp loci file...")
        ccf.write_temp_loci(clusters_to_keep, all_nucleotide_sequences, all_alleles, results_output)

        print("\nCreate graphs for the BLAST results...")
        cds_size_dicts = {'IDs': cds_size.keys(),
                        'Size': cds_size.values()}
        cds_translation_size_dicts = {'IDs': cds_size.keys(),
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

def main(schema_directory, output_directory, allelecall_directory, possible_new_loci, alignment_ratio_threshold, 
        pident_threshold, clustering_sim, clustering_cov, genome_presence,
        absolute_size, translation_table, bsr, problematic_proportion, size_ratio, run_mode, processing_mode, cpu):
    
    temp_paths: List[str] = [os.path.join(allelecall_directory, "temp"), 
                            os.path.join(allelecall_directory, "unclassified_sequences.fasta"),
                            os.path.join(allelecall_directory, "missing_classes.fasta")]
    # Put all constants in one dict in order to decrease number of variables
    # used around.
    constants: List[Any] = [alignment_ratio_threshold, 
                pident_threshold,
                genome_presence,
                clustering_sim,
                clustering_cov,
                absolute_size,
                translation_table,
                bsr,
                problematic_proportion,
                size_ratio]
    
    if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
        sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                    "was run using --no-cleanup and --output-unclassified flag.")

    unclassified_cds_output = os.path.join(output_directory, "unclassified_cds_processing" if run_mode == 'unclassified_cds' else 'schema_processing') 
    identify_spurious_genes(schema_directory,
                unclassified_cds_output,
                allelecall_directory,
                possible_new_loci,
                constants,
                temp_paths,
                run_mode,
                processing_mode,
                cpu)