#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

try:
    from utils import (core_functions as cof,
                       file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       iterable_functions as itf,
                       kmers_functions as kf,
                       blast_functions as bf,
                       linux_functions as lf,
                       classify_cds_functions as ccf)

    from RefineSchema import (UnclassifiedCDS as uc)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (core_functions as cof,
                                        file_functions as ff,
                                        sequence_functions as sf,
                                        clustering_functions as cf,
                                        iterable_functions as itf,
                                        kmers_functions as kf,
                                        blast_functions as bf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf)
    
    from SchemaRefinery.RefineSchema import (SpuriousLoci as sl)

def classify_cds(schema, output_directory, allelecall_directory, constants, temp_paths, cpu):

    temp_folder = temp_paths[0]
    file_path_cds = temp_paths[1]
    #missing_classes_fastas = temp_paths[2]

    # Verify if the dataset is small, if it is, keep minimum genomes in which
    # specific CDS cluster is present to 5 if not to 1% of the dataset size.
    count_genomes_path = os.path.join(temp_folder, '1_cds_prediction')
    if not constants[2]:
        number_of_genomes = len(ff.get_paths_in_directory_with_suffix(count_genomes_path, '.fasta'))
        if number_of_genomes <= 20:
            constants[2] = 5
        else:
            constants[2] = round(number_of_genomes * 0.01)
    # Get all of the genomes IDs.
    genomes_ids = ff.get_paths_in_directory_with_suffix(count_genomes_path, '.fasta')

    print("Identifying CDS not present in the schema...")
    # Get dict with CDS ids as key and sequence as values.
    not_included_cds = sf.fetch_fasta_dict(file_path_cds, True)
    #Make IDS universally usable
    for key, value in list(not_included_cds.items()):
        not_included_cds[itf.replace_by_regex(key, '_', '-')] = not_included_cds.pop(key)

    print("\nFiltering missing CDS in the schema...")
    # Count CDS size
    cds_size = {}
    for key, sequence in not_included_cds.items():
        cds_size.setdefault(key, len(str(sequence)))

    dropped_cds = {}
    total_cds = len(not_included_cds)
    print(f"\nIdentified {total_cds} valid CDS not present in the schema.")
    # Filter by size.
    if constants[5]:
        for key, values in list(not_included_cds.items()):
            if len(values) < constants[5]:
                dropped_cds.setdefault(key, 'Dropped_due_to_cds_size')
                del not_included_cds[key]
        print(f"{len(not_included_cds)}/{total_cds} have size greater or equal to {constants[5]} bp.")
    else:
        constants[5] = 0
        print("No size threshold was applied to the CDS filtering.")

    # Create directories.
    ff.create_directory(output_directory)

    cds_output = os.path.join(output_directory, '1_CDS_processing')
    ff.create_directory(cds_output)
    # This file contains unique CDS.
    cds_not_present_file_path = os.path.join(cds_output, 'CDS_not_found.fasta')

    # Count the number of CDS not present in the schema and write CDS sequence
    # into a FASTA file.
    frequency_cds = {}
    cds_presence_in_genomes = {}

    print("Identifying CDS present in the schema and counting frequency of missing CDSs in the genomes...")
    cds_present = os.path.join(temp_folder,"2_cds_preprocess/cds_deduplication/distinct.hashtable")
    # Get dict of CDS and their sequence hashes.
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)

    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in list(not_included_cds.items()):
            cds_not_found.write(f">{id_}\n{str(sequence)}\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus the first (subtract the first CDS genome).
            if hashed_seq in decoded_sequences_ids:
                #Count frequency.
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq][1:])
                cds_presence_in_genomes.setdefault(id_, decoded_sequences_ids[hashed_seq][1:])
            else:
                frequency_cds[id_] = 0

    print("\nTranslate and deduplicate CDS...")
    # Translate the CDS and find unique proteins using hashes, the CDS with
    # the same hash will be added under that hash in protein_hashes.
    cds_not_present_trans_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    cds_not_present_untrans_file_path = os.path.join(cds_output, "CDS_not_found_untranslated.fasta")
    # Translate and deduplicate protein sequences.
    cds_translation_dict, protein_hashes, _ = sf.translate_seq_deduplicate(not_included_cds,
                                                                           cds_not_present_trans_file_path,
                                                                           cds_not_present_untrans_file_path,
                                                                           constants[5],
                                                                           True,
                                                                           constants[6],
                                                                           True)
    # Count translation sizes.
    cds_translation_size = {}
    for key, sequence in cds_translation_dict.items():
        cds_translation_size.setdefault(key, len(sequence))

    # Print additional information about translations and deduplications.
    print(f"\n{len(cds_translation_dict)}/{len(not_included_cds)} unique protein translations.")
    
    print("\nIdentify problematics CDSs...")
    [niphems_presence_in_genome,
     niphs_in_genomes,
     niphs_presence_in_genomes] = ccf.identify_problematic_cds(cds_presence_in_genomes,
                                                           cds_translation_dict,
                                                           protein_hashes,
                                                           not_included_cds,
                                                           cds_output,
                                                           constants[7],
                                                           dropped_cds,
                                                           constants[9],
                                                           cpu)

    print("\nExtracting minimizers for the translated sequences and clustering...")
    # Create variables to store clustering info.
    reps_groups = {}
    clusters = {}
    reps_sequences = {}

    # Remove dropped_cds
    for key, protein in list(cds_translation_dict.items()):
        if key in dropped_cds:
            protein_hash = itf.identify_string_in_dict_get_key(key, protein_hashes)
            same_protein_id = protein_hashes[protein_hash]
            if key == same_protein_id[0]:
                protein_hashes[protein_hash].remove(key)
                if not protein_hashes[protein_hash]:
                    del protein_hashes[protein_hash]
                    continue
                new_id = protein_hashes[protein_hash][0]
                cds_translation_dict[new_id] = cds_translation_dict.pop(key)
            else:
                protein_hashes[protein_hash].remove(key)

    cds_translation_dict = {k: v for k, v in cds_translation_dict.items() if k not in dropped_cds}

    # Sort by size of proteins.
    cds_translation_dict = {k: v for k, v in sorted(cds_translation_dict.items(),
                                                    key=lambda x: len(x[1]),
                                                    reverse=True)}

    # Cluster by minimizers.
    [clusters, reps_sequences, 
     reps_groups, prot_len_dict] = cf.minimizer_clustering(cds_translation_dict,
                                                           5,
                                                           5,
                                                           True,
                                                           1, 
                                                           clusters,
                                                           reps_sequences, 
                                                           reps_groups,
                                                           1,
                                                           constants[3], 
                                                           constants[4],
                                                           True,
                                                           constants[9])

    # Reformat the clusters output, we are interested only in  the ID of cluster members.
    clusters = {cluster_rep: [value[0] for value in values]
                for cluster_rep, values in clusters.items()}
    # For protein hashes get only those that have more than one CDS.
    filtered_protein_hashes = {hash_prot: cds_ids for hash_prot, cds_ids in protein_hashes.items()
                      if len(cds_ids) > 1}
    # Add also the unique CDS ID to clusters that have the same protein as representative.
    for cluster_rep, values in list(clusters.items()):
        for cds_id in list(values):
            protein_hash = itf.identify_string_in_dict_get_key(cds_id, filtered_protein_hashes)
            if protein_hash:
                clusters[cluster_rep] += filtered_protein_hashes[protein_hash][1:]

    total_number_clusters = len(clusters)
    print(f"{len(cds_translation_dict)} unique proteins have been clustered into {total_number_clusters} clusters.")
    singleton_cluster = len([cluster for cluster in clusters if len(cluster) == 1])
    print(f"\tOut of those clusters, {singleton_cluster} are singletons")
    print(f"\tOut of those clusters, {total_number_clusters - singleton_cluster} have more than one CDS.")
    
    print("\nFiltering clusters...")
    # Get frequency of cluster.
    frequency_in_genomes = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Add reason for filtering out CDS.
    dropped_cds.update({cds_id: 'Dropped_due_to_cluster_frequency_filtering' for cds_id in itf.flatten_list([clusters[rep] for rep in clusters if frequency_in_genomes[rep] < constants[2]])})

    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_in_genomes[rep] >= constants[2]}

    print(f"After filtering by CDS frequency in the genomes (>= {constants[2]}),"
           f" out of {total_number_clusters} clusters, {len(clusters)} remained.")
     
    intial_length = len(clusters)
    
    if intial_length != len(clusters):
         print(f"After filtering by CDS frequency in the genomes (>= {len(genomes_ids)}),"
               f" out of {intial_length} clusters, {len(clusters)} remained.")

    print("\nRetrieving kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    # Get the representatives protein sequence.
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items()
                             if rep_id.split('_')[0] in clusters}
    # Sort the representative translation dict from largest to smallest.
    reps_translation_dict = {k: v for k, v in sorted(reps_translation_dict.items(),
                                                     key=lambda x: len(x[1]),
                                                     reverse=True)}
    # recalculate the sim and cov between reps, get all of the values, so threshold
    # is set to 0.
    for cluster_id in reps_translation_dict:
        kmers_rep = set(kf.determine_minimizers(reps_translation_dict[cluster_id],
                                                5,
                                                5,
                                                1,
                                                True,
                                                True))
        
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep,
                                                               reps_groups,
                                                               0,
                                                               0,
                                                               prot_len_dict,
                                                               cluster_id,
                                                               5,
                                                               False)
    
        reps_kmers_sim[cluster_id] = {match_values[0]: match_values[1:]
                                      for match_values in reps_kmers_sim[cluster_id]}

    print("\nReplacing CDSs IDs with the cluster representative ID...")
    cds_original_ids = ccf.replace_ids_in_clusters(clusters,
                                                   frequency_cds,
                                                   dropped_cds,
                                                   not_included_cds,
                                                    prot_len_dict,
                                                    cds_translation_dict,
                                                    protein_hashes,
                                                    cds_presence_in_genomes,
                                                    niphems_presence_in_genome,
                                                    niphs_presence_in_genomes,
                                                    niphs_in_genomes,
                                                    reps_kmers_sim)

    # Create directories.
    blast_output = os.path.join(output_directory, '2_BLAST_processing')
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, '1_BLASTn_processing')
    ff.create_directory(blastn_output)
    # Create directory and files path where to write FASTAs.
    representatives_blastn_folder = os.path.join(blastn_output,
                                                'cluster_representatives_fastas')
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  'all_clusters.fasta')
    # Write files for BLASTn.
    rep_paths_nuc = {}
    # Write master file for the representatives.
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for members in clusters.values():
            for member in members:
                all_fasta.write(f">{member}\n{str(not_included_cds[member])}\n")

            cluster_rep_id = members[0]
            rep_fasta_file = os.path.join(representatives_blastn_folder,
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths_nuc[cluster_rep_id] = rep_fasta_file
            # Write the representative FASTA file.
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.write(f">{cluster_rep_id}\n{str(not_included_cds[cluster_rep_id])}\n")
    
    # Create BLAST db for the schema DNA sequences.
    print("\nCreating BLASTn database for the unclassified and missed CDSs...")
    # Get the path to the makeblastdb executable.
    makeblastdb_exec = lf.get_tool_path('makeblastdb')
    blast_db = os.path.join(blastn_output, 'blast_db_nucl', 'blast_nucleotide_db')
    bf.make_blast_db(makeblastdb_exec, representatives_all_fasta_file, blast_db, 'nucl')

    # Run the BLASTn and BLASTp
    run_type = 'cds_vs_cds' # Set run type as cds_vs_cds
    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = cof.run_blasts(blast_db,
                        clusters,
                        cds_translation_dict,
                        rep_paths_nuc,
                        blast_output,
                        constants,
                        cpu,
                        clusters,
                        run_type)
    
    # Add various results to the dict
    cof.add_items_to_results(representative_blast_results,
                         reps_kmers_sim,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_in_genomes,
                         [True, True],
                         clusters)

    print("\nFiltering BLAST results into classes...")
    results_output = os.path.join(output_directory, '3_CDS_processing_results')
    ff.create_directory(results_output)
    blast_results = os.path.join(results_output, 'blast_results')
    ff.create_directory(blast_results)
    report_file_path = os.path.join(blast_results, 'blast_all_matches.tsv')
    
    # Separate results into different classes.
    classes_outcome = cof.separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    # Write all of the BLASTn results to a file.
    cof.alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("\nProcessing classes...")
    sorted_blast_dict = cof.sort_blast_results_by_classes(representative_blast_results,
                                                          classes_outcome)
    # Process the results_outcome dict and write individual classes to TSV file.
    [processed_results,
     count_results_by_class,
     count_results_by_class_with_inverse,
     reps_and_alleles_ids,
     drop_mark] = cof.process_classes(sorted_blast_dict,
                                classes_outcome,
                                clusters)

    count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)

    clusters_to_keep, drop_possible_loci = cof.extract_clusters_to_keep(classes_outcome, count_results_by_class, drop_mark)

    clusters_to_keep['1a'] = {values[0]: values for key, values in clusters_to_keep['1a'].items()}
    
    # Add new frequencies in genomes for joined groups
    # Update the changed clusters frequency from Joined CDSs
    updated_frequency_in_genomes = {}
    new_cluster_freq = {}
    for cluster_id, cluster_members in clusters_to_keep['1a'].items():
        new_cluster_freq[cluster_id] = 0
        for member in cluster_members:
            new_cluster_freq[(cluster_id)] += frequency_in_genomes[member]
        for member in cluster_members:
            updated_frequency_in_genomes[member] = new_cluster_freq[cluster_id]
    #Add all the others frequencies.
    updated_frequency_in_genomes.update(frequency_in_genomes)
    updated_frequency_in_genomes.update(new_cluster_freq)

    group_reps_ids = {}
    group_alleles_ids = {}
    cof.count_number_of_reps_and_alleles(clusters_to_keep,
                                         clusters,
                                         drop_possible_loci,
                                         group_reps_ids,
                                         group_alleles_ids)
    
    print("\nAdd remaining cluster that didn't match by BLASTn...")
    # Add cluster not matched by BLASTn
    all_matched_clusters = itf.flatten_list([v for v in {key: value for key, value in clusters_to_keep.items() if key != '1a'}.values()]) + itf.flatten_list([values for values in clusters_to_keep['1a'].values()])
    clusters_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in clusters.keys() if cluster not in all_matched_clusters])

    processed_drop = []
    # Add Ids of the dropped cases due to frequency during classification
    ccf.add_cds_to_dropped_cds(drop_possible_loci,
                               dropped_cds,
                               clusters_to_keep,
                               clusters,
                               'Dropped_due_to_smaller_genome_presence_than_matched_cluster', processed_drop)

    print("\nFiltering problematic probable new loci...")

    [proportion_of_niph_genomes,
     dropped_due_to_niphs_or_niphems,
     clusters_to_keep_all_members,
     clusters_to_keep_all_genomes] = ccf.identify_problematic_loci(niphems_presence_in_genome,
                                                                   niphs_in_genomes,
                                                                   niphs_presence_in_genomes,
                                                                    cds_presence_in_genomes,
                                                                    clusters_to_keep,
                                                                    clusters,
                                                                    constants[8], 
                                                                    dropped_cds,
                                                                    drop_possible_loci,
                                                                    results_output)
    
    # Add Ids of the dropped cases due to frequency during NIPH and NIPHEMs
    # classification
    ccf.add_cds_to_dropped_cds(drop_possible_loci,
                               dropped_cds,
                               clusters_to_keep,
                               clusters,
                               'Dropped_due_high_presence_of_NIPHs_and_NIPHEMs_in_genomes', processed_drop)

    # Remove from all releveant dicts
    ccf.remove_dropped_cds_from_analysis(dropped_cds,
                                         not_included_cds,
                                         niphems_presence_in_genome,
                                         cds_translation_dict,
                                         niphs_presence_in_genomes,
                                         protein_hashes,
                                         niphs_in_genomes)

    print("\nExtracting results...")
    all_relationships, related_clusters, recommendations = cof.extract_results(processed_results,
                                                                          count_results_by_class,
                                                                          frequency_in_genomes,
                                                                          clusters_to_keep,
                                                                          drop_possible_loci,
                                                                          classes_outcome)
    print("\nWritting count_results_by_cluster.tsv and related_matches.tsv files...")
    cof.write_blast_summary_results(related_clusters,
                                count_results_by_class_with_inverse,
                                group_reps_ids,
                                group_alleles_ids,
                                frequency_in_genomes,
                                recommendations,
                                True,
                                results_output)

    print("\nWritting classes and cluster results to files...")
    cof.write_processed_results_to_file(clusters_to_keep,
                                    representative_blast_results,
                                    classes_outcome,
                                    None,
                                    clusters,
                                    None,
                                    None,
                                    [False, False],
                                    blast_results)
    
    print("\nUpdating IDs and saving changes...")
    
    ccf.update_ids_and_save_changes(clusters_to_keep,
                                    clusters,
                                    cds_original_ids,
                                    dropped_cds,
                                    not_included_cds,
                                    results_output)

    cds_cases, loci_cases = cof.print_classifications_results(clusters_to_keep,
                                                              drop_possible_loci,
                                                              False,
                                                              clusters,
                                                              False,
                                                              run_type)

    print("\nWritting possible new loci Fastas...")
    [groups_paths_reps,
     groups_paths,
     groups_trans_reps_paths,
     groups_trans,
     reps_trans_dict_cds,
     trans_dict_cds,
     master_file,
     alleles] = cof.wrap_up_blast_results(clusters_to_keep,
                                              not_included_cds,
                                              clusters,
                                              results_output,
                                              constants,
                                              None,
                                              None,
                                              cds_cases,
                                              loci_cases,
                                              run_type)
    
    print("Identifying new representatives for possible new loci...")
    ccf.find_new_representatives(groups_trans_reps_paths,
                                 groups_trans,
                                 groups_paths_reps,
                                cpu,
                                not_included_cds,
                                constants,
                                results_output)

    print("Writting members file...")
    ccf.write_cluster_members_to_file(results_output,
                                      clusters_to_keep,
                                      clusters, frequency_in_genomes,
                                    drop_possible_loci)

    print("Create graphs for the BLAST results...")
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

    print("\nReading schema loci short FASTA files...")
    # Create directory
    results_output = os.path.join(output_directory, '4_Schema_processing')
    ff.create_directory(results_output)

    allele_ids = [True, True]
    run_type = 'loci_vs_cds' # Set run type as loci_vs_cds
    # Run Blasts for the found loci against schema short
    representative_blast_results = sl.process_schema(schema,
                                                  groups_paths,
                                                  results_output,
                                                  trans_dict_cds,
                                                  alleles,
                                                  updated_frequency_in_genomes,
                                                  allelecall_directory, 
                                                  master_file,
                                                  allele_ids,
                                                  run_type,
                                                  False,
                                                  constants,
                                                  cpu)

def main(schema, output_directory, allelecall_directory, alignment_ratio_threshold_gene_fusions, 
        pident_threshold_gene_fusions, clustering_sim, clustering_cov, genome_presence,
        size_threshold, translation_table, bsr, problematic_proportion, size_ratio, cpu):
    
    temp_paths = [os.path.join(allelecall_directory, "temp"), 
                      os.path.join(allelecall_directory, "unclassified_sequences.fasta"),
                      os.path.join(allelecall_directory, "missing_classes.fasta")]
    # Put all constants in one dict in order to decrease number of variables
    # used around.
    constants = [alignment_ratio_threshold_gene_fusions, 
                pident_threshold_gene_fusions,
                genome_presence,
                clustering_sim,
                clustering_cov,
                size_threshold,
                translation_table,
                bsr,
                problematic_proportion,
                size_ratio]
    
    if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
        sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                    "was run using --no-cleanup and --output-unclassified flag.")

    unclassified_cds_output = os.path.join(output_directory, "unclassified_cds")
    classify_cds(schema,
                 unclassified_cds_output,
                 allelecall_directory,
                constants,
                temp_paths,
                cpu)