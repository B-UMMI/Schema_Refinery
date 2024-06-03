import os

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       core_functions as cof,
                       blast_functions as bf,
                       kmers_functions as kf,
                       iterable_functions as itf,
                       graphical_functions as gf,
                       pandas_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      core_functions as cof,  
                                      blast_functions as bf,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf,
                                      pandas_functions as pf)

def create_graphs(file_path, output_path, filename, other_plots = None):
    """
    Create graphs based on representative_blast_results written inside a TSV file,
    this function creates severall plots related to palign and protein values, with
    the option to create additional plots based on inputs values.
    
    Parameters
    ----------
    file_path : str
        Path to the TSV file.
    output_path : str
        Path to the output directory.
    other_plots : list, optional
        List that contains additional data to create plots.

    Returns
    -------
    Create an HTML file inside the output_path that contains all of the created
    graphs.
    """
    results_output = os.path.join(output_path, "Graph_folder")
    ff.create_directory(results_output)
    
    blast_results_df = ff.import_df_from_file(file_path, '\t')
    
    # Create boxplots
    traces = []
    for column in ['Global_palign_all_min', 'Global_palign_all_max', 'Global_palign_pident_min', 'Global_palign_pident_max', 'Palign_local_min']:
        traces.append(gf.create_violin_plot(y = blast_results_df[column], name = blast_results_df[column].name))
    
    violinplot1 = gf.generate_plot(traces, "Palign Values between BLAST results", "Column", "Palign")
    
    # Create line plot.
    traces = []
    for column in ['Prot_BSR', 'Prot_seq_Kmer_sim', 'Prot_seq_Kmer_cov']:
        traces.append(gf.create_violin_plot(y = blast_results_df[column], name = blast_results_df[column].name))
    
    violinplot2 = gf.generate_plot(traces, "Protein values between BLAST results", "BLAST entries ID", "Columns")
    
    # Create other plots
    extra_plot = []
    if other_plots:
        for plot in other_plots:
            plot_df = pf.dict_to_df(plot[0])
            for column in plot_df.columns.tolist():
                if plot[1] == 'histogram':
                   trace = gf.create_histogram(x = plot_df[column], name = plot_df[column].name)
            
            extra_plot.append(gf.generate_plot(trace, plot[2], plot[3], plot[4]))

    gf.save_plots_to_html([violinplot1, violinplot2] + extra_plot, results_output, filename)

def main(schema, output_directory, allelecall_directory, constants, temp_paths, cpu):

    temp_folder = temp_paths[0]
    file_path_cds = temp_paths[1]
    #missing_classes_fastas = temp_paths[2]

    # Verify if the dataset is small, if it is, keep minimum genomes in which
    # specific CDS cluster is present to 5 if not to 1% of the dataset size.
    if not constants[2]:
        count_genomes_path = os.path.join(temp_folder, "1_cds_prediction")
        number_of_genomes = len(os.listdir(count_genomes_path))
        if number_of_genomes <= 20:
            constants[2] = 5
        else:
            constants[2] = round(number_of_genomes * 0.01)
        
    print("Identifying CDS present in the schema...")
    cds_present = os.path.join(temp_folder,"3_cds_preprocess/cds_deduplication/distinct_cds_merged.hashtable")
    # Get dict of CDS and their sequence hashes.
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)

    print("Identifying CDS not present in the schema...")
    # Get dict with CDS ids as key and sequence as values.
    not_included_cds = sf.fetch_fasta_dict(file_path_cds, True)

    """
    print("Identifying CDS identified as missing classes...")
    missing_classes_fastas = sf.fetch_fasta_dict(missing_classes_fastas, True)
    print("Filtering missing CDS in the schema...")
    missing_classes_fastas = {itf.remove_by_regex(key.split('|')[3], '&.*'): value 
                              for key, value in missing_classes_fastas.items() 
                              if itf.regex_present(['&ASM', '&ALM', '&NIPH', '&NIPHEM'], key)}
    # Deduplicate the FASTA dict.
    missing_classes_fastas = itf.deduplicate_fasta_dict(missing_classes_fastas)

    not_included_cds.update(missing_classes_fastas)
    """
    print("Filtering missing CDS in the schema...")
    # Count CDS size
    cds_size = {}
    for key, sequence in not_included_cds.items():
        cds_size.setdefault(key, len(str(sequence)))

    total_cds = len(not_included_cds)
    print(f"\nIdentified {total_cds} valid CDS not present in the schema.")
    # Filter by size.
    if constants[5]:
        for key, values in list(not_included_cds.items()):
            if len(values) < constants[5]:
                del not_included_cds[key]
        print(f"{len(not_included_cds)}/{total_cds} have size greater or equal to {constants[5]} bp.")
    else:
        constants[5] = 0
        print("No size threshold was applied to the CDS filtering.")

    # Create directories.
    ff.create_directory(output_directory)

    cds_output = os.path.join(output_directory, "1_CDS_processing")
    ff.create_directory(cds_output)
    # This file contains unique CDS.
    cds_not_present_file_path = os.path.join(cds_output, "CDS_not_found.fasta")
    
    # Count the number of CDS present in the schema and write CDS sequence
    # into a FASTA file.
    frequency_cds = {}
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id_+"\n")
            cds_not_found.writelines(str(sequence)+"\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus 1 (subtract the first CDS genome).
            if hashed_seq in decoded_sequences_ids:
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq]) - 1
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
                                                                           True)
    # Count translation sizes.
    cds_translation_size = {}
    for key, sequence in cds_translation_dict.items():
        cds_translation_size.setdefault(key, len(sequence))

    # Print additional information about translations and deduplications.
    print(f"\n{len(cds_translation_dict)}/{len(not_included_cds)} unique protein translations.")

    print("Extracting minimizers for the translated sequences and clustering...")
    # Create variables to store clustering info.
    reps_groups = {}
    clusters = {}
    reps_sequences = {}

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
                                                           True)
    # Print additional information about clustering.
    total_number_clusters = len(clusters)
    print(f"{len(cds_translation_dict)} unique proteins have been clustered into {total_number_clusters} clusters.")
    singleton_cluster = len([cluster for cluster in clusters if len(cluster) == 1])
    print(f"\tOut of those clusters, {singleton_cluster} are singletons")
    print(f"\tOut of those clusters, {total_number_clusters - singleton_cluster} have more than one CDS.")
    
    # Reformat the clusters output, we are interested only in  the ID of cluster members.
    clusters = {cluster_rep: [value[0] for value in values]
                for cluster_rep, values in clusters.items()}
    # For protein hashes get only those that have more than one CDS.
    protein_hashes = {hash_prot: cds_ids for hash_prot, cds_ids in protein_hashes.items()
                      if len(cds_ids) > 1}
    
    # Add also the unique CDS ID that have the same protein as representative.
    for cluster_rep, values in clusters.items():
        for cds_ids in protein_hashes.values():
            # Break since there is only one possible match in protein_hashes.
            if cluster_rep in cds_ids:
                clusters[cluster_rep] + cds_ids[1:]
                break

    print("\nFiltering clusters...")
    # Get frequency of cluster.
    frequency_cds_cluster = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_cds_cluster[rep] >= constants[2]}
    print(f"After filtering by CDS frequency in the genomes (>= {constants[2]}),"
          f" out of {total_number_clusters} clusters, {len(clusters)} remained.")

    print("\nRetrieving kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    # Get the representatives protein sequence.
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items()
                             if rep_id in clusters}
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
                                                True)
                        )
        
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep,
                                                               reps_groups,
                                                               0,
                                                               0,
                                                               prot_len_dict,
                                                               cluster_id,
                                                               5)

        reps_kmers_sim[cluster_id] = {match_values[0]: match_values[1:]
                                      for match_values in reps_kmers_sim[cluster_id]}

    # Create directories.
    blast_output = os.path.join(output_directory, "2_BLAST_processing")
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # Create directory and files path where to write FASTAs.
    representatives_blastn_folder = os.path.join(blastn_output,
                                                "cluster_representatives_fastas")
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  "all_cluster_representatives.fasta")
    # Write files for BLASTn.
    rep_paths_nuc = {}
    # Master file.
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for cluster_rep_id in clusters:

            all_fasta.writelines(">"+cluster_rep_id+"\n")
            all_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")

            rep_fasta_file = os.path.join(representatives_blastn_folder,
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths_nuc[cluster_rep_id] = rep_fasta_file
            # Representative file
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.writelines(">"+cluster_rep_id+"\n")
                rep_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")
    
    # Create BLAST db for the schema DNA sequences.
    print("\nCreating BLASTn database for the unclassified and missed CDSs...")
    blast_db = os.path.join(blastn_output, "blast_db_nuc")
    bf.make_blast_db(representatives_all_fasta_file, blast_db, 'nucl')

    # Run the BLASTn and BLASTp
    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = cof.run_blasts(blast_db,
                        clusters,
                        reps_translation_dict,
                        rep_paths_nuc,
                        blast_output,
                        constants,
                        cpu)
    
    # Add various results to the dict
    cof.add_items_to_results(representative_blast_results,
                         reps_kmers_sim,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster,
                         [False, False])

    print("\nFiltering BLAST results into classes...")
    results_output = os.path.join(output_directory, "3_Classes_processing")
    ff.create_directory(results_output)
    report_file_path = os.path.join(results_output, "blast_all_matches.tsv")
    
    # Separate results into different classes.
    classes_outcome = cof.separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    # Write all of the BLASTn results to a file.
    cof.alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("Processing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    [cds_to_keep, important_relationships, drop_list, _, _] = cof.process_classes(representative_blast_results,
                                                                                classes_outcome)

    cof.report_main_relationships(important_relationships,
                              representative_blast_results,
                              cds_to_keep['1a'],
                              False,
                              results_output)
    
    print("\nAdd remaining cluster that didn't match by BLASTn...")
    # Add cluster not matched by BLASTn
    cds_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in clusters.keys() if cluster not in representative_blast_results.keys()])

    print("\nWritting classes results to files...")
    cof.write_processed_results_to_file(cds_to_keep,
                                    representative_blast_results,
                                    classes_outcome,
                                    None,
                                    None,
                                    None,
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    [groups_paths,
     reps_trans_dict_cds,
     master_file_rep] = cof.wrap_up_blast_results(cds_to_keep,
                                              not_included_cds,
                                              clusters,
                                              results_output,
                                              constants,
                                              drop_list,
                                              None,
                                              None,
                                              frequency_cds_cluster,
                                              False)
    
    # Add new frequencies in genomes for joined groups
    new_cluster_freq = {}
    for cluster_id, cluster_members in cds_to_keep['1a'].items():
        new_cluster_freq[cluster_id] = 0
        for member in cluster_members:
            new_cluster_freq[cluster_id] += frequency_cds_cluster[member]
        for member in cluster_members:
            frequency_cds_cluster[member] = new_cluster_freq[cluster_id]

    print("Create graphs for the BLAST results...")
    cds_size_dicts = {'IDs': cds_size.keys(),
                      'Size': cds_size.values()}
    cds_translation_size_dicts = {'IDs': cds_size.keys(),
                                  'Size': [int(cds/3) for cds in cds_size.values()]}
    create_graphs(report_file_path,
                  results_output,
                  'All_of_CDS_graphs',
                  [[cds_size_dicts, 'histogram', "Nucleotide Size", 'Size', 'CDS'],
                   [cds_translation_size_dicts, 'histogram','Protein Size' , 'Size', 'CDS']])
    
    for file in ff.get_paths_in_directory(os.path.join(results_output, 'blast_results_by_class')):
        create_graphs(file,
                      results_output,
                      f"graphs_class_{os.path.basename(file).split('_')[-1].replace('.tsv', '')}")

    print("\nReading schema loci short FASTA files...")
    # Create directory
    results_output = os.path.join(output_directory, "4_Schema_processing")
    ff.create_directory(results_output)

    loci_ids = [True, False]
    # Run Blasts for the found loci against schema short
    representative_blast_results = cof.process_schema(schema,
                                                  groups_paths,
                                                  results_output,
                                                  reps_trans_dict_cds,
                                                  cds_to_keep,
                                                  frequency_cds_cluster,
                                                  allelecall_directory, 
                                                  master_file_rep,
                                                  loci_ids,
                                                  False,
                                                  constants,
                                                  cpu)