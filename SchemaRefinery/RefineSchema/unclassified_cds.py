import os
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       core_functions as cof,
                       blast_functions as bf,
                       alignments_functions as af,
                       kmers_functions as kf,
                       iterable_functions as itf,
                       graphical_functions as gf,
                       pandas_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf,
                                      pandas_functions as pf)

def write_processed_results_to_file(cds_to_keep, representative_blast_results,
                                    classes_outcome, all_alleles, cds_matched_loci, output_path):
    """
    Write the results from processed_classes into various files.
    
    Parameters
    ----------
    cds_to_keep : dict
        Dict of the CDS to keep by each classification.
    representative_blast_results : dict
        Dict that contains BLAST results of the representatives with all of the additional
        info.
    classes_outcome : list
        List of list that contains class IDS used in the next function.
    all_alleles : dict
        Dict that contains the loci and joined group as keys and their alleles and 
        elements IDS as values.
        Can be None if no loci are involved.
    cds_matched_loci : dict
        Dict that contains with which loci alleles did CDS and joined group match.
        Can be None if no loci are involved.
    output_path : str
        Path were to write files.
        
    Returns
    -------
    No returns, writes files in output path.
    """
    def process_clusters(cds_to_keep, representative_blast_results, all_alleles, cds_matched_loci, output_path):
        """
        Process and write cluster results.

        Parameters
        ----------
        cds_to_keep : dict
            Dictionary of classes and their corresponding CDS.
        representative_blast_results : dict
            Dictionary of representative blast results.
        all_alleles : dict
            Dict that contains the IDs as key and of all alleles related to that ID as values.
        cds_matched_loci : dict
            Dictionary of CDS matched loci.
        output_path : str
            Path to the output directory.

        Returns
        -------
        add_groups_ids : bool
            True if additional group IDs are present, False otherwise.
        """
        # Loop over each class and its corresponding CDS
        for class_, cds in cds_to_keep.items():
            # Loop over each cluster in the CDS
            for id_, cluster in enumerate(cds, 1):
                # Process the cluster and get the necessary details
                id_, cluster, cluster_type, is_cds, add_groups_ids = process_cluster(class_, id_, cluster, all_alleles, cds)
                # Generate a dictionary to be written to the file
                write_dict = generate_write_dict(id_, cluster, is_cds, cds_matched_loci, representative_blast_results)
                # Define the path of the report file
                report_file_path = os.path.join(output_path, f"blast_{cluster_type}_{id_}.tsv")
                # Write the dictionary to the file
                alignment_dict_to_file(write_dict, report_file_path, 'w', add_groups_ids)
        
        return add_groups_ids

    def process_cluster(class_,id_ , cluster, all_alleles, cds):
        """
        Process a single cluster.

        Parameters
        ----------
        class_ : str
            Class of the cluster.
        id_ : str or int
            ID of the cluster.
        cluster : str or int
            ID of the cluster.
        all_alleles : dict
            Dictionary of all alleles.
        cds : dict or str
            If single CDS then contain str if joined cluster then a dict.

        Returns
        -------
        id_ : str or int
            ID of the cluster.
        cluster : list
            List of the clusters.
        cluster_type : str
            Type of the cluster.
        is_cds : bool
            True if it's a CDS and not a loci, False otherwise.
        add_groups_ids : bool
            True if additional group IDs are present, False otherwise.

        """
        # Check the class and process accordingly
        if class_ == '1a':
            cluster_type = 'joined_cluster'
            cluster = cds[id_]
        else:
            id_ = cluster
            cluster = [cluster]
            cluster_type = 'retained'

        # Check if all_alleles exist
        if all_alleles:
            add_groups_ids = True
            is_cds = False
            cluster_alleles = []
            for entry in cluster:
                if entry not in all_alleles or type(entry) == int:
                    if type(entry) == int:
                        cluster = all_alleles[entry]
                    cluster_type = 'CDS_cluster'
                    is_cds = True
                else:
                    cluster_type = 'loci'
                    cluster_alleles += all_alleles[entry]
            if not is_cds:
                cluster = cluster_alleles
        else:
            add_groups_ids = False
            is_cds = True

        return id_, cluster, cluster_type, is_cds, add_groups_ids

    def generate_write_dict(id_, cluster, is_cds, cds_matched_loci, representative_blast_results):
        """
        Generate the dictionary to be written to file.

        Parameters
        ----------
        id_ : str or int
            ID of the cluster.
        cluster : list
            List of the clusters.
        is_cds : bool
            Boolean indicating if it's a CDS.
        cds_matched_loci : dict
            Dictionary of CDS matched loci.
        representative_blast_results : dict
            Dictionary of representative blast results.

        Returns
        -------
        write_dict : dict
            Dictionary to be written to file.
        """
        # Check if it's a CDS and if it matches loci
        if is_cds and cds_matched_loci:
            queries = []
            if type(id_) == int:
                queries = cds_matched_loci[id_]
            else:
                for c in cluster:
                    queries += cds_matched_loci[c]
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items() if subject in cluster}
                        for query, subjects in representative_blast_results.items()
                        if query in queries}
        else:
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items()}
                        for query, subjects in representative_blast_results.items()
                        if query in cluster}
        return write_dict

    def process_classes(classes_outcome, representative_blast_results, output_path, add_groups_ids):
        """
        Process and write class results.

        Parameters
        ----------
        classes_outcome : list
            List of class outcomes.
        representative_blast_results : dict
            Dictionary of representative blast results.
        output_path : str
            Path to the output directory.
        add_groups_ids : bool
            Boolean indicating if additional group IDs are present.

        Returns
        -------
        No returns, writes files in output path.
        """
        # Loop over each class in the outcome
        for class_ in classes_outcome:
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                                for subject, entries in subjects.items()}
                        for query, subjects in representative_blast_results.items()}
            # Define the path of the report file
            report_file_path = os.path.join(output_path, f"blastn_group_{class_}.tsv")
            # Write the dictionary to the file
            cof.alignment_dict_to_file(write_dict, report_file_path, 'w', add_groups_ids)

    # Create directories for output
    blast_by_cluster_output = os.path.join(output_path, 'blast_by_cluster')
    ff.create_directory(blast_by_cluster_output)
    blast_results_by_class_output = os.path.join(output_path, 'blast_results_by_class')
    ff.create_directory(blast_results_by_class_output)

    # Process and write cluster results
    add_groups_ids = process_clusters(cds_to_keep, representative_blast_results, all_alleles, cds_matched_loci, blast_by_cluster_output)

    # Process and write class results
    process_classes(classes_outcome, representative_blast_results, blast_results_by_class_output, add_groups_ids)

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

def process_schema(schema, groups_paths, results_output, reps_trans_dict_cds, 
                   cds_to_keep, frequency_cds_cluster, allelecall_directory, 
                   master_file_rep, not_included_cds, constants, cpu):
    """
    This function processes data related to the schema seed, importing, translating
    and BLASTing against the unclassified CDS clusters representatives groups to
    validate them.
    
    Parameters
    ----------
    schema : str
        Path to the schema seed folder.
    groups_trans_reps : dict
        Dict with the paths to the translations of the unclassified CDS clusters.
    results_output : str
        Path were to write the results of this function.
    reps_trans_dict_cds : dict
        Dict that contains the translations for each CDS.
    cds_to_keep : dict     
        Dict of the CDS to keep by each classification.
    master_file_rep : str
        Path to the maste file containing retained CDS.
    frequency_cds_cluster : dict
        Contains the frequency of each CDS/loci in the genomes.
    allelecall_directory : str
        Path to the allele call directory.
    not_included_cds : dict
        Dict that contains the unclassified CDS IDs as keys and their
        DNA sequences as values.
    constants : list
        Contains the constants to be used in this function.
    cpu : int
        Number of CPUs to use during multi processing.

    Returns
    -------
    representative_blast_results : dict
        Dict that contains BLAST results of the representatives with all of the additional
        info.

    """
    # Get all of the schema loci short FASTA files path.
    schema_short_path = os.path.join(schema, 'short')
    schema_loci_short = {loci_path.replace("_short.fasta", ""): os.path.join(schema_short_path, loci_path) 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    # Create a folder for short translations.
    short_translation_folder = os.path.join(results_output, "short_translation_folder")
    ff.create_directory(short_translation_folder)

    # Find the file in the allele call results that contains the total of each.
    # classification obtained for each loci.
    results_statistics = os.path.join(allelecall_directory, 'loci_summary_stats.tsv')
    # Convert TSV table to dict.
    results_statistics_dict = itf.tsv_to_dict(results_statistics)
    # Add the results for all of the Exact matches to the frequency_cds_cluster dict.
    for key, value in results_statistics_dict.items():
        frequency_cds_cluster.setdefault(key, int(value[0]))
    # Translate each short loci and write to master fasta.
    i = 1
    len_short_folder = len(schema_loci_short)
    all_alleles = {}
    for loci, loci_short_path in schema_loci_short.items():
        print(f"\rTranslated fasta short loci: {i}/{len_short_folder}", end='', flush=True)
        i += 1
        fasta_dict = sf.fetch_fasta_dict(loci_short_path, False)
        
        for loci_id, sequence in fasta_dict.items():
            all_alleles.setdefault(loci, []).append(loci_id)

        loci_short_translation_path = os.path.join(short_translation_folder, f"{loci}.fasta")
        translation_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                              loci_short_translation_path,
                                                              None,
                                                              constants[5],
                                                              False)
        for loci_id, sequence in translation_dict.items():
            reps_trans_dict_cds[loci_id] = sequence

    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = cof.run_blasts(master_file_rep,
                     schema_loci_short,
                     reps_trans_dict_cds,
                     schema_loci_short,
                     results_output,
                     constants,
                     cpu,
                     all_alleles)

    cof.add_items_to_results(representative_blast_results,
                         None,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster,
                         True,
                         cds_to_keep['1a'])

    # Add CDS joined clusters to all_alleles IDS
    cds_joined_cluster = cds_to_keep['1a']
    all_alleles.update(cds_joined_cluster)
    # Separate results into different classes.
    classes_outcome = cof.separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    
    report_file_path = os.path.join(results_output, "blast_all_matches.tsv")
    # Write all of the BLASTn results to a file.
    cof.alignment_dict_to_file(representative_blast_results, report_file_path, 'w', True)
    
    print("\nProcessing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    [cds_to_keep, important_relationships, drop_list] = cof.process_classes(representative_blast_results,
                                                                        classes_outcome,
                                                                        all_alleles)
    # Replace the alleles entries with their loci ID.
    cds_to_keep = {
        class_: set(
            [entry if not itf.identify_string_in_dict(entry, all_alleles) else itf.identify_string_in_dict(entry, all_alleles) for entry in entries]
        )
        if class_ != '1a' else entries for class_, entries in cds_to_keep.items()
    }

    drop_list = set([entry if not itf.identify_string_in_dict(entry, all_alleles) else itf.identify_string_in_dict(entry, all_alleles) for entry in drop_list])
    # Filter repeated entries
    seen = set()
    for class_, entries in list(cds_to_keep.items()):
        for entry in list(entries):
            if entry not in seen:
                seen.add(entry)
            else:
                cds_to_keep[class_].remove(entry)

    cds_matched_loci = {}
    for class_, entries in list(cds_to_keep.items()):
        for entry in list(entries):
            if entry not in schema_loci_short:
                if type(entry) == int:
                    id_ = entry
                    entry = cds_joined_cluster[entry]
                else:
                    id_ = entry
                    entry = [entry]
                cds_matched_loci.setdefault(id_, set([i[0] for i in itf.flatten_list(important_relationships.values()) if i[1] in entry]))

    print("\nWritting classes results to files...")
    write_processed_results_to_file(cds_to_keep,
                                    representative_blast_results,
                                    classes_outcome,
                                    all_alleles,
                                    cds_matched_loci,
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    
    cof.report_main_relationships(important_relationships,
                              representative_blast_results,
                              all_alleles,
                              True,
                              results_output)

    [_, _, _] = cof.wrap_up_blast_results(cds_to_keep,
                                           not_included_cds,
                                           all_alleles,
                                           results_output,
                                           constants,
                                           drop_list,
                                           schema_loci_short,
                                           groups_paths,
                                           None)

    return representative_blast_results

def main(schema, output_directory, allelecall_directory, constants, temp_paths, cpu):

    temp_folder = temp_paths[0]
    file_path_cds = temp_paths[1]

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
                

    print("\nTranslate and deduplicate unclassified CDS...")
    # Translate the CDS and find unique proteins using hashes, the CDS with
    # the same hash will be added under that hash in protein_hashes.
    cds_not_present_trans_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    cds_not_present_untrans_file_path = os.path.join(cds_output, "CDS_not_found_untranslated.fasta")
    # Translate and deduplicate protein sequences.
    cds_translation_dict, protein_hashes, _ = sf.translate_seq_deduplicate(not_included_cds,
                                                                           cds_not_present_trans_file_path,
                                                                           cds_not_present_untrans_file_path,
                                                                           constants[5],
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
    
    # Run the BLASTn and BLASTp
    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     self_score_dict] = cof.run_blasts(representatives_all_fasta_file,
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
                         False)

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
    [cds_to_keep, important_relationships, drop_list] = cof.process_classes(representative_blast_results,
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
                                              frequency_cds_cluster)
    
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
    # Create BLASTn_processing directory
    blastn_output = os.path.join(results_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # Run Blasts for the found loci against schema short
    representative_blast_results = process_schema(schema,
                                                  groups_paths,
                                                  results_output,
                                                  reps_trans_dict_cds,
                                                  cds_to_keep,
                                                  frequency_cds_cluster,
                                                  allelecall_directory, 
                                                  master_file_rep,
                                                  not_included_cds,
                                                  constants,
                                                  cpu)