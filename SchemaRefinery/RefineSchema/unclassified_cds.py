import os
import pickle
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       blast_functions as bf,
                       alignments_functions as af,
                       kmers_functions as kf,
                       iterable_functions as itf,
                       graphical_functions as gf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf)

def fetch_fasta_dict(file_path_cds, count_seq):
    """
    Compares the hashes list with the hashes obtained from cds.

    Parameters
    ----------
    hash_list : list
        List containing all of the sequences hashes present in the input files.
    count_seq : bool
        If count the number of processed sequences inside the file

    Returns
    -------
    fasta_dict : dict
        Returns dict with key as fasta header and value as fasta sequence.
    """
    
    fasta_dict = {}
    i = 1
    # Read FASTA files
    for rec in sf.read_fasta_file_iterator(file_path_cds):
        if count_seq:
            print(f"Processed {i} CDS")
            i += 1
        fasta_dict[rec.id] = rec.seq

    return fasta_dict

def decode_CDS_sequences_ids(path_to_file):
    """
    Function to read a dict contained in pickle file and decode its values based
    on polyline.
    
    Parameters
    ----------
    path_to_file : str
        Path to the pickle file.
        
    Returns
    -------
    decoded_dict : dict
        Contains hashed CDS as keys, and number id of the genome which
        that CDS is present
    
    """
    # Load pickle file
    with open(path_to_file, "rb") as infile:
        hash_table = pickle.load(infile)
    # Decode the file
    decoded_dict = {}
    for key, value in hash_table.items():
        decoded_dict[key] = itf.polyline_decoding(value)
        
    return decoded_dict

def alignment_dict_to_file(blast_results_dict, file_path, write_type):
    """
    Writes alignments strings to file.

    Parameters
    ----------
    alignment_string_dict : dict
        dict containing alignments string to write to file
    file_path : str
        File path to create to write the file
    write_type : str
        If to create new file and write or to append to existing file

    Returns
    -------
    No return, writes or appends a file at the file_path
    """
    # Write first column into TSV file
    with open(file_path, write_type) as report_file:
        if write_type == 'w':
            report_file.writelines(["Query\t",
                                    "Subject\t",
                                    "Query_length\t",
                                    "Subject_length\t",
                                    "Query_start\t",
                                    "Query_end\t",
                                    "Subject_start\t",
                                    "Subject_end\t",
                                    "Length\t",
                                    "Score\t",
                                    "Number_of_gaps\t",
                                    "Pident\t",
                                    "Prot_BSR\t",
                                    "Prot_seq_Kmer_sim\t",
                                    "Prot_seq_Kmer_cov\t",
                                    "Cluster_frequency_in_genomes_query_cds\t",
                                    "Cluster_frequency_in_genomes_subject_cds\t",
                                    "Global_palign_all\t",
                                    "Global_palign_pident_min\t",
                                    "Global_palign_pident_max\t",
                                    "Palign_local\t",
                                    "Class\n"])
        # Write all of the matches into TSV file
        for results in blast_results_dict.values():
            for result in results.values():
                for r in result.values():
                    report_file.writelines('\t'.join([str(r) for r in r.values()]) + '\n')

def translate_seq_deduplicate(seq_dict, path_to_write, count_seq):
    """
    Translates the DNA sequence to protein and verifies if that protein is alredy
    present in the dict, thus ensuring that the dict contains deduplicated sequences,
    it writes the sequences to a FASTA files and return the dict.
    
    Parameters
    ----------
    seq_dict : dict
        Dict that contains sequence ID as key and the sequence as value.
    path_to_write : str
        Path to the file to create and write.
    count_seq : bool
        If there is need to print into stdout the number of processed sequences.
        
    Returns
    -------
    translation_dict : dict
        Dict that contais sequence ID as key and translated sequence as value.
    protein_hashes : dict
        Dict that contais sequence hash as key and sequences IDs as values.
    """
    translation_dict = {}
    protein_hashes = {}
    if count_seq:
        i = 1
        total = len(seq_dict)
    with open(path_to_write, 'w+') as translation:
        for id_s, sequence in seq_dict.items():
            if count_seq:
                print(f"Translated {i}/{total} CDS")
                i += 1
            # Translate
            protein_translation = str(sf.translate_dna(str(sequence),
                                                       11,
                                                       0,
                                                       True)[0][0])
            # Hash the sequence
            prot_hash = sf.seq_to_hash(protein_translation)
            # Find unique proteins
            if prot_hash not in protein_hashes:
                protein_hashes[prot_hash] = [id_s]
                translation_dict[id_s] = protein_translation
                translation.writelines('>'+id_s+"\n")
                translation.writelines(protein_translation+"\n")
            # Remember CDS with that protein hash for future
            else:
                protein_hashes[prot_hash].append(id_s)
                
    return translation_dict, protein_hashes

def separate_blastn_results_into_classes(representative_blast_results, constants):
    """
    Separates one BLASTn dict into various classes and adds them into one dict

    Parameters
    ----------
    representative_blast_results : dict
        dict containing BLAST results

    Returns
    -------
    results_outcome : dict
        Dict that contains the outcomes for the results when they were filtered
        by classes.
    classes_outcome : list
        List of list that contains class IDS used in the next function
    """
        
    def add_class_to_dict(class_name):
        """
        Adds class as last item in representative_blast_results dict
        
        Parameters
        ----------
        class_name : str
            Class name, which is a key in cluster_classes to where to add entries.
            
        Returns
        -------
        Modifies the representative_blast_results dict in the parent function.
        """
        
        representative_blast_results[query][id_subject][id_].update({'class': class_name})
        
    results_outcome = {}
    
    # Create all of the classes
    classes_outcome = ['1a',
                       '3a',
                       '1b',
                       '1c',
                       '2a',
                       '2b',
                       '3b',
                       'drop']

    for class_ in classes_outcome:
        results_outcome[class_] = []
    # Process results into classes
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                
                if blastn_entry['global_palign_all'] >= 0.8:
                    # Based on BSR
                    if blastn_entry['bsr'] >= 0.6:
                        add_class_to_dict('1a')
                        results_outcome['1a'].append([query, id_subject])
                    # If BSR <0.6 verify if query cluster is the most prevalent
                    elif blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds'] * 10:
                        add_class_to_dict('1b')
                        results_outcome['1b'].append([query, id_subject])
                        
                    # Ignore cases were 1b was already assigned or will be assigned
                    elif blastn_entry['frequency_in_genomes_subject_cds'] >= blastn_entry['frequency_in_genomes_query_cds'] * 10:
                        add_class_to_dict('drop')
                        results_outcome['drop'].append([query, id_subject])
                    # Add two as separate
                    else:
                        add_class_to_dict('1c')
                        results_outcome['1c'].append([query, id_subject])
                # Palign < 0.8        
                else:
                    if blastn_entry['pident'] >= constants[1]:
                        # verify if query cluster is the most prevalent
                        if blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds'] * 10:
                            add_class_to_dict('2a')
                            results_outcome['2a'].append([query, id_subject])
                        # Ignore cases were 2a was already assigned or will be assigned
                        elif blastn_entry['frequency_in_genomes_subject_cds'] >= blastn_entry['frequency_in_genomes_query_cds'] * 10:
                            add_class_to_dict('drop')
                            results_outcome['drop'].append([query, id_subject])
                        else:
                            add_class_to_dict('2b')
                            results_outcome['2b'].append([query, id_subject])
                            
                    else:
                        if blastn_entry['global_palign_pident_max'] >= 0.8:
                            add_class_to_dict('3a')
                            results_outcome['3a'].append([query, id_subject])
                        else:
                            add_class_to_dict('3b')
                            results_outcome['3b'].append([query, id_subject])

    return results_outcome, classes_outcome

def process_classes(results_outcome):
    """
    Identifies the relationships between representatives classified as other classes
    that matched by BLASTn with members join member of the same cluster classified
    as 1a, those relationships are written to a file for the end user to see if
    they are relevant.
    
    Parameters
    ----------
    results_outcome : dict
        Dict that contains the outcomes for the results when they were filtered
        by classes.
        
    Returns
    -------
    results_outcome : dict
        Dict that contains the outcomes for the results when they were filtered
        by classes.
    relationships : dict
        Dict that contains relationships between various clusters.
    cluster_dict_1a : dict
        Joined clusters dict, these clusters contain various CDS representatives.
    """

    # Join the various CDS groups into single group based on ids matches and remove Join from results_outcome
    cluster_dict_1a = {i+1: join for i, join in enumerate(cf.cluster_by_ids(results_outcome.pop('1a')))}
            
    # and remove duplicates
    for class_, results in results_outcome.items():
        results_outcome[class_] = itf.get_unique_sublists(results_outcome[class_])
        
    relationships = {}
    for class_ in results_outcome:
        relationships[class_] = {}
    
    # Process all elements in other classes against each other
    for class_, results in results_outcome.items():
        if class_ == '1a':
            continue
        # Iterate over results
        for result in results:
            # If cluster is classified as 1a add the id of the cluster
            cluster_id = itf.identify_string_in_dict(result[1], cluster_dict_1a)
            if not cluster_id:
                cluster_id = result[1]
            # If entry exists
            if cluster_id not in relationships[class_]:
                relationships[class_].update({cluster_id: [result]})
            # Add to entry
            else:
                relationships[class_][cluster_id].append(result)

    # Remove entries that were added to Joined cluster
    drop = itf.flatten_list(list(cluster_dict_1a.values()))
    for class_, results in results_outcome.items():
        if class_ in ['2a','2b','1c','3a','3b']:
            results_outcome[class_] = [result for result in results
                                       if result[0] not in drop]
    # Remove duplicates entries
    for class_, results in results_outcome.items():
        results_outcome[class_] = set([result[0] for result in results])

    # Get unique relationships
    for class_, relationship in list(relationships.items()):
        for id_, r in list(relationship.items()):
            relationships[class_][id_] = itf.get_unique_sublists(r)
            
    # Add the joined cluster again to the dict
    results_outcome['1a'] = list(cluster_dict_1a.values())
    return results_outcome, relationships, cluster_dict_1a

def write_processed_results_to_file(results_outcome, relationships, representative_blast_results,
                                    cluster_dict_1a, classes_outcome, output_path):
    """
    Write the results from processed_classes into various files.
    
    Parameters
    ----------
    results_outcome : dict
        Dict that contains the outcomes for the results when they were filtered
        by classes.
    relationships : dict
        Dict that contains relationships between various clusters
    representative_blast_results : dict
        Dict containg representatibes BLAST results with all of the additional
        info.
    cluster_dict_1a : dict
        Joined clusters dict, these clusters contain various CDS representatives.
    classes_outcome : list
        List of list that contains class IDS used in the next function   
    output_path : str
        Path were to write files.
        
    Returns
    -------
    No returns, writes files in output path.
    """
    
    # Create directory
    blast_by_cluster_output = os.path.join(output_path, "blast_results_by_cluster")
    ff.create_directory(blast_by_cluster_output)
    # Get all of the BLAST entries for that cluster
    for class_, cluster in results_outcome.items():
        for i, cluster in enumerate(cluster):
            if class_ == '1a':
                id_ = i + 1
                cluster_type = 'joined_cluster'
            else:
                id_ = cluster
                cluster = [cluster]
                cluster_type = 'retained'
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                   for subject, entries in subjects.items()}
                          for query, subjects in representative_blast_results.items()
                          if query in cluster}
        
            report_file_path = os.path.join(blast_by_cluster_output, f"blast_{cluster_type}_{id_}.tsv")
            alignment_dict_to_file(write_dict, report_file_path, 'w')
    
    # Create directory 
    joined_cluster_relationships_output = os.path.join(blast_by_cluster_output, "1_blast_results_by_cluster_relationships")
    ff.create_directory(joined_cluster_relationships_output)
    report_relationships_output = os.path.join(blast_by_cluster_output, "2_relationships_to_joined_clusters")
    ff.create_directory(report_relationships_output)
    
    for class_, relationship in relationships.items():
        if len(relationship) == 0:
            continue
        for cluster_id, r_ids in relationship.items():
            # Search the entries if query is in relationship_ids set and if the
            # subject is member of the joined representatives cluster (we want
            # only relevant BLAST matches) (substract -1) to adjust for cluster
            # number id and to fetch from list and lastly get the right class,
            # so it is ordered by class.
            
            query_ids = [query_id[0] for query_id in r_ids]
            subject_ids = [subject_id[1] for subject_id in r_ids]
            
            # create the lists
            transformed_query_ids = []
            # Verify if the ids are present inside joined clusters
            # if they are then add the joined cluster id to the transformed_query_ids list
            if itf.any_match_lists(query_ids, itf.flatten_list(list(cluster_dict_1a.values()))):
                # Iterate over current ids
                for id_ in query_ids:
                    # If id_ is present in some joined cluster this function returns
                    # the key otherwise it returns None
                    value = itf.identify_string_in_dict(id_, cluster_dict_1a)
                    if value:
                        transformed_query_ids.append(value)
                    # If it is None then just add the original id
                    else:
                        transformed_query_ids.append(id_)
            # If no elements of query_ids matched with any of the joined clusters
            else:
                transformed_query_ids = query_ids
            
            # Get entries based on BLAST
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()
                                             if entry['class'] == class_}
                                   for subject, entries in subjects.items() if subject in subject_ids}
                          for query, subjects in representative_blast_results.items() if query in query_ids}
                
            report_file_path = os.path.join(joined_cluster_relationships_output, f"blast_relationships_to_{cluster_id}.tsv")
            if os.path.exists(report_file_path):
                write_type = 'a'
            else:
                write_type = 'w'
            # Write BLAST results to file
            alignment_dict_to_file(write_dict, report_file_path, write_type)
                
            relationships_report_file_path = os.path.join(report_relationships_output, f"relationships_to_cluster_{cluster_id}_report.txt")
            
            # Write all of the report files
            with open(relationships_report_file_path, write_type) as relationships_report_file:
                # Class that has CDS that are partially contained or cantains other CDS
                if class_ == '3a':

                    relationships_report_file.writelines("The following CDS may partially contain the following "
                                                         "CDS from this cluster:\n")
                else:
                    relationships_report_file.writelines("The following CDS matched with BLASTn to the following"
                                                         f" elements of this cluster and have the classification '{class_}'"
                                                         " to this elements however they are probably different loci:\n")
                # Write the cluster id
                relationships_report_file.writelines(f"{cluster_id}:\n")
                seen = ""
                for i, query_id in enumerate(transformed_query_ids):
                    # If query_id is a string value
                    if type(query_id) == str:
                        # If cluster has a string value, meaning that it only contains one CDS
                        # it doesn´t matter to add subjects since they will all be the same
                        if type(cluster_id) == str:
                            relationships_report_file.writelines(f"\t{query_id}\n")
                        # If query was already seen, add white spaces to add another subject under it to facilitate readibility
                        elif query_ids[i] == seen:
                            white_spaces = itf.create_whitespace_string(f"CDS {query_ids[i]} against ")
                            relationships_report_file.writelines("\t" + white_spaces + f"{subject_ids[i]}\n")
                        # Add subjects to the report to know to what CDS in the cluster did the query match
                        else:
                            relationships_report_file.writelines(f"\tCDS {query_ids[i]} against {subject_ids[i]}\n")
                    # If query was already seen, add white spaces to add another subject under it to facilitate readibility
                    elif query_ids[i] == seen:
                        white_spaces = itf.create_whitespace_string(f"Cluster {query_id} entry: {query_ids[i]} against ")
                        relationships_report_file.writelines("\t" + white_spaces + f"{subject_ids[i]}\n")
                    # Add subjects to the report to know to what CDS in the cluster did the query match
                    else:
                        relationships_report_file.writelines(f"\tCluster {query_id} entry: {query_ids[i]} against {subject_ids[i]}\n")
                    seen = query_ids[i]

    # Write all of the ids inside Joined cluster
    # Create directory 
    cluster_members_output = os.path.join(output_path, "joined_cluster_members")
    ff.create_directory(cluster_members_output)
    # Write files
    for cluster_id, cluster in cluster_dict_1a.items():
        cluster_output_path = os.path.join(cluster_members_output, f"Joined_cluster_{cluster_id}.txt")
        with open(cluster_output_path, 'w') as output:
            for c in cluster:
                output.writelines(c + '\n')

    # Create directory 
    joined_cluster_relationships_output = os.path.join(output_path, "blast_results_by_class")
    ff.create_directory(joined_cluster_relationships_output)
    # Write classes to file
    for class_ in classes_outcome:
        # Fetch all entries with the desired class
        write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                               for subject, entries in subjects.items()}
                      for query, subjects in representative_blast_results.items()}

        report_file_path = os.path.join(joined_cluster_relationships_output, f"blastn_group_{class_}.tsv")
        # Write individual class to file
        alignment_dict_to_file(write_dict, report_file_path, 'w')

def wrap_up_results(schema, results_outcome, not_included_cds, clusters,
                    blastn_processed_results_path, output_path, cpu):
    """
    This function wraps up the results for this module by writing FASTAs files
    for the possible new loci to include into the schema and creates graphs for
    each results group. It also translates schema short FASTAs into proteins
    and the possible new loci, to calculate the BSR values to see if those
    possible new loci are already present in the schema.
    
    Parameters
    ----------
    results_outcome : dict
        Dict that contains the outcomes for the results when they were filtered
        by classes.
    not_included_cds : dict
        Dict that contains all of the DNA sequences for all of the CDS.
    clusters : dict
        Dict that contains the cluster representatives as keys and similar CDS
        as values.
    blastn_processed_results_path : str
        Path to the folder where classes TSV files were saved.
    output_path : str
        Path to were write the FASTA files.
    cpu : int
        Number of cores to use in BLAST multiprocessing
        
    Returns
    -------
    Writes TSV and HTML files
    """
    # Create directories
    cds_outcome_results = os.path.join(output_path, "Blast_results_outcomes_graphs")
    ff.create_directory(cds_outcome_results)
    
    cds_outcome_results_fastas_folder = os.path.join(output_path, "results_outcomes_fastas")
    ff.create_directory(cds_outcome_results_fastas_folder)
    
    # Write FASTA files for each CDS group to join or retain
    print("Writting FASTA file for possible new loci...")
    outcome_paths = {}
    for outcome in results_outcome:
        i = 1
        for group in results_outcome[outcome]:
            if outcome == '1a':
                cds_outcome_results_fastas_file = os.path.join(cds_outcome_results_fastas_folder, f"Joined_outcome_{i}.fasta")
                outcome_paths[f"Joined_{outcome}_{i}"] = cds_outcome_results_fastas_file
            else:
                cds_outcome_results_fastas_file = os.path.join(cds_outcome_results_fastas_folder, f"Retained_outcome_{group}.fasta")
                outcome_paths[f"Retained_{outcome}_{group}"] = cds_outcome_results_fastas_file
                group = [group]
            i += 1
            with open(cds_outcome_results_fastas_file, 'w') as fasta_file:
                for rep_id in group:
                    cds_ids = [cds_id for cds_id in clusters[rep_id]]
                    for cds_id in cds_ids:
                        fasta_file.writelines(f">{cds_id}\n")
                        fasta_file.writelines(str(not_included_cds[cds_id])+"\n")
    # Create directories
    cds_outcome_translation = os.path.join(output_path, "cds_outcome_translation")
    ff.create_directory(cds_outcome_translation)
    # Translate possible new loci
    outcomes_translations = {}
    for key, o_path in outcome_paths.items():
        translation_path = os.path.join(cds_outcome_translation, key + ".fasta")
        outcomes_translations[key] = translation_path
        fasta_dict = fetch_fasta_dict(o_path, False)
        translation_dict, _ = translate_seq_deduplicate(fasta_dict, translation_path, False)
        
    print("Reading schema loci short FASTA files...")
    # Get all of the schema loci short FASTA files path
    schema_short_path = os.path.join(schema, 'short')
    schema_loci_short = {loci_path.replace(".fasta", ""): os.path.join(schema_short_path, loci_path) 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    
    short_translation_folder = os.path.join(output_path, "short_translation_folder")
    ff.create_directory(short_translation_folder)
    master_loci_short_translation_path = os.path.join(short_translation_folder, "master_short_loci.fasta")
    i = 1
    len_short_folder = len(schema_loci_short)
    with open(master_loci_short_translation_path, 'w') as master_fasta:
        
        for loci, loci_short_path in schema_loci_short.items():
            print(f"Translated {i}/{len_short_folder} CDS")
            loci_short_translation_path = os.path.join(short_translation_folder, f"{loci}.fasta")
            i += 1
            fasta_dict = fetch_fasta_dict(loci_short_path, False)
            translation_dict, _ = translate_seq_deduplicate(fasta_dict, loci_short_translation_path, False)
            
            for loci_id, sequence in translation_dict.items():
                master_fasta.writelines(">"+loci_id+"\n")
                master_fasta.writelines(str(sequence)+"\n")
    
    # Run BLASTp between new possible loci vs existing loci (all new loci sequences vs loci short FASTA)
    blastp_results_vs_loci_results = os.path.join(output_path, "blastp_results_vs_loci_results")
    ff.create_directory(blastp_results_vs_loci_results)
    bsr_values = {}
    # Create query entries
    for query in outcomes_translations:
        bsr_values[query] = {}
    i = 1
    total = len(outcomes_translations)
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                outcomes_translations, 
                                repeat('blastp'),
                                repeat(blastp_results_vs_loci_results),
                                repeat(outcomes_translations),
                                repeat(master_loci_short_translation_path)):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, False)
            # Since BLAST may find several local aligments choose the largest one to calculate BSR
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    largest_alignment = 0
                    for entry_id, result in results.items():
                        if result['query_length'] > largest_alignment:
                            largest_alignment = result['query_length']
                            bsr_values[query].update({subject_id: bf.compute_bsr(result['score'], self_score)})
                            
            print(
                f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total}")
            i += 1             
    # Create graphs for all results and for each class
    for tsv_file_path in os.listdir(blastn_processed_results_path):
        abs_path = os.path.join(blastn_processed_results_path, tsv_file_path)
        file_name = os.path.basename(tsv_file_path)
        # Create directories
        graphs_path = os.path.join(cds_outcome_results, f"{file_name.replace('tsv','')}")
        ff.create_directory(graphs_path)
        # Render histograms
        gf.render_histogram(abs_path,
                            graphs_path,
                            ['Query_length', 'Subject_length'],
                            ['Length', 'Count'])
        # Render line charts
        gf.render_line_chart(abs_path,
                             graphs_path,
                             ['Pident', 'Prot_BSR', 'Prot_seq_Kmer_sim',
                              'Prot_seq_Kmer_cov'],
                             ['Entries', 'Values'], False)

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
    # Get dict of CDS and their sequence hashes
    decoded_sequences_ids = decode_CDS_sequences_ids(cds_present)

    print("Identifying CDS not present in the schema...")
    # Get dict with CDS ids as key and sequence as values
    not_included_cds = fetch_fasta_dict(file_path_cds, True)

    print(f"Identified {len(not_included_cds)} valid CDS not present in the schema")
    # Create directories
    ff.create_directory(output_directory)

    cds_output = os.path.join(output_directory, "1_CDS_processing")
    ff.create_directory(cds_output)
    # This file contains unique CDS
    cds_not_present_file_path = os.path.join(cds_output, "CDS_not_found.fasta")
    
    # Count the number of CDS present in the schema and write CDS sequence
    # into a FASTA file
    frequency_cds = {}
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id_+"\n")
            cds_not_found.writelines(str(sequence)+"\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus 1 (subtract the first CDS genome)
            if hashed_seq in decoded_sequences_ids:
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq]) - 1
            else:
                frequency_cds[id_] = 0
                

    print("Translate unclassified CDS...")
    # Translate the CDS and find unique proteins using hashes, the CDS with
    # the same hash will be added under that hash in protein_hashes
    cds_not_present_translation_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    # Translate and deduplicate protein sequences
    cds_translation_dict, protein_hashes = translate_seq_deduplicate(not_included_cds,
                                                                     cds_not_present_translation_file_path,
                                                                     True)

    print("Extracting minimizers for the translated sequences and clustering...")
    # Create variables to store clustering info
    reps_groups = {}
    clusters = {}
    reps_sequences = {}

    # Sort by size of proteins
    cds_translation_dict = {k: v for k, v in sorted(cds_translation_dict.items(),
                                                    key=lambda x: len(x[1]),
                                                    reverse=True)}
    # Cluster by minimizers
    [clusters, reps_sequences, 
     reps_groups, prot_len_dict] = cf.minimizer_clustering(cds_translation_dict,
                                                           5, 5, True, 1, 
                                                           clusters,
                                                           reps_sequences, 
                                                           reps_groups,
                                                           1, constants[3], 
                                                           constants[4],
                                                           True)
    # Reformat the clusters output, we are interested only in  the ID of cluster members
    clusters = {cluster_rep: [value[0] for value in values]
                for cluster_rep, values in clusters.items()}
    # For protein hashes get only those that have more than one CDS
    protein_hashes = {hash_prot: cds_ids for hash_prot, cds_ids in protein_hashes.items()
                      if len(cds_ids) > 1}
    
    # Add also the unique CDS ID that have the same protein as representative
    for cluster_rep, values in clusters.items():
        for cds_ids in protein_hashes.values():
            # Break since there is only one possible match in protein_hashes
            if cluster_rep in cds_ids:
                clusters[cluster_rep] + cds_ids[1:]
                break

    print("Filtering clusters...")
    # Get frequency of cluster
    frequency_cds_cluster = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_cds_cluster[rep] >= constants[2]}

    print("Retrieving kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    # Get the representatives protein sequence
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items()
                             if rep_id in clusters}
    # Sort the representative translation dict from largest to smallest
    reps_translation_dict = {k: v for k, v in sorted(reps_translation_dict.items(),
                                                     key=lambda x: len(x[1]),
                                                     reverse=True)}
    # recalculate the sim and cov between reps, get all of the values, so threshold
    # is set to 0.
    for cluster_id in reps_translation_dict:
        kmers_rep = set(kf.determine_minimizers(reps_translation_dict[cluster_id],
                                                5, 5, 1, True, True))
        
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep,
                                                               reps_groups,
                                                               0,
                                                               0,
                                                               prot_len_dict,
                                                               cluster_id, 5)

        reps_kmers_sim[cluster_id] = {match_values[0]: match_values[1:]
                                      for match_values in reps_kmers_sim[cluster_id]}

    print("Running BLASTn between cluster representatives...")
    # Create directories
    blast_output = os.path.join(output_directory, "2_BLAST_processing")
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # Create directory and files path where to write FASTAs
    representatives_blastn_folder = os.path.join(blastn_output,
                                                "cluster_representatives_fastas")
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  "all_cluster_representatives.fasta")
    # Write files for BLASTn
    rep_paths_nuc = {}
    # Master file
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
    # Create directory
    blastn_results_folder = os.path.join(blastn_output, "blastn_results")
    ff.create_directory(blastn_results_folder)
    # Run BLASTn for all representatives (rep vs all)
    total_reps = len(rep_paths_nuc)
    representative_blast_results = {}
    representative_blast_results_coords_all = {}
    representative_blast_results_coords_pident = {}
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                clusters,
                                repeat('blastn'),
                                repeat(blastn_results_folder),
                                repeat(rep_paths_nuc),
                                repeat(representatives_all_fasta_file)):

            filtered_alignments_dict, _, alignment_coords_all, alignment_coords_pident = af.get_alignments_dict_from_blast_results(
                res[1], constants[1], True, False)
            # Save the BLASTn results
            representative_blast_results.update(filtered_alignments_dict)
            representative_blast_results_coords_all.update(alignment_coords_all)
            representative_blast_results_coords_pident.update(alignment_coords_pident)

            print(
                f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps}")
            i += 1

    print("Running BLASTp based on BLASTn results matches...")
    # Obtain the list for what BLASTp runs to do, no need to do all vs all as previously
    # Based on BLASTn results
    blastp_runs_to_do = {query: itf.flatten_list([[query],[subject[1]['subject']
                                            for subject in subjects.values()]]) 
                         for query, subjects in representative_blast_results.items()}
    
    # Create directories
    blastp_results = os.path.join(blast_output,
                                  "BLASTp_processing")
    ff.create_directory(blastp_results)
    
    blastn_results_matches_translations = os.path.join(blastp_results,
                                                       "blastn_results_matches_translations")
    ff.create_directory(blastn_results_matches_translations)

    representatives_blastp_folder = os.path.join(blastn_results_matches_translations,
                                                "cluster_rep_translation")
    ff.create_directory(representatives_blastp_folder)
    
    blastp_results_folder = os.path.join(blastp_results,
                                         "BLASTp_results")
    ff.create_directory(blastp_results_folder)
    # Write the protein FASTA files
    rep_paths_prot = {}
    rep_matches_prot = {}    
    for query_id, subjects_ids in blastp_runs_to_do.items():
        # First write the representative protein sequence
        rep_translation_file = os.path.join(representatives_blastp_folder,
                                            f"cluster_rep_translation_{query_id}.fasta")
        rep_paths_prot[query_id] = rep_translation_file
        with open(rep_translation_file, 'w') as trans_fasta:
            trans_fasta.writelines(">"+query_id+"\n")
            trans_fasta.writelines(str(reps_translation_dict[query_id])+"\n")
        # Then write in another file all of the matches for that protein sequence
        # including the representative itself
        rep_matches_translation_file = os.path.join(blastn_results_matches_translations,
                                                    f"cluster_matches_translation_{query_id}.fasta")
        
        rep_matches_prot[query_id] = rep_matches_translation_file
        with open(rep_matches_translation_file, 'w') as trans_fasta:            
            for subject_id in subjects_ids:
                trans_fasta.writelines(">"+subject_id+"\n")
                trans_fasta.writelines(str(reps_translation_dict[subject_id])+"\n")


    # Calculate BSR based on BLASTp.
    total_blasts = len(blastp_runs_to_do)
    bsr_values = {}
    # Create query entries
    for query in blastp_runs_to_do:
        bsr_values[query] = {}
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches)        
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                blastp_runs_to_do, 
                                repeat('blastp'),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):
            
            filtered_alignments_dict, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True)
            # Since BLAST may find several local aligments choose the largest one to calculate BSR
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    largest_alignment = 0
                    for entry_id, result in results.items():
                        if result['query_length'] > largest_alignment:
                            largest_alignment = result['query_length']
                            bsr_values[query].update({subject_id: bf.compute_bsr(result['score'], self_score)})

            print(f"Running BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts}")
            i += 1

    # Add kmer cov, kmer sim and frequency of the cds in the genomes
    for query, subjects_dict in list(representative_blast_results.items()):
        for subject, blastn_results in list(subjects_dict.items()):
            # Some results may not have kmers matches or BSR values so put them
            # as 0
            if subject in reps_kmers_sim[query]:
                sim, cov = reps_kmers_sim[query][subject]
            else:
                sim = 0
                cov = 0
            if subject in bsr_values[query]:
                bsr = bsr_values[query][subject]
                # For some reason some isolates have bsr slighty higher than 1
                # related to blast database and sequences used.
                if bsr > 1.0:
                    bsr = float(round(bsr))
            # If subject not in inside queries alignments. This means that even
            # though there was an BLASTn align, they didn´t align when BLASTp
            # was employed.
            else:
                bsr = 0
                
            # Calculate total alignment for all of the fragments of BLASTn
            # if there more than one BLASTn alignments
            # For query and subject
            if len(blastn_results) > 1:
                total_length = {}
                gaps_all = representative_blast_results_coords_all[query][subject].pop('gaps')
                gaps_pident = representative_blast_results_coords_pident[query][subject].pop('gaps')
                for ref, intervals in representative_blast_results_coords_all[query][subject].items():
                    # Sort by start position
                    sorted_intervals = [interval for interval in sorted(intervals,
                                                                        key=lambda x: x[0])]
                    # Merge alignments and calculate the total length of BLASTn alignments.
                    length = sum([(interval[1] - interval[0] + 1) for interval in af.merge_intervals(sorted_intervals)])
                    total_length[ref] = length
                # Calculate global palign
                global_palign_all = min([(total_length['query'] + gaps_all) / blastn_results[1]['query_length'],
                                     (total_length['subject']) / blastn_results[1]['subject_length']])
                
                for ref, intervals in representative_blast_results_coords_pident[query][subject].items():
                    if len(intervals) != 0:
                        # Sort by start position
                        sorted_intervals = [interval for interval in sorted(intervals,
                                                                            key=lambda x: x[0])]
                        # Merge alignments and calculate the total length of BLASTn alignments.
                        length = sum([(interval[1] - interval[0] + 1) for interval in af.merge_intervals(sorted_intervals)])
                        total_length[ref] = length
                    else:
                        # Since we are considering soem values that may not pass pident threshold
                        # it means some may not have any interval based on pident
                        total_length[ref] = 0
                # Calculate global palign
                global_palign_pident_min = min([(total_length['query'] - gaps_pident) / blastn_results[1]['query_length'],
                                                (total_length['subject']) / blastn_results[1]['subject_length']])
                
                global_palign_pident_max = max([(total_length['query'] - gaps_pident) / blastn_results[1]['query_length'],
                                                (total_length['subject']) / blastn_results[1]['subject_length']])
            # If there is only one results
            else:
                global_palign_all = min([(blastn_results[1]['query_end'] - blastn_results[1]['query_start'] + 1 - result['gaps']) / blastn_results[1]['query_length'],
                                    (blastn_results[1]['subject_end'] - blastn_results[1]['subject_start'] + 1) / blastn_results[1]['subject_length']])
                
                global_palign_pident_min = min([(blastn_results[1]['query_end'] - blastn_results[1]['query_start'] + 1 - result['gaps']) / blastn_results[1]['query_length'],
                                    (blastn_results[1]['subject_end'] - blastn_results[1]['subject_start'] + 1) / blastn_results[1]['subject_length']])
                
                global_palign_pident_max = max([(blastn_results[1]['query_end'] - blastn_results[1]['query_start'] + 1 - result['gaps']) / blastn_results[1]['query_length'],
                                    (blastn_results[1]['subject_end'] - blastn_results[1]['subject_start'] + 1) / blastn_results[1]['subject_length']])

            # Create update dict with the values to add
            update_dict = {'bsr' : bsr,
                           'kmers_sim': sim,
                           'kmers_cov': cov,
                           'frequency_in_genomes_query_cds' : frequency_cds_cluster[query],
                           'frequency_in_genomes_subject_cds' : frequency_cds_cluster[subject],
                           'global_palign_all': global_palign_all,
                           'global_palign_pident_min': global_palign_pident_min,
                           'global_palign_pident_max': global_palign_pident_max}
            
            for entry_id, result in list(blastn_results.items()):
                # Calculate local Palign particular to that BLASTn match
                local_palign = min([(result['query_end'] - result['query_start'] + 1 + result['gaps']) / result['query_length'],
                                    (result['subject_end'] - result['subject_start'] + 1) / result['subject_length']])
                # Verify if inverse alignment was made
                if local_palign >= 0:
                    # update the update dict
                    update_dict.update({'local_palign' : local_palign})
                    # Add everything to the dict
                    representative_blast_results[query][subject][entry_id].update(update_dict)
                # Remove palign that is negative, meaning tha reverse blast alignment was made
                else:
                    del representative_blast_results[query][subject][entry_id]
            # Remove empty subjects dicts
            if len(representative_blast_results[query][subject]) == 0:
                del representative_blast_results[query][subject]
        # Remove empty query dicts
        if len(representative_blast_results[query]) == 0:
            del representative_blast_results[query]

    print("Filtering BLAST results into classes...")
    results_output = os.path.join(output_directory, "3_Results_files")
    ff.create_directory(results_output)
    blastn_processed_results_path = os.path.join(results_output, "Processed_Blast")
    ff.create_directory(blastn_processed_results_path)
    report_file_path = os.path.join(blastn_processed_results_path, "blast_all_matches.tsv")
    
    # Separate results into different classes
    results_outcome, classes_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                                                       constants)
    # Write all of the BLASTn results to a file
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("Processing classes...")
    # Process the results_outcome dict and write individual classes to TSV file
    [results_outcome,
     relationships,
     cluster_dict_1a] = process_classes(results_outcome)
    print("Writting classes results to files...")
    write_processed_results_to_file(results_outcome, relationships,
                                    representative_blast_results, cluster_dict_1a,
                                    classes_outcome, blastn_processed_results_path)
    
    print("Wrapping up results...")
    wrap_up_results(schema, results_outcome, not_included_cds, clusters, 
                    blastn_processed_results_path, blastn_processed_results_path, cpu)