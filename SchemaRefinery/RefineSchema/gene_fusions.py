import os
import sys
import concurrent.futures
import copy
from itertools import repeat

try:
    from utils import (file_functions as ff, 
                       sequence_functions as sf,
                       clustering_functions as cf, 
                       blast_functions as bf, 
                       alignments_functions as af,
                       kmers_functions as kf,
                       list_functions as lf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff, 
                                      sequence_functions as sf, 
                                      clustering_functions as cf, 
                                      blast_functions as bf, 
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      list_functions as lf)

def fetch_not_included_cds(file_path_cds):
    """
    Compares the hashes list with the hashes obtained from cds.

    Parameters
    ----------
    hash_list : list
        List containing all of the sequences hashes present in the input files.

    Returns
    -------
    not_included_cds : dict
        Returns dict with key as fasta header and value as fasta sequence.
    """
    not_included_cds = {}
    i = 1

    for rec in sf.read_fasta_file_iterator(file_path_cds):
        print(f"Processed {i} CDS")
        i += 1
        not_included_cds[rec.id] = rec.seq

    return not_included_cds

def alignment_string_dict_to_file(alignment_string_dict, file_path):
    """
    Writes alignments strings to file.

    Parameters
    ----------
    alignment_string_dict : dict
        dict containing alignments string to write to file
    file_path : str
        File path to create to write the file
    """
    
    with open(file_path, 'w') as report_file:
        report_file.writelines(["Query\t", "Subject\t", "Query Start-End\t", "Subject Start-End\t", "Query Biggest Alignment Ratio\t", 
                            "Subject Biggest Alignment Ratio\t", "Query Length\t", "Subject Length\t", "Number of Gaps\t", 
                            "Pident - Percentage of identical matches\t", "Kmer sim\t", "Kmer cov\n"])

        for res in alignment_string_dict.values():
            report_file.writelines(res.values())
            
def separate_blastn_results_into_classes(representative_blast_results, representative_alignment_strings, path):
    """
    Separates one BLASTn dict into various dict depending on the class.

    Parameters
    ----------
    representative_blast_results : dict
        dict containing blast results to filter
    representative_alignment_strings : dict
        dict containing alignments string to filter
    path : str
        Dir path to write the files for each class of results
        
    Returns
    -------
    c1 : dict
        dict containing filtered blast results to class 1
    c2 : dict
        dict containing filtered blast results to class 2
    c1_strings : dict
        dict containing alignments string for class 1
    c2_strings : dict
        dict containing alignments string for class 1
    """
    
    representative_blast_results_filtered = copy.deepcopy(representative_blast_results)
    
    #Create various dicts and split the results into different classes
    cluster_classes = {}
    cluster_classes_strings = {}
    
    for class_ in ["similar_size","different_size","None_assigned"]:
        cluster_classes[class_] = {}
        cluster_classes_strings[class_] = {}
        #create dict for each rep inside cluster_classes
        for query in representative_blast_results_filtered.keys():
            cluster_classes[class_][query] = {}
            cluster_classes_strings[class_][query] = {}
    
    for query, rep_b_result in representative_blast_results_filtered.items():
        for id_entry, reps in rep_b_result.items():
            reps = reps[0]
            length_threshold = 0.05
            high = reps['subject_length'] * (1+length_threshold)
            low = reps['subject_length'] * (1-length_threshold)
            
            pident = reps['pident']

            if type(pident) == str:
                pident_list = pident.split(';')
                pident = sum(pident_list)/len(pident_list)
                    
            if (low > reps['query_length'] or reps['query_length'] > high) and pident >= 90:
                cluster_classes["different_size"][query].update({id_entry : reps})
                cluster_classes_strings["different_size"][query].update({id_entry : representative_alignment_strings[query][id_entry]})
                del representative_blast_results[query][id_entry]
                del representative_alignment_strings[query][id_entry]

            elif low <= reps['query_length'] <= high and pident >= 90:
                cluster_classes["similar_size"][query].update({id_entry : reps})
                cluster_classes_strings["similar_size"][query].update({id_entry : representative_alignment_strings[query][id_entry]})
                
                del representative_blast_results[query][id_entry]
                del representative_alignment_strings[query][id_entry]
            
            else:
                cluster_classes["None_assigned"][query].update({id_entry : reps})
                cluster_classes_strings["None_assigned"][query].update({id_entry : representative_alignment_strings[query][id_entry]})
                
                del representative_blast_results[query][id_entry]
                del representative_alignment_strings[query][id_entry]        
        
        for class_ in cluster_classes.keys():
            if len(cluster_classes[class_][query]) == 0:
                del cluster_classes[class_][query]
                del cluster_classes_strings[class_][query]
                
        #If after separating BLASTn results into different classes the cluster rep
        #is empty, then delete ir from dict containing all of the results
        if len(representative_blast_results[query]) == 0:
            del representative_blast_results[query]
            del representative_alignment_strings[query]

                
    for k, v in cluster_classes_strings.items():
        report_file_path = os.path.join(path, f"blastn_group_{k}.tsv")
        #write individual group to file
        alignment_string_dict_to_file(v, report_file_path)

    after_filter = os.path.join(path, "blastn_filtered.tsv")
    alignment_string_dict_to_file(representative_alignment_strings, after_filter)
    
    return cluster_classes
    
def main(schema, output_directory, allelecall_directory, clustering_sim, 
         clustering_cov, alignment_ratio_threshold_gene_fusions, 
         pident_threshold_gene_fusions, cpu):
    
    constants_threshold = [alignment_ratio_threshold_gene_fusions, pident_threshold_gene_fusions]
    
    ff.create_directory(output_directory)
    file_path_cds = os.path.join(allelecall_directory, "unclassified_sequences.fasta")
    if not os.path.exists(file_path_cds):
        sys.exit(f"Error: {file_path_cds} must exist, make sure that AlleleCall "
                 "was run using --no-cleanup and --output-unclassified flag.")

    print("Identifying CDS not present in the schema...")

    not_included_cds = fetch_not_included_cds(file_path_cds)
        
    print(f"Identified {len(not_included_cds)} valid CDS not present in the schema")

    cds_output = os.path.join(output_directory, "CDS_processing")
    ff.create_directory(cds_output)
    cds_not_present_file_path = os.path.join(cds_output, "CDS_not_found.fasta")
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id+"\n")
            cds_not_found.writelines(str(sequence)+"\n")


    print("Translate unclassified CDS...")
    cds_translation_dict = {}
    cds_not_present_translation_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    protein_hashes = {}
    i = 1
    total = len(not_included_cds)
    with open(cds_not_present_translation_file_path, 'w+') as translation:
        for id_s, sequence in not_included_cds.items():
            print(f"Translated {i}/{total} CDS")
            i += 1
            protein_translation = str(sf.translate_dna(str(sequence), 11, 0, True)[0][0])
            prot_hash = sf.seq_to_hash(protein_translation)
            if prot_hash not in protein_hashes:
                protein_hashes[prot_hash] = [id_s]
                cds_translation_dict[id_s] = protein_translation
                translation.writelines(id_s+"\n")
                translation.writelines(protein_translation+"\n")
            else:
                protein_hashes[prot_hash].append(id_s)

    schema_short = os.path.join(schema, "short")
    schema_short_files_paths = {f.replace("_short.fasta", ""): os.path.join(schema_short, f) 
                                for f in os.listdir(schema_short) 
                                if not os.path.isdir(f) and f.endswith(".fasta")}

    print("Extracting minimizers for the translated sequences and clustering...")

    reps_groups = {}
    clusters = {}
    reps_sequences = {}
    
    #sort by size of proteins
    cds_translation_dict = {k: v for k, v in sorted(cds_translation_dict.items(),
                                                     key=lambda x: len(x[1]),
                                                     reverse=True)}
                                                                                                       
    clusters, reps_sequences, reps_groups, prot_len_dict  = cf.minimizer_clustering(cds_translation_dict, 
                                                                                    5, 5, True, 1, clusters, 
                                                                                    reps_sequences, reps_groups, 
                                                                                    20, clustering_sim, clustering_cov, 
                                                                                    True)
        
    print("Filtering clusters...")
    singleton_clusters = {}
    filtered_clusters = {}
    
    #Separate singletons and clusters with more than one protein and get size
    cluster_size = {}
    for k, v in clusters.items():
        sizes = lf.get_max_min_values([entry_data[2] for entry_data in v])
        cluster_size[k] = sizes[1]/sizes[0]
        if len(v) > 1:
            filtered_clusters[k] = v
        else:
            singleton_clusters[k] = v
    
    print("Retrivieng kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items() 
                             if rep_id in clusters}
    
    reps_translation_dict = {k: v for k, v in sorted(reps_translation_dict.items(),
                                                     key=lambda x: len(x[1]),
                                                     reverse=True)}
    #recalculate the sim and cov between reps
    for cluster_id in reps_translation_dict:
        kmers_rep = set(kf.determine_minimizers(reps_translation_dict[cluster_id], 
                                                5, 5, 1 ,True, True))
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep, 
                                                               reps_groups, 
                                                               0, 
                                                               0,
                                                               prot_len_dict, 
                                                               cluster_id, 5)
        
        reps_kmers_sim[cluster_id] = {match_values[0]:match_values[1:] 
                                      for match_values in reps_kmers_sim[cluster_id]}
            
    print("Building BLASTn database...")
    blastn_output = os.path.join(output_directory, "BLASTn_processing")
    ff.create_directory(blastn_output)
    #write fasta file with all of the representatives sequences   
    representatives_blast_folder = os.path.join(blastn_output,
                                                "cluster_representatives_fastas")
    ff.create_directory(representatives_blast_folder)

    representatives_all_fasta_file = os.path.join(representatives_blast_folder,
                                                  "all_cluster_representatives.fasta")
    
    rep_paths = {}
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for cluster_rep_id in clusters:
            
            all_fasta.writelines(">"+cluster_rep_id+"\n")
            all_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")
            
            rep_fasta_file = os.path.join(representatives_blast_folder, 
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths[cluster_rep_id] = rep_fasta_file
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.writelines(">"+cluster_rep_id+"\n")
                rep_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")
            
    blastn_results_folder = os.path.join(blastn_output, "blast_results")
    ff.create_directory(blastn_results_folder)
    
    total_reps = len(rep_paths)
    representative_blast_results = {}
    representative_alignment_strings = {}
    i=1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing, 
                                clusters.keys(), repeat('blastn'), 
                                repeat(blastn_results_folder), 
                                repeat(rep_paths), 
                                repeat(representatives_all_fasta_file)):
            
            alignment_strings, filtered_alignments_dict = af.process_blast_results(res[1], constants_threshold)
            
            if len(alignment_strings) > 1:
                representative_blast_results[res[0]] = filtered_alignments_dict
                representative_alignment_strings[res[0]] = alignment_strings

            print(f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps}")
            i+=1
            
    #Reformat output of af.process_blast_results[1] and add kmer cov and sim
    for key, alignment_string in representative_alignment_strings.items():
        dict_alignment_string = {}
        for string in alignment_string:
            split_string = string.split('\t')
            dict_entry = split_string[0]+';'+split_string[1]
            dict_alignment_string[dict_entry] = string
            #get sim and cov
            if split_string[1] in reps_kmers_sim[key]:
                sim, cov = reps_kmers_sim[key][split_string[1]]
            else:
                sim = 0
                cov = 0
                
            update_dict = {'kmers_sim' : sim,
                           'kmers_cov' : cov}
            #Add kmer cov and sim to strings for so it also writes into a file
            dict_alignment_string[dict_entry] = string.replace('\n','\t') + '\t'.join([str(sim),str(cov)]) + '\n'
            representative_blast_results[key][dict_entry][0].update(update_dict)
                
        representative_alignment_strings[key] = dict_alignment_string
        
    #write all relevant rep BLASTn matches to file
    blastn_processed_results = os.path.join(blastn_output, "Processed_Blastn")
    ff.create_directory(blastn_processed_results)
    report_file_path = os.path.join(blastn_processed_results, "blastn_all_matches.tsv")
    alignment_string_dict_to_file(representative_alignment_strings, report_file_path)
            
    
    print("Filtering BLASTn results into subclusters...")
    cluster_classes = separate_blastn_results_into_classes(representative_blast_results, 
                                                            representative_alignment_strings, 
                                                            blastn_processed_results)
            
            
