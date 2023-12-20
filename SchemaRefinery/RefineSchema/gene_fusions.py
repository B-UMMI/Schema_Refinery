import os
import sys
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff, 
                       sequence_functions as sf,
                       clustering_functions as cf, 
                       blast_functions as bf, 
                       alignments_functions as af)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff, 
                                      sequence_functions as sf, 
                                      clustering_functions as cf, 
                                      blast_functions as bf, 
                                      alignments_functions as af)

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
                                                                                                       
    clusters, reps_sequences, reps_groups  = cf.minimizer_clustering(cds_translation_dict, 
                                                                  5, 5, True, 1, clusters, 
                                                                  reps_sequences, reps_groups, 
                                                                  20, clustering_sim, clustering_cov, 
                                                                  True)
    print("Filtering clusters...")
    singleton_clusters = {}
    filtered_clusters = {}
    
    #Separate singletons and clusters with more than one protein
    for k, v in clusters.items():
        if len(v) > 1:
            filtered_clusters[k] = v
        else:
            singleton_clusters[k] = v
    
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
        for cluster_rep_id in filtered_clusters:
            
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
    representative_blast_results = []
    i=1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing, 
                                filtered_clusters, repeat('blastn'), 
                                repeat(blastn_results_folder), 
                                repeat(rep_paths), 
                                repeat(representatives_all_fasta_file)):
            
            alignment_strings, filtered_alignments_dict = af.process_blast_results(res[1], constants_threshold)
            representative_blast_results.append([res[0], [alignment_strings, filtered_alignments_dict]])

            print(f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps}")
            i+=1



