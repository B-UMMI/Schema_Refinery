import os
import sys
import concurrent.futures
from itertools import repeat

try:
    from utils.file_functions import check_and_delete_file, create_directory
    from utils.sequence_functions import seq_to_hash, read_fasta_file_iterator
    from utils.kmers_functions import determine_minimizers
    from utils.sequence_functions import translate_dna
    from utils.clustering_functions import minimizer_clustering
except ModuleNotFoundError:
    from SchemaRefinery.utils.file_functions import check_and_delete_file, create_directory
    from SchemaRefinery.utils.sequence_functions import seq_to_hash, read_fasta_file_iterator
    from SchemaRefinery.utils.kmers_functions import determine_minimizers
    from SchemaRefinery.utils.sequence_functions import translate_dna
    from SchemaRefinery.utils.clustering_functions import minimizer_clustering

def hash_sequences(file_path):
    """
    Hashes sequences in fasta file based on input file_path.

    Parameters
    ----------
    file_paths : str
        Contains file path to the fasta files.

    Returns
    -------
    hash_list : set
        Returns a list containing all of the sequences hashes present in the input files.
    """

    hash_set = set()
    for rec in read_fasta_file_iterator(file_path):
        hash_set.add(seq_to_hash(str(rec.seq)))

    return hash_set

def fetch_not_included_cds(all_schema_hashes, file_path_cds):
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

    for rec in read_fasta_file_iterator(file_path_cds):
        print(f"Processed {i} CDS")
        i += 1
        seq_hash = seq_to_hash(str(rec.seq))
        if seq_hash not in all_schema_hashes:
            not_included_cds[rec.id] = rec.seq

    return not_included_cds

def main(schema, output_directory, allelecall_directory,clustering_sim, 
         clustering_cov, cpu):
    create_directory(output_directory)
    file_path_cds = os.path.join(allelecall_directory, "temp", "3_cds_preprocess", "cds_deduplication", "distinct_cds_merged.fasta")
    if not os.path.exists(file_path_cds):
        sys.exit(f"Error: {file_path_cds} must exist, make sure that AlleleCall was run using --no-cleanup flag.")

    print("Identifying CDS not present in the schema...")
    schema_file_paths = {f.replace(".fasta", ""): os.path.join(schema, f) for f in os.listdir(schema) if not os.path.isdir(f) and f.endswith(".fasta")}
    all_schema_hashes = set()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(hash_sequences, schema_file_paths.values()):
            all_schema_hashes.update(res)

    not_included_cds = fetch_not_included_cds(all_schema_hashes, file_path_cds)
    
    file_invalid_cds = os.path.join(allelecall_directory,"invalid_cds.txt")
    with open(file_invalid_cds, 'r') as invalid:
        for invalid_line in invalid.read().splitlines():
            del not_included_cds[invalid_line.split(":")[0]]
        print(f"Removing invalid CDS.")
        
    print(f"Identified {len(not_included_cds)} valid CDS not present in the schema")

    cds_not_present_file_path = os.path.join(output_directory, "CDS_not_found.fasta")
    with open(cds_not_present_file_path, 'w') as cds_not_found:
        for id, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id+"\n")
            cds_not_found.writelines(str(sequence)+"\n")


    print("Translate not found CDS...")
    cds_translation_dict = {}
    cds_not_present_translation_file_path = os.path.join(output_directory, "CDS_not_found_translation.fasta")
    protein_hashes = []
    i = 1
    total = len(not_included_cds)
    with open(cds_not_present_translation_file_path, 'w') as translation:
        for id_s, sequence in not_included_cds.items():
            print(f"Translated {i}/{total} CDS")
            i += 1
            protein_translation = str(translate_dna(str(sequence), 11, 0, True)[0][0])
            prot_hash = seq_to_hash(protein_translation)
            if prot_hash not in protein_hashes:
                protein_hashes.append(prot_hash)
                cds_translation_dict[id_s] = protein_translation
                translation.writelines(id_s+"\n")
                translation.writelines(protein_translation+"\n")

    schema_short = os.path.join(schema, "short")
    schema_short_files_paths = {f.replace("_short.fasta", ""): os.path.join(schema_short, f) for f in os.listdir(schema_short) if not os.path.isdir(f) and f.endswith(".fasta")}

    print("Extracting minimizers for the translated sequences and clustering...")

    reps_groups = {}
    clusters = {}
    reps_sequences = {}
    clusters, reps_sequences, reps_groups  = minimizer_clustering(cds_translation_dict, 
                                                                  5, 5, True, 1, clusters, 
                                                                  reps_sequences, reps_groups, 
                                                                  20, clustering_sim, clustering_cov, 
                                                                  True)
    print("Filtering clusters...")
    for cluster in clusters:
        if cluster
        


