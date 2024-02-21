import os
import sys
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
                       list_functions as lf,
                       graphical_functions as gf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      list_functions as lf,
                                      graphical_functions as gf)

def main(schema_path_folder, output_directory, allelecall_directory, temp_paths, cpu):
        
    print("Reading loci FASTA files...")
    # Get all of the schema loci FASTA files path
    schema_loci = {loci_path.replace(".fasta", ""): loci_path 
                   for loci_path in os.listdir(schema_path) 
                   if loci_path.endswith('.fasta')}
    
    # Get all of the schema loci short FASTA files path
    schema_short_path = os.path.join(schema_path, 'short')
    schema_loci_short = {loci_path.replace(".fasta", ""): loci_path 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    print(f"Found {len(schema_loci_short)} loci.")
    
    
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                blastp_runs_to_do.keys(), 
                                repeat('blastp'),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):