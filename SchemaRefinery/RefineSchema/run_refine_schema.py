import os
import sys
from typing import Any

try:
    from RefineSchema import (IdentifySpuriousGenes)

except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema import (IdentifySpuriousGenes)


import os
import sys
from typing import List, Any

def main(schema: str, 
         output_directory: str, 
         allelecall_directory: str, 
         alignment_ratio_threshold_paralagous: float, 
         pident_threshold_paralagous: float, 
         alignment_ratio_threshold_gene_fusions: float, 
         pident_threshold_gene_fusions: float, 
         clustering_sim: float, 
         clustering_cov: float, 
         genome_presence: float,
         size_threshold: int, 
         cpu: int) -> None:

    if allelecall_directory:
        temp_paths: List[str] = [os.path.join(allelecall_directory, "temp"), 
                                 os.path.join(allelecall_directory, "unclassified_sequences.fasta"),
                                 os.path.join(allelecall_directory, "missing_classes.fasta")]
        # Put all constants in one dict in order to decrease number of variables
        # used around.
        constants: List[Any] = [alignment_ratio_threshold_gene_fusions, 
                                pident_threshold_gene_fusions,
                                genome_presence,
                                clustering_sim,
                                clustering_cov,
                                size_threshold]
        
        if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
            sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                     "was run using --no-cleanup and --output-unclassified flag.")
        
        print("Identifying spurious loci in the schema...")
        spurious_loci_output: str = os.path.join(output_directory, "spurious_loci")
        IdentifySpuriousGenes.identify_spurious_genes(schema_path, output_directory, allelecall_directory, possible_new_loci, constants, temp_paths, run_mode, processing_mode, cpu)