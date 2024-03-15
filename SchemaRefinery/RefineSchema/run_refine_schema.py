import os
import sys
try:
    from RefineSchema import (paralagous_loci, 
                              unclassified_cds,
                              spurious_loci_identification)

except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema import (paralagous_loci, 
                                             unclassified_cds,
                                             spurious_loci_identification)


def main(schema, output_directory, allelecall_directory, alignment_ratio_threshold_paralagous, 
         pident_threshold_paralagous, alignment_ratio_threshold_gene_fusions, 
         pident_threshold_gene_fusions, clustering_sim, clustering_cov, genome_presence,
         size_threshold, cpu):

    if True:
        print("Identifying paralagous loci...")
        paralagous_loci_output = os.path.join(output_directory, "paralagous_loci")
        paralagous_loci.main(schema, paralagous_loci_output, alignment_ratio_threshold_paralagous, 
                             pident_threshold_paralagous, cpu)

    if allelecall_directory:
        temp_paths = [os.path.join(allelecall_directory, "temp"), 
                      os.path.join(allelecall_directory, "unclassified_sequences.fasta")]
        # Put all constants in one dict in order to decrease number of variables
        # used around.
        constants = [alignment_ratio_threshold_gene_fusions, 
                     pident_threshold_gene_fusions,
                     genome_presence,
                     clustering_sim,
                     clustering_cov,
                     size_threshold]
        
        if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
            sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                     "was run using --no-cleanup and --output-unclassified flag.")
            
        print("Identifying genes fusions...")
        unclassified_cds_output = os.path.join(output_directory, "unclassified_cds")
        unclassified_cds.main(schema, unclassified_cds_output, allelecall_directory,
                              constants,temp_paths, cpu)
        
        print("Identifying spurious loci in the schema...")
        unclassified_cds_output = os.path.join(output_directory, "spurious_loci")
        spurious_loci_identification.main(schema, output_directory,
                                          allelecall_directory, cpu)