import os
try:
    from RefineSchema import (paralagous_loci, 
                              gene_fusions)

except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema import (paralagous_loci, 
                                             gene_fusions)


def main(schema, output_directory, allelecall_directory, alignment_ratio_threshold_paralagous, 
         pident_threshold_paralagous, alignment_ratio_threshold_gene_fusions, 
         pident_threshold_gene_fusions, clustering_sim, clustering_cov, cpu):

    if True:
        print("Identifying paralagous loci...")
        paralagous_loci_output = os.path.join(output_directory, "paralagous_loci")
        paralagous_loci.main(schema, paralagous_loci_output, alignment_ratio_threshold_paralagous, 
                             pident_threshold_paralagous, cpu)

    if allelecall_directory:
        print("Identifying genes fusions...")
        gene_fusions_output = os.path.join(output_directory, "gene_fusions")
        gene_fusions.main(schema, gene_fusions_output, allelecall_directory, 
                          clustering_sim, clustering_cov, alignment_ratio_threshold_gene_fusions, 
                          pident_threshold_gene_fusions, cpu)