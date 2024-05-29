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
                                      core_functions as cof,  
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf,
                                      pandas_functions as pf)

def main(schema, allelecall_directory, output_directory, constants, cpu):
    frequency_in_genomes = {}
    loci_ids = [True, True]
    cof.process_schema(schema,
                       [],
                       output_directory,
                       None,
                       None,
                       frequency_in_genomes,
                       allelecall_directory,
                       None,
                       loci_ids,
                       True,
                       constants,
                       cpu)