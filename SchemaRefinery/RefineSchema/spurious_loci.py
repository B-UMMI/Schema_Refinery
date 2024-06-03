try:
    from utils import (core_functions as cof,)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (core_functions as cof)

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