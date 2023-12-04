try:
    from RefineSchema import paralagous_loci
except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema import paralagous_loci


def main(schema, output_directory, alignment_ratio_threshold, pident_threshold, cpu):

    paralagous_loci.main(schema, output_directory, alignment_ratio_threshold, pident_threshold, cpu)