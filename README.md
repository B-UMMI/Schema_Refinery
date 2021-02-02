# Schema_Refinery

Set of scripts and instructions to refine wg/cgMLST schemas.

## annotations_merger.py

This script merges the annotations from different sources
(reference proteomes, Uniprot and Genbank files) into one single file.

## download_ncbi_assemblies.py

This script accepts a Genome Assembly and Annotation report
table from the NCBI and downloads the genomes/assemblies of
the samples listed in the table.

## genbank_annotations.py

This script extracts annotations for a schema from protein
genbank files.

## inter_loci_validation.py

This script verifies if a set of FASTA files corresponding to different loci
respect the specified BSR and length thresholds applied by chewBBACA's algorithm
during schema creation or allele calling. It outputs files with relevant information
that help identify problematic cases in which representative alleles from different
loci/genes match with a BSR >= 0.6.

## intra_loci_validation.py

This script checks if the set of alleles for each gene/locus file given as
input respects the specified BSR and length thresholds applied by chewBBACA's
algorithm during schema creation or allele calling.

## loci_positions.py

This script determines the genome position of loci first representatives.

## locus_pairwise_graph.py

This script creates a locus graph based on BSR.

## match_schemas.py

This script compares two schemas and outputs the loci matches.

## merge_loci.py

This script merges loci from a schema that are very similar.

## proteome_matcher.py

This script matches TrEMBL and SwissProt annotations to schemas.

## proteome_splitter.py

This script splits a set of Uniprot's reference proteomes into a pair of files,
one file with sequences from Swiss-Prot and another with sequences from TrEMBL.

## remove_loci.py

This script removes loci from a schema.

## schema_validation_functions.py

This script contains a collection of fucntions that are used to validade schemas.
