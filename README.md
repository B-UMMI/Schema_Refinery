# Schema_Refinery

Set of scripts and instructions to refine wg/cgMLST schemas.

## alleles_frequency.py

This script calculates alleles frequencies.

## annotations_merger.py

This script merges the annotations from different sources
(reference proteomes, Uniprot and Genbank files) into one single file.

## download_ncbi_assemblies.py

### Description

Accepts a Genome Assembly and Annotation report table from the NCBI and downloads the genome assemblies of the samples listed
in the table. The default behavior will download files from RefSeq. If the genome assembly or file that has to be downloaded
cannot be found in RefSeq, it will try to download from Genbank. The full list of species can be viewed at
https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/ and the script accepts a table with a list of assemblies for a single
species (e.g.: https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/186/). Files are download in GZIP format. It is
possible to specify a file extension/suffix to download genome assembly data in a specific file format. You can download the
following file types from the FTP server:

- "assembly_report.txt"
- "assembly_stats.txt"
- "cds_from_genomic.fna.gz"
- "feature_count.txt.gz"
- "feature_table.txt.gz"
- "genomic.fna.gz"
- "genomic.gbff.gz"
- "genomic.gff.gz"
- "genomic.gtf.gz"
- "protein.faa.gz"
- "protein.gpff.gz"
- "rna_from_genomic.fna.gz"
- "translated_cds.faa.gz"

#### Usage

Download genome assemblies in fasta.gz format:

```py
python download_ncbi_assemblies.py -t <input_table> -o <output_directory>
```

Download compressed GFF files:

```py
python download_ncbi_assemblies.py -t <input_table> -o <output_directory> --fe genomic.fna.gz
```

Download only the genome assemblies that are available in RefSeq:

```py
python download_ncbi_assemblies.py -t <input_table> -o <output_directory> --ftp refseq
```

Download only the genome assemblies that are available in RefSeq using 4 threads

```py
python download_ncbi_assemblies.py -t <input_table> -o <output_directory> --ftp refseq --threads 4
```

## download_reference_proteomes.py

Accepts a Proteomes report table from Uniprot and downloads
the proteomes of the samples listed in the table.

## fasta_header_integer.py

Converts fasta headers to consecutive integers, useful for the creation of BLAST databases.

## genbank_annotations_2.py

This script extracts annotations for a schema from genbank files.

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

## protein_diversity.py

Given an input schema, and for each gene file, translates DNA
sequences of all alleles in each locus and determines
which alleles code for the exact same protein sequence.

## proteome_matcher.py

This script matches TrEMBL and SwissProt annotations to schemas.

## proteome_splitter.py

This script splits a set of Uniprot's reference proteomes into a pair of files,
one file with sequences from Swiss-Prot and another with sequences from TrEMBL.

## remove_loci.py

This script removes loci from a schema.

## schema_validation_functions.py

This script contains a collection of fucntions that are used to validade schemas.

## annotation_merger.py

This script merges annotations from different sources into one single comprehensible file.

## assembly_statistics.py

This script calculates statistics for input assemblies.

## paralogous_interloci_validation_merger.py

This script merges paralagous loci into groups obtained from inter_loci_validation.py.

## allele_call_paralogous_merger.py

This script merges paralagous loci into group from allele call module output.

## split_locus.py

This script splits loci based on sequence size intervals.

## reassign_ids.py

This script reassigns allele identifiers to ensure that the set of alleles for a locus have sequential identifiers.

## new_loci_annotations.py

This script uses ids identifiers to match with older annotations creating annotation table for newer schemas.

## missing_ids.py

This script determines if any locus in a schema does not have sequential allele identifiers.

## mask_matrix.py

Masks elements in a matrix with allelic profiles created with chewBBACA.

## locus_length_variability.py

This scripts computes allele length statistics for all loci in a schema, identifying loci with multimodal allele length distributions and loci with alleles that deviate from the mode value.

## gc_content.py

This script outputs a more detailed GC content information than assembly_statistics.py.

## Extract_cgAlleles.py

This module determines the set of genes in the core genome based on a matrix with allelic profiles and a threshold that defines the proportion of genomes a gene must be present in to be included in the core genome.
