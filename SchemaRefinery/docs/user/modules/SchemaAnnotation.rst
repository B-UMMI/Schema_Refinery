SchemaAnnotation - Annotate schemas in a directory
==================================================

Description
-----------

The `SchemaAnnotation` module parses command-line arguments and initiates the schema annotation process. This module sets up an argument parser to handle various command-line options for annotating schemas and then calls the main function of the `SchemaAnnotation` class with the parsed arguments.

Features
--------

- Annotating schemas in a directory.
- Configurable parameters for the annotation process.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `SchemaAnnotation` module can be used as follows:

.. code-block:: bash

    SR SchemaAnnotation -s /path/to/schema -o /path/to/output -ao uniprot-proteomes genbank-files -pt path/to/proteome/table -gf path/to/genbank/files

Command-Line Arguments
----------------------

-s, --schema-directory
    (Required) Path to the schema's directory.

-o, --output-directory
    (Required) Path to the output directory where to save the files.

-ao, --annotation-options
    (Required) Annotation options to run.
    Choices: uniprot-proteomes, genbank-files, uniprot-sparql, match-schemas.

-pt, --proteome-table
    (Optional) TSV file downloaded from UniProt that contains the list of proteomes.

-gf, --genbank-files
    (Optional) Path to the directory that contains Genbank files with annotations to extract.

-ca, --chewie-annotations
    (Optional) File with the results from chewBBACA UniprotFinder module.

-ss, --subject-schema
    (Optional) Path to the subject schema directory. This argument is needed by the Match Schemas sub-module.

--bsr
    (Optional) Minimum BSR value to consider aligned alleles as alleles for the same locus.
    Default: 0.6

-t, --threads
    (Optional) Number of threads for concurrent download.
    Default: 1

-c, --cpu
    (Optional) Number of CPU cores for multiprocessing.
    Default: 1

-r, --retry
    (Optional) Maximum number of retries when a download fails.
    Default: 7

-tt, --translation-table
    (Optional) Translation table to use for the CDS translation.
    Default: 11

-cs, --clustering-sim
    (Optional) Similarity value for kmers representatives (float: 0-1).
    Default: 0.9

-cc, --clustering-cov
    (Optional) Coverage value for kmers representatives (float: 0-1).
    Default: 0.9

-sr, --size_ratio
    (Optional) Size ratio to consider alleles as the same locus.
    Default: 0.8

-rm, --run-mode
    (Optional) Mode to run the module.
    Choices: reps, alleles.
    Default: reps

-pm, --processing-mode
    (Optional) Mode to run the module for Schema match.
    Choices: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.
    Default: None

-egtc, --extra_genbank_table_columns
    (Optional) List of columns to add to annotation file.
    Default: []

-gia, --genbank-ids-to-add
    (Optional) List of GenBank IDs to add to final results.
    Default: []

-pia, --proteome-ids-to-add
    (Optional) List of Proteome IDs to add to final results.
    Default: []

--nocleanup
    (Optional) Flag to indicate whether to skip cleanup after running the module.

Examples
--------

Here are some example commands to use the `SchemaAnnotation` module:

.. code-block:: bash

    # Annotate schema using default parameters
    SR SchemaAnnotation -s /path/to/schema -o /path/to/output -ao uniprot-proteomes -pt path/to/proteome/table

    # Annotate schema with custom parameters
    SR SchemaAnnotation -s /path/to/schema -o /path/to/output -ao uniprot-proteomes genbank-files -pt path/to/proteome/table -gf path/to/genbank/files -c 4 -t 4 -b 0.7 -tt 1 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `SchemaAnnotation` module, consider the following troubleshooting steps:

- Verify that the paths to the schema and output directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
