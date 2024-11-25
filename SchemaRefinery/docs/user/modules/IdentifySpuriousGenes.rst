IdentifySpuriousGenes - Identify spurious genes in a schema
===========================================================

Description
-----------

The `IdentifySpuriousGenes` module parses command-line arguments and initiates the process to identify spurious genes in a schema. This module sets up an argument parser to handle various command-line options for identifying spurious genes and then calls the main function of the `IdentifySpuriousGenes` class with the parsed arguments.

Features
--------

- Identification of spurious genes in a schema.
- Configurable parameters for the identification process.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `IdentifySpuriousGenes` module can be used as follows:

.. code-block:: bash

    python schema_refinery.py -s /path/to/schema -o /path/to/output -a /path/to/allelecall -pnl /path/to/possible_new_loci -at 0.9 -pt 90 -cs 0.9 -cc 0.9 -gp 10 -as 201 -tt 11 -b 0.6 -sr 0.8 -m schema -pm reps_vs_alleles -c 4 --nocleanup

Command-Line Arguments
----------------------

-s, --schema-directory
    (Required) Path to the created schema directory.

-o, --output-directory
    (Required) Path to the directory to which files will be stored.

-a, --allelecall-directory
    (Required) Path to the directory that contains allele call directory that was run with --no-cleanup.

-pnl, --possible-new-loci
    (Optional) Path to the directory that contains possible new loci.

-at, --alignment_ratio_threshold
    (Optional) Threshold value for alignment used to identify spurious CDS (float: 0-1).
    Default: 0.9

-pt, --pident_threshold
    (Optional) Threshold value for pident values used to identify spurious CDS (int 0-100).
    Default: 90

-cs, --clustering-sim
    (Optional) Similarity value for kmers representatives (float: 0-1).
    Default: 0.9

-cc, --clustering-cov
    (Optional) Coverage value for kmers representatives (float: 0-1).
    Default: 0.9

-gp, --genome_presence
    (Optional) The minimum number of genomes specific cluster of CDS must be present in order to be considered.

-as, --absolute_size
    (Optional) Size of the CDS to consider processing.
    Default: 201

-tt, --translation_table
    (Optional) Translation table to use for the CDS translation.
    Default: 11

-b, --bsr
    (Optional) BSR value to consider alleles as the same locus.
    Default: 0.6

-sr, --size_ratio
    (Optional) Size ratio to consider alleles as the same locus.
    Default: 0.8

-m, --run-mode
    (Optional) Run mode for identifying spurious loci.
    Choices: schema, loci_vs_cds
    Default: schema

-pm, --processing-mode
    (Optional) Mode to run the module.
    Choices: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.
    Default: reps_vs_alleles

-c, --cpu
    (Optional) Number of CPUs to run BLAST instances.
    Default: 1

--nocleanup
    (Optional) Flag to indicate whether to skip cleanup after running the module.

Examples
--------

Here are some example commands to use the `IdentifySpuriousGenes` module:

.. code-block:: bash

    # Identify spurious genes using default parameters
    python schema_refinery.py -s /path/to/schema -o /path/to/output -a /path/to/allelecall

    # Identify spurious genes with custom parameters
    python schema_refinery.py -s /path/to/schema -o /path/to/output -a /path/to/allelecall -pnl /path/to/possible_new_loci -at 0.9 -pt 90 -cs 0.9 -cc 0.9 -gp 10 -as 201 -tt 11 -b 0.6 -sr 0.8 -m schema -pm reps_vs_alleles -c 4 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `IdentifySpuriousGenes` module, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and allele call directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
