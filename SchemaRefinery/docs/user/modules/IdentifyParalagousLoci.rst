IdentifyParalogousLoci - Identify paralogous loci in a schema
=============================================================

Description
-----------

The `IdentifyParalogousLoci` module parses command-line arguments and initiates the process to identify paralogous loci in a schema. This module sets up an argument parser to handle various command-line options for identifying paralogous loci and then calls the main function of the `IdentifyParalogousLoci` class with the parsed arguments.

Features
--------

- Identification of paralogous loci in a schema.
- Configurable parameters for the identification process.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `IdentifyParalogousLoci` module can be used as follows:

.. code-block:: bash

    python schema_refinery.py -s /path/to/schema -o /path/to/output -c 4 -b 0.6 -tt 11 -st 0.2 -pm alleles_vs_alleles --nocleanup

Command-Line Arguments
----------------------

-s, --schema-directory
    (Required) Folder that contains the schema to identify paralogous loci.

-o, --output-directory
    (Required) Path to the directory to which files will be stored.

-c, --cpu
    (Optional) Number of CPUs to run BLAST instances.
    Default: 1

-b, --bsr
    (Optional) BSR value to consider alleles as the same locus.
    Default: 0.6

-tt, --translation_table
    (Optional) Translation table to use for the CDS translation.
    Default: 11

-st, --size_threshold
    (Optional) Size threshold to consider two paralogous loci as similar.
    Default: 0.2

-pm, --processing-mode
    (Optional) Mode to run the module.
    Choices: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.
    Default: alleles_vs_alleles

--nocleanup
    (Optional) Flag to indicate whether to skip cleanup after running the module.

Outputs
-------
Folder and file structure for the output directory of the `IdentifyParalogousLoci` module is shown below. The output directory contains the following files and folders:
.. code-block:: bash
    OutputFolderName
    ├── Blast # --nocleanup
        ├── Blast_db_prot
        ├── Blast_output
            ├── blast_results_x.tsv
            ├── blast_results_y.tsv
            ├── blast_results_z.tsv
            └── ...
        ├── master_file.fasta
        ├── self_score_folder
            ├── blast_results_x.tsv
            ├── blast_results_y.tsv
            ├── blast_results_z.tsv
            └── ...
        └── Translation
            ├── x_translation.fasta
            ├── y_translation.fasta
            ├── z_translation.fasta
            └── ...
    ├── paralogous_loci_report.tsv
    ├── paralogous_loci_report_cluster_by_id.tsv
    └── paralogous_loci_report_passed_all_checks.tsv

Examples
--------

Here are some example commands to use the `IdentifyParalogousLoci` module:

.. code-block:: bash

    # Identify paralogous loci using default parameters
    python schema_refinery.py -s /path/to/schema -o /path/to/output

    # Identify paralogous loci with custom parameters
    python schema_refinery.py -s /path/to/schema -o /path/to/output -c 4 -b 0.7 -tt 1 -st 0.3 -pm reps_vs_reps --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `IdentifyParalogousLoci` module, consider the following troubleshooting steps:

- Verify that the paths to the schema and output directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.