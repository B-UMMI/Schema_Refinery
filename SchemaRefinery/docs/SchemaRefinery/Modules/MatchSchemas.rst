MatchSchemas - Match schemas in a directory
===========================================

Description
-----------

The `MatchSchemas` module parses command-line arguments and initiates the process to match schemas. This module sets up an argument parser to handle various command-line options for matching schemas and then calls the main function of the `MatchSchemas` class with the parsed arguments.

Features
--------

- Configurable parameters for the matching process.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------

- Python 3.9 or higher
- BLAST (`https://www.ncbi.nlm.nih.gov/books/NBK279690/ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- Install requirements using the following command:

.. code-block:: bash

    pip install -r requirements.txt

Usage
-----

The `MatchSchemas` module can be used as follows:

.. code-block:: bash

    SR MatchSchema -fs /path/to/first_schema -ss /path/to/second_schema -o /path/to/output_folder -c 6 -b 0.5 --nocleanup

Command-Line Arguments
----------------------

::

    -fs, --first-schema-directory
        (Required) Folder that contains the first schema to identify paralogous loci.

    -ss, --second-schema-directory
        (Required) Folder that contains the second schema to identify paralogous loci.

    -o, --output-directory
        (Required) Path to the directory to which files will be stored.

    -c, --cpu
        (Optional) Number of CPUs to run BLAST instances.
        Default: 1

    -b, --bsr
        (Optional) BSR value to consider alleles as the same locus.
        Default: 0.6

    -tt, --translation-table
        (Optional) Translation table to use for the CDS translation.
        Default: 11

    -ra, --rep-vs-alleles
        (Optional) If True then after the rep vs rep Blast the program will run a second Blast with rep vs alleles.
        Default: False

    --nocleanup
        (Optional) Flag to indicate whether to skip cleanup after running the module.

    --debug
        (Optional) Flag to indicate whether to run the module in debug mode.
        Default: False

    --logger
        (Optional) Path to the logger file.
        Default: None

Algorithm Explanation
---------------------

MatchSchemas Algorithm Flowchart:

.. image:: source/match_schemas.png
   :alt: MatchSchemas Algorithm
   :width: 80%
   :align: center

The algorithm goes throught three matching process and a fourth optional one.
    - DNA hashes comparision
    - Protein hashes comparision
    - Blast with reps vs reps
    - (optional) Blast with rep vs alleles

Outputs
-------
Folder and file structure for the output directory of the `MatchSchemas` module is shown below. The output directory contains the following files and folders:

::

    OutputFolderName
    ├── blast_processing # --nocleanup
    |   ├── subject_master_file.fasta
    |   ├── subject_master_rep_file.fasta
    │   ├── Blast_db_prot
    │   │   ├── Blast_db_protein.pdb
    │   │   ├── Blast_db_protein.phr
    │   │   ├── Blast_db_protein.pin
    │   │   ├── Blast_db_protein.pog
    │   │   ├── Blast_db_protein.pos
    │   │   ├── Blast_db_protein.pot
    │   │   ├── Blast_db_protein.psq
    │   │   ├── Blast_db_protein.ptf
    │   │   └── Blast_db_protein.pto
    │   ├── blastp_results
    │   │   ├── blast_results_x.tsv
    │   │   ├── blast_results_y.tsv
    │   │   ├── blast_results_z.tsv
    │   │   └── ...
    │   ├── Query_Translation
    │   │   ├── x_translation.fasta
    │   │   ├── y_translation.fasta
    │   │   ├── z_translation.fasta
    │   │   └── ...
    |   ├── Query_Translation_Rep
    │   │   ├── x_rep_translation.fasta
    │   │   ├── y_rep_translation.fasta
    │   │   ├── z_rep_translation.fasta
    │   │   └── ...
    |   ├── Query_Not_Translated
    │   │   ├── x_not_translated.fasta
    │   │   ├── y_not_translated.fasta
    │   │   ├── z_not_translated.fasta
    │   │   └── ...
    |   ├── Query_Not_Translated_Rep
    │   │   ├── x_rep_not_translated.fasta
    │   │   ├── y_rep_not_translated.fasta
    │   │   ├── z_rep_not_translated.fasta
    │   │   └── ...
    │   ├── Subject_Translation
    │   │   ├── x_translation.fasta
    │   │   ├── y_translation.fasta
    │   │   ├── z_translation.fasta
    │   │   └── ...
    |   ├── Subject_Translation_Rep
    │   │   ├── x_rep_translation.fasta
    │   │   ├── y_rep_translation.fasta
    │   │   ├── z_rep_translation.fasta
    │   │   └── ...
    │   └── self_score_folder
    │       ├── blast_results_x.tsv
    │       ├── blast_results_y.tsv
    │       ├── blast_results_z.tsv
    │       └── ...
    ├── existing_matches.txt
    └── Match_Schemas_Results.tsv

.. toctree::
   :maxdepth: 1

   MatchSchemasReportExplanation

Report files description
------------------------

.. csv-table:: **Match_Schemas_Results.tsv**
    :header: "Query", "Subject", "BSR", "Process"
    :widths: 30, 30, 10, 30

    x, y, 0.8, hashes_vs_hashes
    z, a, 1.0, hashes_vs_hashes
    b, c, 0.965, rep_vs_rep
    d, Not Matched, NA, rep_vs_rep
    Not Matched, e, NA, rep_vs_alleles
    ...

columns description:

::
    
    Query: The locus from the query schema.
    Subject: The best match for the locus from the subject schema.
    BSR: The BSR value for the best match.
    Process: From which type of comparision is this match from.

For the hashes comparisions the BSR value will always be set to 1.0 as this will be considered the same locus.

Examples
--------

Here are some example commands to use the `MatchSchemas` module:

.. code-block:: bash

    # Match schemas using default parameters
    SR MatchSchemas -fs /path/to/first_schema -ss /path/to/second_schema -o /path/to/output_folder

    # Match schemas with custom parameters
    SR MatchSchemas -fs /path/to/first_schema -ss /path/to/second_schema -o /path/to/output_folder -ra True -c 4 -b 0.7 -tt 4 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `MatchSchemas` module, consider the following troubleshooting steps:

- Verify that the paths to the query and subject schema directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
