.. _IdentifyParalogousLoci:

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

- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `IdentifyParalogousLoci` module can be used as follows:

.. code-block:: bash

    SR IdentifyParalogousLoci -s /path/to/schema -o /path/to/output -c 4 -b 0.6 -tt 11 -st 0.2 -pm alleles_vs_alleles --nocleanup

Command-Line Arguments
----------------------
::

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

Algorithm Explanation
---------------------

Outputs
-------
Folder and file structure for the output directory of the `IdentifyParalogousLoci` module is shown below. The output directory contains the following files and folders:

::

    OutputFolderName
    ├── Blast # --nocleanup
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
    │   ├── Blast_output
    │   │   ├── blast_results_x.tsv
    │   │   ├── blast_results_y.tsv
    │   │   ├── blast_results_z.tsv
    │   │   └── ...
    │   ├── master_file.fasta
    │   ├── self_score_folder
    │   │   ├── blast_results_x.tsv
    │   │   ├── blast_results_y.tsv
    │   │   ├── blast_results_z.tsv
    │   │   └── ...
    │   └── Translation
    │       ├── x_translation.fasta
    │       ├── y_translation.fasta
    │       ├── z_translation.fasta
    │       └── ...
    ├── paralogous_loci_report.tsv
    ├── paralogous_loci_report_cluster_by_id.tsv
    └── paralogous_loci_report_passed_all_checks.tsv

Output files and folders description:
-------------------------------------

**OutputFolderName**: The folder where the output files are stored.

    Blast: Folder containing BLASTp database, BLASTp output files, master file, self-score folder, and translation files.
        Blast_db_prot: Folder containing the BLASTp database.
            Blast_db_protein.pdb: Position-specific Data Base file. Contains position-specific scoring matrices (PSSMs) used in PSI-BLAST searches.
            Blast_db_protein.phr: Protein Header Record file. Contains the header information for each sequence in the protein database.
            Blast_db_protein.pin: Protein Index file. Contains the index of the sequences in the protein database.
            Blast_db_protein.pog: Protein Organism Group file. Contains information about the taxonomic grouping of the sequences in the protein database.
            Blast_db_protein.pos: Protein Organism Sequence file. Contains the actual sequence data for the protein database.
            Blast_db_protein.pot: Protein Organism Taxonomy file. Contains taxonomic information for the sequences in the protein database.
            Blast_db_protein.psq: Protein Sequence Query file. Contains the sequence data in a format optimized for BLAST searches.
            Blast_db_protein.ptf: Protein Taxonomy File. Contains taxonomy information for the sequences in the protein database.
            Blast_db_protein.pto: Protein Taxonomy Organism file. Contains organism-specific taxonomy information for the sequences in the protein database.
        Blast_output: Folder containing the BLASTp output files.
            blast_results_x.tsv: TSV file containing the BLASTp results for the locus x.
            blast_results_y.tsv: TSV file containing the BLASTp results for the locus y.
            blast_results_z.tsv: TSV file containing the BLASTp results for the locus z.
            ...: All of the other TSV BLASTp results files.
        master_file.fasta: FASTA file containing all of the protein sequences used in the analysis (used to create BLAST DB).
        self_score_folder: Folder containing the self-score BLAST results.
            blast_results_x.tsv: TSV file containing the BLASTp results for self-score for the locus x.
            blast_results_y.tsv: TSV file containing the BLASTp results for self-score for the locus y.
            blast_results_z.tsv: TSV file containing the BLASTp results for self-score for the locus z.
            ...: All of the other TSV BLASTp for self-score results files.
        Translation: Folder containing the translation files.
            x_translation.fasta: FASTA file containing the translation for the locus x.
            y_translation.fasta: FASTA file containing the translation for the locus y.
            z_translation.fasta: FASTA file containing the translation for the locus z.
            ...: All of the other translation files.

    **paralogous_loci_report.tsv**: TSV file containing the report of the paralogous loci.
    **paralogous_loci_report_cluster_by_id.tsv**: TSV file containing the report of the paralogous loci clustered by ID.
    **paralogous_loci_report_passed_all_checks.tsv**: TSV file containing the report of the paralogous loci clustered by ID that passed all checks.

Report files description
------------------------

.. csv-table:: **paralogous_loci_report.tsv**
   :header: "Query_loci_id", "Subject_loci_id", "BSR", "if_loci_intersect", "if_close_distance", "Loci_min_allele_size", "Loci_max_allele_size", "Loci_mode_allele_size", "Loci_mean_allele_size"
   :widths: 20, 20, 10, 10, 10, 20, 20, 20, 20

   x, a, 0.7360028348688873, False, False, 416.0|1349.0, 544.0|1645.0, 515.0|1628.0, 476.75|1615.125
   x, b, 0.6146651702207258, False, True, 416.0|599.0, 544.0|738.0, 515.0|738.0, 476.75|720.625
   x, c, 0.6523642732049036, True, True, 416.0|466.0, 544.0|547.0, 515.0|547.0, 476.75|512.5714285714286
   ...

columns description:

::
    Query_loci_id: The locus ID of the query.
    Subject_loci_id: The locus ID of the subject.
    BSR: The BSR value between the query and subject loci.
    if_loci_intersect: If the loci intersect.
    if_close_distance: If the loci are close in distance.
    Loci_min_allele_size: The minimum allele size of the loci, query and subject values are separated by '|'.
    Loci_max_allele_size: The maximum allele size of the loci, query and subject values are separated by '|'.
    Loci_mode_allele_size: The mode allele size of the loci, query and subject values are separated by '|'.
    Loci_mean_allele_size: The mean allele size of the loci, query and subject values are separated by '|'.

.. csv-table:: **paralogous_loci_report_cluster_by_id.tsv**
   :header: "Joined_loci_id", "Clustered_loci_ids"
   :widths: 20, 80

   x, x,a,b,c
   y, y,d,e,f
   z, z,g,h
   ...

columns description:

.. csv-table:: **paralogous_loci_report_passed_all_checks.tsv**
   :header: "Joined_loci_id", "Clustered_loci_ids"
   :widths: 20, 80

   x, x,b,c
   y, y,d
   z, z,h
   ...

Examples
--------

Here are some example commands to use the `IdentifyParalogousLoci` module:

.. code-block:: bash

    # Identify paralogous loci using default parameters
    SR IdentifyParalogousLoci -s /path/to/schema -o /path/to/output

    # Identify paralogous loci with custom parameters
    SR IdentifyParalogousLoci -s /path/to/schema -o /path/to/output -c 4 -b 0.7 -tt 4 -st 0.3 -pm reps_vs_reps --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `IdentifyParalogousLoci` module, consider the following troubleshooting steps:

- Verify that the paths to the schema and output directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.