IdentifySpuriousGenes - Identify spurious genes in a schema
===========================================================

Description
-----------

The `IdentifySpuriousGenes` module parses command-line arguments and initiates the process to identify spurious genes in a schema. This module sets up an argument parser to handle various command-line options for identifying spurious genes and then calls the main function of the `IdentifySpuriousGenes` class with the parsed arguments. This module is essential for researchers and bioinformaticians who need to detect and analyze spurious loci in order to create a more refined and concise final schema.

This module takes unclassified CDS or schema loci and matches them against each other, providing a classification for each match. Based on the classification, the best class is chosen to represent the relationship between the two loci. Using these relationships, the user can select the appropriate loci or unclassified CDS group to be included/joined in the final schema removing those that are spurious.

Features
--------

- Identification of spurious genes in a schema.
- Optional annotation of the recommendations.
- Configurable parameters for the identification process.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------

- Python between 3.9 and 3.11
- BLAST (`https://www.ncbi.nlm.nih.gov/books/NBK279690/ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- ChewBBACA (https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html or using bioconda)
- Install requirements using the following command:

.. code-block:: bash

    pip install -r requirements.txt

Usage
-----

The `IdentifySpuriousGenes` module can be used as follows:

.. code-block:: bash

    SR IdentifySpuriousGenes -s /path/to/schema -o /path/to/output -a /path/to/allelecall -tt 11 -b 0.6 -m schema -pm reps_vs_alleles -c 4 --nocleanup


Command-Line Arguments
----------------------

::

    -s, --schema-directory
        (Required) Path to the created schema directory.

    -o, --output-directory
        (Required) Path to the directory to which files will be stored.

    -a, --allelecall-directory
        (Required) Path to the directory that contains allele call directory from chewBBACA that was run with --no-cleanup and --output-unclassified.

    -ann, --annotations
        (Optional) Path to the tsv file with the schema annotations.
        This file needs to have one column with loci with the same IDs as the ones in the schema.

    -at, --alignment-ratio-threshold
        (Optional) Threshold value for alignment used to identify spurious CDS (float: 0-1).
        Default: 0.9

    -pt, --pident-threshold
        (Optional) Threshold value for pident values used to identify spurious CDS (int 0-100).
        Default: 90

    -cs, --clustering-sim
        (Optional) Similarity value for kmers representatives (float: 0-1).
        Default: 0.9

    -cc, --clustering-cov
        (Optional) Coverage value for kmers representatives (float: 0-1).
        Default: 0.9

    -gp, --genome-presence
        (Optional) The minimum number of genomes specific cluster of CDS must be present in order to be considered.

    -as, --absolute-size
        (Optional) Size of the CDS to consider processing.
        Default: 201

    -tt, --translation-table
        (Optional) Translation table to use for the CDS translation.
        Default: 11

    -b, --bsr
        (Optional) BSR value to consider alleles as the same locus.
        Default: 0.6

    -sr, --size-ratio
        (Optional) Size ratio to consider alleles as the same locus.
        Default: 0.8

    -m, --run-mode
        (Optional) Run mode for identifying spurious loci.
        Choices: unclassified_cds, schema, schema_vs_schema
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

    --debug
        (Optional) Flag to indicate whether to run the module in debug mode.
        Default: False

    --logger
        (Optional) Path to the logger file.
        Default: None

The `--allelecall-directory` argument must be the folder obtained using the chewBBACA module AlleleCall with the arguments --no-cleanup and --output-unclassified. Without these arguments the needed files to run the `IdentifySpuriousGenes` module will not run.

The options `--schema-directory` and `--allelecall-directory` must have two paths each if the `run_mode` is `schema_vs_schema`. One path per schema and each schema must have a respective AlleleCall directory.

.. Note::
    In `processing_mode` the option `reps_vs_reps` is the fastest and covers most of the cases.
    The option `alleles_vs_alleles` takes much more time and changes the recommendation of around more 1% of the total loci. 
    Choose the processing mode according to your computer capabilities and research needs.

Algorithm Explanation
---------------------

Algorithm to identify new loci based on the CDS that are not in the schema:

.. image:: source/IdentifySpuriousGenes_unclassifiedCDS.png
   :alt: Algorithm for unclassified CDS
   :width: 80%
   :align: center

Algorithm to indentify spurious loci based on schema inputs for run modes schema and schema_vs_schema:

.. image:: source/IdentifySpuriousGenes_schema.png
   :alt: Algorithm to identify spurious loci
   :width: 80%
   :align: center


This classification algorithm goes through 2 rounds of BLAST. First a BLASTp is done, and the results filtered based on BSR and pident thresholds. The loci that are not matched are then passed to a BLASTn processing. By running the BLASTn after the BLASTp we ensure that valid CDS are still caught and matched regardless if they contain frameshifts or are pseudo or partial genes. In this way we obtain more complete results and matches. 

The BLAST tool has a limit of hits written per loci in the output file. Using a value for this limit that will not overload the space alocated for our tool, it can happen that not all matches are written down if the loci have a high number of alleles. For that reason, all BLAST processes have the loci being aligned as a query, including all its alleles, removed from the databased used for that run. This is done to ensure that for all loci that will have a match, that match won't be with itself.

Each BLAST results are parsed and given a class based on the following rules:

.. image:: source/algorithm_classification.png
   :alt: Classification algorithm
   :width: 80%
   :align: center

---------------------------------------------------------------------------------

Between the two loci, the best class is chosen based on the following order of the classes to represent the relationship between the two loci.

classification order: 1a, 1b, 2a, 3a, 2b, 1c, 3b, 4a, 4b, 4c, 5, 6, 7

Depending on the classification, each locus will have a specific action recommended:
::
    Join: 1a
    Choice: 1c, 2b, 3b, 4b, 6
    Drop: 1b, 2a, 3a, 4a
    Add: 4c, 5, 7


Choice is a recommendation for the user to decide if the locus should just be added with no alteration ("Add"), dropped ("Drop") or joined to the cluster it is in ("Join"). 

The other loci will have the action "Add" and will just be added, without alteration, to the new schema.

All matches resukting from the BLASTn will be classified as 6 which corresponds to the action "Choice". As a BLASTn does not involve BSR and frequency values it can not be classified as "Drop" or "Add" from the algorithm. However, for this file to be used in the `CreateSchemaStructure` module these "Choice" actions will have to be change into one of the other actions. This implies that the user must review and annotate the recommendations anc choose the action it finds best.

The recommendation associated to each locus can change in the final output. For example, a locus that was part of a 1a cluster can be marked as "Drop" if it does not pass certain thresholds of frequency.

The annotation option will use the `consolidate` mode from the `SchemaAnnoation` module, so the input format should comform with the rules set in the SchemaAnnotation documentation.


Outputs
-------
Folder and file structure for the output directory of the `IdentifySpuriousGenes` module is shown below. The output directory contains the following files and folders:

Since there are two run modes, the output directory structure will vary based on the run mode selected.

**For --run-mode schema and schema_vs_schema:**

::

    OutputFolderName
    ├── 1_schema_processing # --nocleanup
    │   ├── master_nucleotide.fasta
    │   ├── schema
    │   │   ├── .genes_list
    │   │   ├── .schema_config
    │   │   ├── loci_x.fasta
    │   │   ├── new_loci_y.fasta
    │   │   ├── ...
    │   │   └── short
    │   │       ├── loci_x.fasta
    │   │       ├── new_loci_y.fasta
    │   │       └── ...
    │   └── schema_translation_folder
    │       ├── loci_x.fasta
    │       ├── new_loci_y.fasta
    │       └── ...
    ├── 2_BLAST_processing # --nocleanup
    │   ├── 1_BLASTn_processing
    │   │   ├── blast_db_nucl
    │   │   │   ├── Blast_db_nucleotide.ndb
    │   │   │   ├── Blast_db_nucleotide.nhr
    │   │   │   ├── Blast_db_nucleotide.nin
    │   │   │   ├── Blast_db_nucleotide.nog
    │   │   │   ├── Blast_db_nucleotide.nsd
    │   │   │   ├── Blast_db_nucleotide.nsi
    │   │   │   ├── Blast_db_nucleotide.nsq
    │   │   │   ├── Blast_db_nucleotide.ntf
    │   │   │   └── Blast_db_nucleotide.nto
    │   │   └── BLASTn_results
    │   │       ├── blast_results_x.tsv
    │   │       ├── blast_results_y.tsv
    │   │       └── ...
    │   └── 2_BLASTp_processing
    │       ├── blastn_results_matches_translations
    │       │   │   └── cluster_rep_translation
    │       │   │       ├── cluster_rep_translation_x.fasta
    │       │   │       ├── cluster_rep_translation_y.fasta
    │       │   │       └── ...
    │       │   ├── cluster_matches_translation_x.tsv
    │       │   ├── cluster_matches_translation_y.tsv
    │       │   └── ...
    │       ├── BLASTp_results
    │       │   ├── blast_results_x.tsv
    │       │   ├── blast_results_y.tsv
    │       │   └── ...
    │       └── BLASTp_results_self_score_results
    │           ├── blast_results_x.tsv
    │           ├── blast_results_y.tsv
    │           └── ...
    ├── 3_processing_results # --nocleanup
    │   └── blast_results
    │       ├── blast_all_matches.tsv
    │       ├── blast_by_cluster
    │       │   ├── blast_joined_cluster_x.tsv
    │       │   ├── blast_retained_y.tsv
    │       │   └── ...
    │       └── blast_results_by_class
    │           ├── class_1a.tsv
    │           ├── class_2a.tsv
    │           └── ...
    ├── count_results_by_cluster.tsv
    ├── drop_loci_reason.tsv
    ├── recommendations.tsv
    ├── recommendations_annotations.tsv # -ann
    └── related_matches.tsv

**For --run-mode unclassified_cds:**

::

    OutputFolderName
    ├── 1_CDS_processing # --nocleanup
    │   ├── CDS_not_found.fasta
    │   └── CDS_not_found_translation.fasta
    ├── 2_BLAST_processing # --nocleanup
    │   ├── 1_BLASTn_processing
    │   │   ├── cluster_representatives_fastas_dna
    │   │   │   ├── cluster_rep_x.fasta
    │   │   │   ├── cluster_rep_y.fasta
    │   │   │   └── ...
    │   │   ├── Blast_db_nucl
    │   │   │   ├── Blast_db_nucleotide.ndb
    │   │   │   ├── Blast_db_nucleotide.nhr
    │   │   │   ├── Blast_db_nucleotide.nin
    │   │   │   ├── Blast_db_nucleotide.nog
    │   │   │   ├── Blast_db_nucleotide.nsd
    │   │   │   ├── Blast_db_nucleotide.nsi
    │   │   │   ├── Blast_db_nucleotide.nsq
    │   │   │   ├── Blast_db_nucleotide.ntf
    │   │   │   └── Blast_db_nucleotide.nto
    │   │   └── BLASTn_results
    │   │       ├── blast_results_x.tsv
    │   │       ├── blast_results_y.tsv
    │   │       └── ...
    │   └── 2_BLASTp_processing
    │       ├── blastn_results_matches_translations
    │       │   │   └── cluster_rep_translation
    │       │   │       ├── cluster_rep_translation_x.fasta
    │       │   │       ├── cluster_rep_translation_y.fasta
    │       │   │       └── ...
    │       │   ├── cluster_matches_translation_x.tsv
    │       │   ├── cluster_matches_translation_y.tsv
    │       │   └── ...
    │       ├── BLASTp_results
    │       │   ├── blast_results_x.tsv
    │       │   ├── blast_results_y.tsv
    │       │   └── ...
    │       └── BLASTp_results_self_score_results
    │           ├── blast_results_x.tsv
    │           ├── blast_results_y.tsv
    │           └── ...
    ├── 3_processing_results # --nocleanup
    │   ├── blast_results
    │   │   ├── blast_all_matches.tsv
    │   │   ├── blast_by_cluster
    │   |   │   ├── blast_joined_cluster_x.tsv
    │   |   │   ├── blast_retained_y.tsv
    │   │   │   └── ...
    │   │   └── blast_results_by_class
    │   │       ├── class_1a.tsv
    │   │       ├── class_2a.tsv
    │   │       └── ...
    │   ├── cds_id_changes.tsv
    │   ├── dropped_cds.tsv
    │   └── Graph_folder
    │       ├── All_of_CDS_graphs.html
    │       ├── graphs_class_1a.html
    │       └── ...
    ├── count_results_by_cluster.tsv
    ├── drop_loci_reason.tsv
    ├── recommendations.tsv
    ├── recommendations_annotations.tsv # -ann
    ├── related_matches.tsv
    └── temp_fastas
        ├── x.fasta
        ├── y.fasta
        └── ...
  

.. toctree::
   :maxdepth: 1

   IdentifySpuriousGenesOutputDescription


Report files description
------------------------

.. csv-table:: **count_results_by_cluster.tsv**
    :header: "Query", "Subject", "1a", "1b", "2a", "3a", "2b", "1c", "3b", "4a", "4b", "4c", "5", "Representatives_count", "Alelles_count", "Frequency_in_genomes_query", "Frequency_in_genomes_subject"
    :widths: 15, 15, 20, 5, 5, 5, 5, 15, 5, 5, 5, 5, 5, 20, 20, 25, 25

    x, y, 378|1024|-|1024, -, -, -, -, 646|1024|1024|1024, -, -, -, -, -, 16|64, 16|64, 223, 133
    x, z, -, -, -, -, -, 128|128|128|128, -, -, -, -, -, 16|8, 16|8, 223, 99
    x, w, 6|224|1|224, -, -, -, -, 218|224|223|224, -, -, -, -, -, 16|14, 16|14, 223, 221
    ...

Columns description:
::

    Query: The query locus.
    Subject: The subject locus.
    1a-5: The count of the loci in the cluster, interpret the values as this, for x query and y subject class 1a '378|1024|-|1024', x has 378 matches out of 1024 to y that are class 1a and while y has no matches '-' out of 1024 to x.
    alleles_used_to_blast_count: The count of alleles used to blast.
    alleles_blasted_against_count: The count of alleles blasted against.
    Frequency_in_genomes_query: The frequency of the query locus in genomes.
    Frequency_in_genomes_subject: The frequency of the subject locus in genomes.

.. csv-table:: **drop_loci_reason.tsv**
    :header: "Possible_new_loci_ID", "Drop_Reason"
    :widths: 30, 60

    x, Dropped_due_to_smaller_genome_presence_than_matched_cluster
    y, Dropped_due_to_smaller_genome_presence_than_matched_cluster
    z, Dropped_due_to_smaller_genome_presence_than_matched_cluster
    ...

Columns description:

::

    Possible_new_loci_ID: The identifier for the possible new locus.
    Drop_Reason: The reason for dropping the locus.

.. csv-table:: **recommendations.tsv**
    :header: "Locus", "Action", "Class"
    :widths: 15, 20

    x, Join, 1a
    y, Join, 1a
    z, Choice, 3b
    #
    a, Choice, 1c
    b, Choice, 1c
    #
    c, Drop, 1b
    #
    ...

Columns description:

::

    Locus: Name of the locus to be joined in the clustered.
    Action: Action to be taken (Join, Choice, Drop or Add).
    #: Separates each cluster of loci.

This is the main output file as it can be used as the input of the `CreateSchemaStructure`.
It is adviced to annotate this file usinf the argument `--annotations` in order to more easily review and confirm the resulting clusters and actions. For this process the user can consult the previews files as well as proteome and functional annotations.

Some loci are in single clusters and with action "Drop". In these cases these loci were inserted in a cluster with one other loci both marked "Join" but due to low frequency this one will be dropped, and the other loci will be changed into "Add". To check the full cluster relationships consult the `related_matches.tsv` file.

.. Note:: 
    Before passing the recommendation file to the `CreateSchemaStructure` module make sure there are no "Choice" in the action column. After review the clusters the user should change these actions into "Join", "Drop" or "Add".



.. csv-table:: **recommendations_annotations.tsv**
   :header: "Loci", "Action", "Locus_annotation", "Annotation"
   :widths: 15, 20, 15, 40

   x, Join, x, annotation
   y, Join, y, annotation
   z, Choice, z, annotation
   #
   a, Choice, a, annotation
   b, Choice, b, annotation
   #
   c, Drop, c, annotation
   #
   ...

Columns description:

::

    Loci: Name of the locus to be joined in the clustered.
    Action: Action to be taken (Join, Choice, Drop or Add).
    Locus_annotation: Name of the Locus in the annotation file (should be the same as Locus).
    Annotation: Column with annotation from the annotation file.
        (The name and numer of columns with annotations will depend on the file with the annotation, follows the structure of the output of the consolidate module).
    #: Separates each cluster of loci.


.. csv-table:: **related_matches.tsv**
    :header: "Query", "Subject", "Class", "Class_count", "Inverse_class", "Inverse_class_count", "Frequency_in_genomes_query", "Frequency_in_genomes_subject", "alleles_used_to_blast_count", "alleles_blasted_against_count"
    :widths: 20, 20, 10, 10, 10, 10, 20, 20, 20, 20

    x, y, 1a, 378/1024, 1c, 1024/1024, 223, 133, 16|64, 16|64
    x, z, 1c, 128/128, 1c, 128/128, 223, 99, 16|8, 16|8
    x, w, 1a, 6/224, 1a, 1/224, 223, 221, 16|14, 16|14
    #
    a, b, 1a, 378/1024, 1c, 1024/1024, 223, 133, 16|64, 16|64
    a, c, 1c, 128/128, 1c, 128/128, 223, 99, 16|8, 16|8
    ...
    
Columns description:

::
    
    Query: The query locus.
    Subject: The subject locus.
    Class: The best class of the matches for that loci.
    Class_count: The count of matches of the loci in the class.
    Inverse_class: The best class of the inverse match for those loci.
    Inverse_class_count: The count of inverse matches of the loci in that class.
    Frequency_in_genomes_query: The frequency of the query locus in genomes.
    Frequency_in_genomes_subject: The frequency of the subject locus in genomes.
    alleles_used_to_blast_count: The count of alleles used to blast.
    alleles_blasted_against_count: The count of alleles blasted against.

This file contains a simple overview of the clusters and the pairs within them alongside some important values for later consultation, such as the frquency of each locus. 

Examples
--------

Here are some example commands to use the `IdentifySpuriousGenes` module:

.. code-block:: bash

    # Identify spurious genes using default parameters
    SR IdentifySpuriousGenes -s /path/to/schema -o /path/to/output -a /path/to/allelecall

    # Identify spurious genes with custom parameters
    SR IdentifySpuriousGenes -s /path/to/schema /path/to/second_schema -o /path/to/output -a /path/to/allelecall /path/to/allelecall_second_schema  -ann /path/to/annotation-files -at 0.9 -pt 90 -cs 0.9 -cc 0.9 -gp 10 -as 201 -tt 11 -b 0.6 -sr 0.8 -m schema_vs_schema -pm reps_vs_alleles -c 4 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `IdentifySpuriousGenes` module, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and allele call directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
