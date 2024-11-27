.. _IdentifySpuriousGenes:

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

- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `IdentifySpuriousGenes` module can be used as follows:

.. code-block:: bash

    SR IdentifySpuriousGenes -s /path/to/schema -o /path/to/output -a /path/to/allelecall -pnl /path/to/possible_new_loci -at 0.9 -pt 90 -cs 0.9 -cc 0.9 -gp 10 -as 201 -tt 11 -b 0.6 -sr 0.8 -m schema -pm reps_vs_alleles -c 4 --nocleanup

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
    Choices: unclassified_cds, schema
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

Outputs
-------
Folder and file structure for the output directory of the `IdentifySpuriousGenes` module is shown below. The output directory contains the following files and folders:

Since there are two run modes, the output directory structure will vary based on the run mode selected.

**For --run-mode schema:**

::

    OutputFolderName
    ├── 1_schema_processing # --nocleanup
    │   ├── master.fasta
    │   ├── schema
    │   │   ├── loci_x.fasta
    │   │   ├── new_loci_y.fasta
    │   │   ├── ...
    │   │   ├── loci_modes
    │   │   ├── species.trn
    │   │   ├── pre_computed
    │   │   │   ├── DNAtable1
    │   │   │   └── PROTEINtable1
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
    │   │   │   ├── cluster_x.tsv
    │   │   │   ├── cluster_y.tsv
    │   │   │   └── ...
    │   │   └── blast_results_by_class
    │   │       ├── class_1.tsv
    │   │       ├── class_2.tsv
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
    └── related_matches.tsv

**For --run-mode unclassified_cds:**

::

    OutputFolderName
    ├── 1_CDS_processing # --nocleanup
    │   ├── CDS_not_found.fasta
    │   └── CDS_not_found_translation.fasta
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
    │   │   │   ├── cluster_x.tsv
    │   │   │   ├── cluster_y.tsv
    │   │   │   └── ...
    │   │   └── blast_results_by_class
    │   │       ├── class_1.tsv
    │   │       ├── class_2.tsv
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
    ├── related_matches.tsv
    ├── temp_fastas
    │   ├── cluster_x.fasta
    │   ├── cluster_y.fasta
    │   └── ...
    └── temp_fastas_path.txt

Output files and folders description:
-------------------------------------
**For --run-mode schema:**

OutputFolderName: The folder where the output files are stored.

    1_schema_processing: Folder containing schema processing results.
        master.fasta: Master FASTA file.
        schema: Folder containing schema files.
            loci_x.fasta: FASTA file for locus x.
            new_loci_y.fasta: FASTA file for new locus y.
            ...: Other loci files.
            loci_modes: Folder containing loci modes.
            species.trn: Species translation file.
            pre_computed: Folder containing pre-computed data.
                DNAtable1: Pre-computed DNA table.
                PROTEINtable1: Pre-computed protein table.
            short: Folder containing short loci files.
                loci_x.fasta: Short FASTA file for locus x.
                new_loci_y.fasta: Short FASTA file for new locus y.
                ...: Other short loci files.
        schema_translation_folder: Folder containing schema translations.
            loci_x.fasta: Translation for locus x.
            new_loci_y.fasta: Translation for new locus y.
            ...: Other translations.

    2_BLAST_processing: Folder containing BLAST processing results.
        1_BLASTn_processing: Folder containing BLASTn processing results.
            blast_db_nucl: Folder containing BLASTn database.
                Blast_db_nucleotide.ndb: BLASTn nucleotide database file.
                Blast_db_nucleotide.nhr: BLASTn nucleotide header file.
                Blast_db_nucleotide.nin: BLASTn nucleotide index file.
                Blast_db_nucleotide.nog: BLASTn nucleotide organism group file.
                Blast_db_nucleotide.nsd: BLASTn nucleotide sequence data file.
                Blast_db_nucleotide.nsi: BLASTn nucleotide sequence index file.
                Blast_db_nucleotide.nsq: BLASTn nucleotide sequence query file.
                Blast_db_nucleotide.ntf: BLASTn nucleotide taxonomy file.
                Blast_db_nucleotide.nto: BLASTn nucleotide taxonomy organism file.
            BLASTn_results: Folder containing BLASTn results.
                blast_results_x.tsv: BLASTn results for x.
                blast_results_y.tsv: BLASTn results for y.
                ...: Other BLASTn results.
        2_BLASTp_processing: Folder containing BLASTp processing results.
            blastn_results_matches_translations: Folder containing BLASTn results matches translations.
                cluster_matches_translation_x.tsv: Cluster matches translation for x.
                cluster_matches_translation_y.tsv: Cluster matches translation for y.
                ...: Other cluster matches translations.
            BLASTp_results: Folder containing BLASTp results.
                blast_results_x.tsv: BLASTp results for x.
                blast_results_y.tsv: BLASTp results for y.
                ...: Other BLASTp results.
            BLASTp_results_self_score_results: Folder containing BLASTp self-score results.
                blast_results_x.tsv: BLASTp self-score results for x.
                blast_results_y.tsv: BLASTp self-score results for y.
                ...: Other BLASTp self-score results.

    3_processing_results: Folder containing processing results.
        blast_results: Folder containing BLAST results.
            blast_all_matches.tsv: TSV file containing all BLAST matches.
            blast_by_cluster: Folder containing BLAST results by cluster.
                cluster_x.tsv: BLAST results for cluster x.
                cluster_y.tsv: BLAST results for cluster y.
                ...: Other cluster results.
            blast_results_by_class: Folder containing BLAST results by class.
                class_1.tsv: BLAST results for class 1.
                class_2.tsv: BLAST results for class 2.
                ...: Other class results.
        cds_id_changes.tsv: TSV file containing changes in CDS IDs.
        dropped_cds.tsv: TSV file containing dropped CDS.
        Graph_folder: Folder containing graphs.
            All_of_CDS_graphs.html: HTML file containing all CDS graphs.
            graphs_class_1a.html: HTML file containing class 1a graphs.
            ...: Other graph files.

    **count_results_by_cluster.tsv**: TSV file containing count results by cluster.
    **drop_loci_reason.tsv**: TSV file containing reasons for dropping loci.
    **recommendations.tsv**: TSV file containing recommendations.
    **related_matches.tsv**: TSV file containing related matches.

**For --run-mode unclassified_cds:**

OutputFolderName: The folder where the output files are stored.

    1_CDS_processing: Folder containing CDS processing results.
        CDS_not_found.fasta: FASTA file containing CDS not found.
        CDS_not_found_translation.fasta: FASTA file containing translations of CDS not found.

    2_BLAST_processing: Folder containing BLAST processing results.
        1_BLASTn_processing: Folder containing BLASTn processing results.
            blast_db_nucl: Folder containing BLASTn database.
                Blast_db_nucleotide.ndb: BLASTn nucleotide database file.
                Blast_db_nucleotide.nhr: BLASTn nucleotide header file.
                Blast_db_nucleotide.nin: BLASTn nucleotide index file.
                Blast_db_nucleotide.nog: BLASTn nucleotide organism group file.
                Blast_db_nucleotide.nsd: BLASTn nucleotide sequence data file.
                Blast_db_nucleotide.nsi: BLASTn nucleotide sequence index file.
                Blast_db_nucleotide.nsq: BLASTn nucleotide sequence query file.
                Blast_db_nucleotide.ntf: BLASTn nucleotide taxonomy file.
                Blast_db_nucleotide.nto: BLASTn nucleotide taxonomy organism file.
            BLASTn_results: Folder containing BLASTn results.
                blast_results_x.tsv: BLASTn results for x.
                blast_results_y.tsv: BLASTn results for y.
                ...: Other BLASTn results.
        2_BLASTp_processing: Folder containing BLASTp processing results.
            blastn_results_matches_translations: Folder containing BLASTn results matches translations.
                cluster_matches_translation_x.tsv: Cluster matches translation for x.
                cluster_matches_translation_y.tsv: Cluster matches translation for y.
                ...: Other cluster matches translations.
            BLASTp_results: Folder containing BLASTp results.
                blast_results_x.tsv: BLASTp results for x.
                blast_results_y.tsv: BLASTp results for y.
                ...: Other BLASTp results.
            BLASTp_results_self_score_results: Folder containing BLASTp self-score results.
                blast_results_x.tsv: BLASTp self-score results for x.
                blast_results_y.tsv: BLASTp self-score results for y.
                ...: Other BLASTp self-score results.

    3_processing_results: Folder containing processing results.
        blast_results: Folder containing BLAST results.
            blast_all_matches.tsv: TSV file containing all BLAST matches.
            blast_by_cluster: Folder containing BLAST results by cluster.
                cluster_x.tsv: BLAST results for cluster x.
                cluster_y.tsv: BLAST results for cluster y.
                ...: Other cluster results.
            blast_results_by_class: Folder containing BLAST results by class.
                class_1.tsv: BLAST results for class 1.
                class_2.tsv: BLAST results for class 2.
                ...: Other class results.
        cds_id_changes.tsv: TSV file containing changes in CDS IDs.
        dropped_cds.tsv: TSV file containing dropped CDS.
        Graph_folder: Folder containing graphs.
            All_of_CDS_graphs.html: HTML file containing all CDS graphs.
            graphs_class_1a.html: HTML file containing class 1a graphs.
            ...: Other graph files.

    **count_results_by_cluster.tsv**: TSV file containing count results by cluster.
    **drop_loci_reason.tsv**: TSV file containing reasons for dropping loci.
    **recommendations.tsv**: TSV file containing recommendations.
    **related_matches.tsv**: TSV file containing related matches.
    **temp_fastas**: Folder containing temporary FASTA files.
        **cluster_x.fasta**: Temporary FASTA file for cluster x.
        **cluster_y.fasta**: Temporary FASTA file for cluster y.
        **...**: Other temporary FASTA files.
    **temp_fastas_path.txt**: Text file containing paths to temporary FASTA files.

Examples
--------

Here are some example commands to use the `IdentifySpuriousGenes` module:

.. code-block:: bash

    # Identify spurious genes using default parameters
    SR IdentifySpuriousGenes -s /path/to/schema -o /path/to/output -a /path/to/allelecall

    # Identify spurious genes with custom parameters
    SR IdentifySpuriousGenes -s /path/to/schema -o /path/to/output -a /path/to/allelecall -pnl /path/to/possible_new_loci -at 0.9 -pt 90 -cs 0.9 -cc 0.9 -gp 10 -as 201 -tt 11 -b 0.6 -sr 0.8 -m schema -pm reps_vs_alleles -c 4 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `IdentifySpuriousGenes` module, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and allele call directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
