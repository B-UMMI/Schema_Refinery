DownloadAssemblies - Download assembliesand their metadata from specified databases
====================================================================================

Description
-----------

The `DownloadAssemblies` module parses command-line arguments and initiates the downloads assemblies
from the specified databases and stores them in the output directory. It downloads the assemblies in
parallel using the number of threads specified by the user. The database options are NCBI or ENA661.
Assemblies can be download based on the taxon name or IDs table provided by the user. The user can also provide a
filtering criteria file to filter the assemblies before downloading. The user can also download BioSample metadata
for the assemblies if the flag is provided.

Overview
--------

The `DownloadAssemblies` module is designed to facilitate the download of genomic assemblies from specified databases.
It supports parallel downloads, filtering based on user-defined criteria, and the retrieval of associated metadata.

Features
--------

- Parallel downloading of assemblies using multiple threads.
- Support for downloading from NCBI and ENA661 databases.
- Filtering of assemblies based on criteria such as genome size, contig number, and assembly level.
- Retrieval of BioSample metadata for downloaded assemblies.
- Option to download assemblies based on taxon name or a provided IDs table.

Dependencies
------------

- Python 3.9 or higher
- Requests library (`pip install requests`)
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Usage
-----

The `DownloadAssemblies` module can be used as follows:

.. code-block:: bash

    SR DownloadAssemblies -t "Streptococcus pyogenes" -db NCBI ENA661K -o /path/to/output -e email@example -th 4 -fm --download

Command-Line Arguments
----------------------

-db, --database
    (Required) Databases from which assemblies will be downloaded.
    Choices: NCBI, ENA661

-o, --output-directory
    (Required) Path to the output directory.

-e, --email
    (Required) Email provided to Entrez.

-t, --taxon
    (Optional) Scientific name of the taxon. Note: This option works only for genus and species for ENA661K while for NCBI can be any taxon.
    Type: str

-th, --threads
    (Optional) Number of threads used for download. You should provide an API key to perform more requests through Entrez.
    Default: 1

-r, --retry
    (Optional) Maximum number of retries when a download or request fails.
    Default: 7

-k, --api-key
    (Optional) Personal API key provided to the NCBI. If not set, only 3 requests per second are allowed through Entrez. With a valid API key the limit increases to 10 requests per second.

-fm, --fetch-metadata
    (Optional) If provided, the process downloads metadata for the assemblies.
    Default: False

-f, --filtering-criteria
    (Optional) TSV file containing filtering parameters applied before assembly download.

--download
    (Optional) If the assemblies that passed the filtering criteria should be downloaded.

-i, --input-table
    (Optional, specific for NCBI) Text file with a list of accession numbers for the NCBI Assembly database.

Filtering criteria example
--------------------------
Filtering criteria file should be a TSV file with the following columns:
.. code-block:: tsv

    abundance\t0.8
    genome_size\t2000000
    size_threshold\t0.2
    max_contig_number\t150
    known_st\tFalse
    any_quality\tFalse
    ST_list_path\tNone
    assembly_level\tchromosome,complete,contig,scaffold
    reference\tFalse
    assembly_source\tall
    file_to_include\tgenome,gbff
    verify_status\tTrue
    exclude_atypical\tTrue

Note: The filtering criteria file is only applicable to certain databases e.g ST_list_path to ENA661K since it is known at the ENA661K table.

Outputs
-------
Folder and file structure for the output directory of the `DownloadAssemblies` module is shown below. The output directory contains the following files and folders:
.. code-block:: bash
    OutputFolderName
    ├── assemblies_ncbi.zip # -db NCBI --download
    ├── ena661k_assemblies # -db ENA661 --download
        ├── x.contigs.fa.gz
        ├── y.contigs.fa.gz
        |── z.contigs.fa.gz
        └── ...
    ├── metadata_all # -fm
        ├── biosamples_ids.tsv
        ├── id_matches.tsv
        ├── all_ids_fetched.tsv
        └── metadata_biosamples.tsv
    |── assemblies_metadata_ena661k.tsv # -db ENA661k
    |── assemblies_metadata_ncbi.tsv # -db NCBI
    ├── metadata_ncbi # -db NCBI --nocleanup
        |── assemblies_ids_to_download.tsv
        └── id_failed_criteria.tsv
    └── metadata_ena661k # -db ENA661k --nocleanup
        |── assemblies_ids_to_download.tsv
        |── failed_to_download.tsv
        └── id_failed_criteria.tsv

Output files and folders description:
-------------------------------------

**assemblies_ncbi.zip**
    Zip file containing all the assemblies and extra information that user wants downloaded from NCBI.

**ena661k_assemblies:** Folder containing the assemblies downloaded from ENA661K.
    **x.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly.
    **y.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly.
    **z.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly.
    **...**

**metadata_all:** Folder containing all the metadata downloaded from NCBI and ENA661K.
    **biosamples_ids.tsv**
        TSV file containing the BioSample IDs for the assemblies.
    **id_matches.tsv**
        TSV file containing the matches between the BioSample IDs and the assembly IDs and SRA IDs.
    **all_ids_fetched.tsv**
        TSV file containing all the IDs fetched from the database.
    **metadata_biosamples.tsv**
        TSV file containing the metadata for the BioSamples.

**assemblies_metadata_ena661k.tsv**
    TSV file containing the selected samples from the ENA661K database.

**assemblies_metadata_ncbi.tsv**
    TSV file containing the metadata for the assemblies downloaded from NCBI.

**metadata_ncbi:** Folder containing metadata related NCBI run.
    **assemblies_ids_to_download.tsv**
        TSV file containing the assembly IDs to download.
    **id_failed_criteria.tsv**
        TSV file containing the assembly IDs that failed the filtering criteria.

**metadata_ena661k:** Folder containing metadata related to ENA661K run.
    **assemblies_ids_to_download.tsv**
        TSV file containing the assembly IDs to download.
    **failed_to_download.tsv**
        TSV file containing the assembly IDs that failed to download.
    **id_failed_criteria.tsv**
        TSV file containing the assembly IDs that failed the filtering criteria.
    
Examples
--------

Here are some example commands to use the `DownloadAssemblies` module:

.. code-block:: bash

    # Download assemblies from NCBI for a specific taxon
    SR DownloadAssemblies -t "Escherichia coli" -db NCBI -o /path/to/output -e email@example.com -th 4 --download

    # Download assemblies from ENA661K using an IDs table
    SR DownloadAssemblies -db ENA661K -o /path/to/output -e email@example.com -th 4 --download -i ids_table.tsv

    # Download assemblies from both NCBI and ENA661K with filtering criteria
    SR DownloadAssemblies -t "Streptococcus pyogenes" -db NCBI ENA661K -o /path/to/output -e email@example.com -th 4 -fm --download

Troubleshooting
---------------

If you encounter issues while using the `DownloadAssemblies` module, consider the following troubleshooting steps:

- Ensure that you have a stable internet connection.
- Verify that your email and API key (if provided) are correct.
- Check the output directory for any error logs or messages.
- Increase the number of retries using the `-r` or `--retry` option if downloads are failing.