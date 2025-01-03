DownloadAssemblies - Download assemblies and their metadata from specified databases
====================================================================================

Description
-----------

The `DownloadAssemblies` module is a module designed to facilitate the download of genomic assemblies from specified databases. This module parses command-line arguments to initiate the download process, allowing users to efficiently retrieve assemblies from either the NCBI or ENA661 databases. The downloaded assemblies are stored in the specified output directory.

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

::

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

Algorithm Explanation
---------------------

The `DownloadAssemblies` Workflow is shown in the flowchart below:

Workflow for downloading assemblies from NCBI:

.. image:: source/DownloadAssemblies_ncbi.png
   :alt: SchemaAnnotation Flowchart
   :width: 80%
   :align: center

Workflow for downloading assemblies from ENA661K:

.. image:: source/DownloadAssemblies_ena661k.png
   :alt: SchemaAnnotation Flowchart
   :width: 80%
   :align: center

Workflow for downloading metadata:

.. image:: source/DownloadAssemblies_fetch_metadata.png
   :alt: SchemaAnnotation Flowchart
   :width: 80%
   :align: center

Filtering criteria example
--------------------------
Filtering criteria file should be a TSV file with the following columns:

.. code-block:: tsv

    abundance   0.8
    genome_size 2000000
    size_threshold  0.2
    max_contig_number   150
    known_st    False
    any_quality False
    ST_list_path    None
    assembly_level  chromosome,complete,contig,scaffold
    reference   False
    assembly_source all
    file_to_include genome,gbff
    verify_status   True
    exclude_atypical    True

Note: The filtering criteria file is only applicable to certain databases e.g ST_list_path to ENA661K since it is known at the ENA661K table.

Outputs
-------
Folder and file structure for the output directory of the `DownloadAssemblies` module is shown below. The output directory contains the following files and folders:

::

    OutputFolderName
    ├── assemblies_ncbi.zip # -db NCBI --download
    ├── ena661k_assemblies # -db ENA661 --download
    │   ├── x.contigs.fa.gz
    │   ├── y.contigs.fa.gz
    │   ├── z.contigs.fa.gz
    │   └── ...
    ├── metadata_all # -fm
    │   ├── biosamples_ids.tsv
    │   ├── id_matches.tsv
    │   ├── all_ids_fetched.tsv
    │   └── metadata_biosamples.tsv
    ├── assemblies_metadata_ena661k.tsv # -db ENA661k
    ├── assemblies_metadata_ncbi.tsv # -db NCBI
    ├── metadata_ncbi # -db NCBI --nocleanup
    │   ├── assemblies_ids_to_download.tsv
    │   └── ids_failed_criteria.tsv
    └── metadata_ena661k # -db ENA661k --nocleanup
        ├── assemblies_ids_to_download.tsv
        ├── failed_to_download.tsv
        └── id_failed_criteria.tsv

.. toctree::
   :maxdepth: 1

   DownloadAssembliesOutputExplanation

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