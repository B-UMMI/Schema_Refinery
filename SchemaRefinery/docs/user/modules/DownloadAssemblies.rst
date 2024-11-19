DownloadAssemblies - Download assemblies from specified databases
===============================================================

Description
-----------

The `DownloadAssemblies` module parses command-line arguments and initiates the downloads assemblies
from the specified databases and stores them in the output directory. It downloads the assemblies in
parallel using the number of threads specified by the user. The database options are NCBI or ENA661.
Assemblies can be download based on the taxon name or IDs table provided by the user. The user can also provide a
filtering criteria file to filter the assemblies before downloading. The user can also download BioSample metadata
for the assemblies if the flag is provided.

Usage
-----

The `DownloadAssemblies` module can be used as follows:

.. code-block:: bash

    SR DownloadAssemblies -db NCBI ENA661K -o /path/to/output -e email@example -th 4 -fm

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
    (Optional) Scientific name of the taxon.
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

Outputs
-------

.. code-block:: bash
    OutputFolderName
    ├── assemblies.zip
    ├── ena661k_assemblies
        ├── x.contigs.fa.gz
        ├── y.contigs.fa.gz
        └── z.contigs.fa.gz
    ├── metadata_all # -fm
        |── biosamples.tsv
        |── id_matches.tsv
        |── all_ids_fetched.tsv
        └── biosample_biosamples.tsv
    |── selected_samples_ena661k.tsv
    ├── metadata_ncbi # --nocleanup
        |── assemblies_ids_to_download.tsv
        └── id_failed_criteria.tsv
    └── metadata_ena661k # --nocleanup
        |── assemblies_ids_to_download.tsv
        └── id_failed_criteria.tsv
        
