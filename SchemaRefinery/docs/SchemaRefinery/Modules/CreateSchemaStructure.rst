CreateSchemaStructure - Create a new schema based on recommendations
=====================================================================

Description
------------
The `CreateSchemaStructure` module facilitates the creation of a schema structure from a given schema or FASTA files based on a file with a set of recommendations from the :doc:`IdentifyParalogousLoci </SchemaRefinery/Modules/IdentifyParalogousLoci>` and :doc:`IdentifySpuriousGenes </SchemaRefinery/Modules/IdentifySpuriousGenes>` modules. 

This module is used in the final step of the `Schema Refinery` workflow to create a schema that reflects the recommendations obtained with other modules and reviewed by the users.

Features
--------
- Creation of a schema structure from a given schema or set of FASTA files.
- Accepts files with recommendations that serve as guidelines to create the new schema.
- Support for parallel processing using multiple CPUs.
- Option to skip intermediate file cleanup after running the module.

Dependencies
------------

- BLAST (manual `here <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- chewBBACA 3.3.10 or higher (`chewBBACA's installation instructions <https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html>`_).

Usage
-----
The `CreateSchemaStructure` module can be used as follows:

.. code-block:: bash

    SR CreateSchemaStructure -s /path/to/schema -o /path/to/output -rf /path/to/recommendation.tsv -c 4 --nocleanup

Command-Line Arguments
----------------------

::

    -rf, --recommendations-file
        (Required) Path to the file containing the recommendations.

    -ff, --fastas-folder
        (Required) Path to the folder containing the FASTA files (Just fastas or schema).

    -o, --output-directory
        (Required) Path to the directory where the output files will be saved.
    
    -tf, --training-file
        (Optional) Path to the Prodigal training file that will be included in the directory of the adapted schema.

    -c, --cpu
        (Optional) Number of CPU cores for multiprocessing.
        Default: 1

    -bsr, --blast-score-ratio
        (Optional) BSR value to consider alleles as the same locus.
        Default: 0.6

    -tt, --translation-table
        (Optional) Translation table to use for the CDS translation.
        Default: 11

    --nocleanup
        (Optional) Flag to indicate whether to skip cleanup after running the module.
        Default: False

    --debug
        (Optional) Flag to indicate whether to run the module in debug mode.
        Default: False

    --logger
        (Optional) Path to the logger file.
        Default: None

.. Note::
    Always verify it the translation table value passed to the `-tt` parameter is the correct one for the species.

The input files with recommendations must have the columns "Locus" and "Action". Additionally, if from the `IdentifySpuriousGenes` module, the file can also have the column "Class". These columns must be the 1st, 2nd and 3rd columns in the file, respectively. This input can have more columns, such as annotations, but these will be ignored. This file needs to have all the loci that are meant to be in the final schema. If a locus is missing from this file it won't be included in the final version of the schema


Algorithm Explanation
---------------------

The `CreateSchemaStructure` module uses the following algorithm to create a schema structure from a given schema or fasta files:

.. image:: source/CreateSchemaStructure.png
   :alt: Algorithm for creating and restructure a schema
   :width: 80%
   :align: center


.. Note::
    Make sure to review the list of recommendations before running the module. The **Choice** recommendations should be carefully reviewed and changed into one of the other options. The module will perform no action for the loci for which the action is **Choice**.

The action **Join** will merge all loci in a cluster with that action into a single FASTA file. The FASTA file will contain all the distinct alleles in the original loci FASTA files and be named based on the first locus in the cluster that was merged with the otehrs.

The action **Drop** will exclude a locus from the final schema.

The action **Add** will copy the FASTA file of a locus from the input schema into the final schema.

.. Important::
	If class 6 is assigned to a locus and the action is **Join** or **Choice**, the process will halt and ask the user to change the action, as those actions are not allowed when class 6 is assigned.

Since this module uses the :doc:`AdaptLoci </SchemaRefinery/Modules/AdaptLoci>` module to format the schema in the end, the schema created will conform with the structure of the schemas used by chewBBACA.

Outputs
-------
The structure of the output directory created by the `CreateSchemaStructure` module is shown below.

::

    OutputFolderName
    ├── schema
    │   ├── x.fasta
    │   ├── y.fasta
    │   ├── z.fasta
    │   ├── ...
    │   └── short
    │       ├── x_short.fasta
    │       ├── y_short.fasta
    │       ├── z_short.fasta
    │       └── ...
    ├── temp_fasta # --nocleanup
    │   ├── x.fasta
    │   ├── y.fasta
    │   ├── z.fasta
    │   └── ...
    ├── schema_invalid_alleles.txt # --nocleanup
    ├── schema_invalid_loci.txt # --nocleanup
    └── schema_summary_stats.tsv # --nocleanup

.. toctree::
    :maxdepth: 1
    
    CreateSchemaStructureOutputDescription

Report files description
------------------------

This file includes summary statistics for the final schema.

.. csv-table:: **schema_summary_stats.tsv**
   :header: "Gene", "Total_alleles", "Valid_alleles", "Number_representatives"
   :widths: 10, 10, 10, 10

   x, 11, 11, 1
   y, 1, 1, 1
   z, 11, 11, 1
   ...


Columns description:

::

    Gene: The locus ID of the query.
    Total_alleles: The number of alleles in the input schema.
    Valid_alleles: The number of alleles chosen to be in the final schema.
    Number_representatives: The number of alleles chosen to be the representatives of that loci.

Additionally, the process creates the `schema_invalid_alleles.txt` and `schema_invalid_loci.txt` files, which include the list of the alleles and loci that did not pass the validation thresholds, respectively. If these are empty then all alleles/loci were valid and moved into the new schema.


Examples
--------

Here are some example commands to use the `CreateSchemaStructure` module:

.. code-block:: bash

    # Create schema structure from a given schema file
    SR CreateSchemaStructure -s /path/to/schema -o /path/to/output -rf /path/to/recommendation.tsv -c 4 --nocleanup

    # Create schema structure from a folder containing FASTA files
    SR CreateSchemaStructure -ff /path/to/fastas -o /path/to/output -rf /path/to/recommendation.tsv -c 4 --nocleanup

    # Create schema structure with custom BSR value and translation table and training file
    SR CreateSchemaStructure -s /path/to/schema -o /path/to/output -rf /path/to/recommendation.tsv -tf /path/to/file/training-file.trn -c 4 --t 4 -bsr 0.8 --nocleanup

Troubleshooting
---------------

If you encounter issues while using the `CreateSchemaStructure` module, consider the following troubleshooting tips:

- Verify that the paths to the schema, output directory, and file with recommendations are valid.
- Check the output directory for any error logs or messages.
- Check if the versions fo the dependencies used by the module are all compatible.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
