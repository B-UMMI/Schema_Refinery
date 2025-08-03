CreateSchemaStructure - Create a new schema based on recommendations
=====================================================================

Description
------------
The `CreateSchemaStructure` module is a module designed to facilitate the creation of a schema structure from a given schema or fasta files based on user recommendation file from the `IdentifyParalogousLoci` and `IdentifySpuriousGenes` modules. This module parses command-line arguments to initiate the schema creation process, allowing users to efficiently generate a schema structure from a provided schema file. The generated schema structure is stored in the specified output directory.

This is the final step in the `SchemaRefinery` workflow providing as output a schema reflecting the user reviewed changes that were recommended by other modules.


Features
--------
- Creation of a schema structure from a given schema or fasta files.
- Takes recommendation files, annotated or not, as guidelines for the new schema.
- Support for parallel processing using multiple CPUs.
- Option to skip cleanup after running the module.

Dependencies
------------
- Python 3.9 or higher
- BLAST (`https://www.ncbi.nlm.nih.gov/books/NBK279690/ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- ChewBBACA (https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html or using bioconda)
- Install requirements using the following command:

.. code-block:: bash

    pip install -r requirements.txt

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
    Always verify it the translation table (argument -tt) being used is the correct one for the species.

The input recommendation file must have the columns "Locus" and "Action". Additionally, if from the `IdentifySpuriousGenes` module, the file can also have the column "Class". These columns must be the 1st, 2nd and 3rd columns in the file, respectively. This input can have more columns, such as annotations, but these will be ignored. This file needs to have all the loci that are meant to be in the final schema. If a locus is missing from this file it won't be included in the final version of the schema


Algorithm Explanation
---------------------

The `CreateSchemaStructure` module uses the following algorithm to create a schema structure from a given schema or fasta files:

.. image:: source/CreateSchemaStructure.png
   :alt: Algorithm for creating and restructure a schema
   :width: 80%
   :align: center


.. Note::
    Make sure that before running this module the input file `recommendations` has been check and that you agree with all the changes proposed. **No action should have the option 'Choice'.** These should be changed into one of the other three actions.

The action **'Join'** will move all the alleles of that cluster of loci into the fasta file of the first locus of the cluster. The other Loci files will be removed from the final schema.

The action **'Drop'** will remove that locus fasta file from the final schema.

The action **'Add'** will just copy the fasta file of that locus from the input schema into the final schema with no alterations.

If the locus is within the class 6, the algorithm will check what the action is. If the action is "Join" the run will end and the user will be asked to change it, as this pairing is not allowed. The same will happen if an action "Choice" is found.

Since this module uses the AdaptLoci module to format the schema in the end, the schema created will conform with the structure of the schemas used by chewBBACA.

Outputs
-------
Folder and file structure for the output directory of the `CreateSchemaStructure` module is shown below. The output directory contains the following files and folders:

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

This file will have the main statistics of the final schema. It gives the number of alleles in each locus and how many of those are representatives.

`schema_invalid_alleles.txt` and `schema_invalid_loci.txt` will have a list of the alleles/loci that did not pass validation thresholds. If these are empty then all alleles/loci were correctly moved into the new schema and validated.


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

If you encounter issues while using the `CreateSchemaStructure` module, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and recommendations directories are correct.
- Check the output directory for any error logs or messages.
- Check if the dependencies versions are all compatible. 
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
