AdaptLoci - Adapt fastas into a chewBBACA schema
================================================

Description
-----------
The `AdaptLoci` module is designed to adapt a set of loci in FASTA format into a chewBBACA-compatible schema. This is the schema structure that is used by Schema Refinery. This module uses chewBBACA's PrepExternalSchema module to validate and create schemas. It also checks the format of the loci names to find names that can lead to issues when running BLAST.


Features
--------

- Adapting sets of loci into a chewBBACA-compatible schemas.
- Checks and warns users about loci names that can lead to issues when running BLAST processes.
- Configurable parameters for the adaptation process.
- Support for parallel processing using multiple CPUs.

Dependencies
------------

- chewBBACA 3.3.10 or higher (`chewBBACA's installation instructions <https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html>`_).

Usage
-----

The `AdaptLoci` module can be used as follows:

.. code-block:: bash

    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output -c 4 -b 0.6 -tt 11

Command-Line Arguments
----------------------

::

    -i, --input-file
        (Required) Path to the folder with the fasta files.

    -o, --output-directory
        (Required) Path to the directory to which files will be stored.
    
    -tf, --training-file
        (Optional) Path to the Prodigal training file that will be included in the directory of the adapted schema.

    -c, --cpu
        (Optional) Number of CPUs to run BLAST instances.
        Default: 1

    -b, --bsr
        (Optional) BSR value to consider alleles as the same locus.
        Default: 0.6

    -tt, --translation_table
        (Optional) Translation table to use for the CDS translation.
        Default: 11

    --debug
        (Optional) Flag to indicate whether to run the module in debug mode.
        Default: False

    --logger
        (Optional) Path to the logger file.
        Default: None

.. Important::
    Always verify it the translation table value passed to the `-tt` parameter is the correct one for the species.


Algorithm Explanation
---------------------

The `AdaptLoci` module adapts a set of loci in FASTA format into a chewBBACA-compatible schema, which includes a folder names `short` with FASTA files containing the representative alleles of each locus. To adapt the loci files, the module calls chewBBACA's PrepExternalSchema module with a simplified list of arguments.
For further details about the adaptation process, visit the `PrepExternalSchema documentation <https://chewbbaca.readthedocs.io/en/latest/user/modules/PrepExternalSchema.html#>`_.

The input should be a folder with a fasta file per loci.


Outputs
-------
The structure of the output directory is shown below.

::

    OutputFolderName
    ├── x.fasta
    ├── y.fasta
    ├── z.fasta
    ├── ...
    └── short
        ├── x_short.fasta
        ├── y_short.fasta
        ├── z_short.fasta
        └── ...

.. toctree::
   :maxdepth: 1

   AdaptLociOutputDescription

Examples
--------

Simple usage commands to run the `AdaptLoci` module:

.. code-block:: bash

    # Adapt loci using the default parameters
    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output

    # Adapt loci using custom parameters
    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output -tf /path/to/file/training-file.trn -c 4 -b 0.7 -tt 4

Troubleshooting
---------------

If you encounter issues while using the `AdaptLoci` module, consider the following troubleshooting tips:

- Verify that the paths to the input files and output directory are valid.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
