AdaptLoci - Adapt fastas into a chewBBACA schema
================================================

Description
-----------
The `AdaptLoci` module is designed to adapt loci from FASTA files into a chewBBACA-compatible schema. This is teh schema structure that is assumed all throughout the tool. This module processes input loci and runs the module PrepExternalSchema from chewBBACA.

This module is essential for researchers and bioinformaticians working on genomic schema refinement, providing a robust and flexible tool for adapting loci into a standardized schema format.

Features
--------

- Adapting loci in a schema.
- Configurable parameters for the adaptation process.
- Support for parallel processing using multiple CPUs.

Dependencies
------------

- Python 3.9 or higher
- BLAST (`https://www.ncbi.nlm.nih.gov/books/NBK279690/ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- ChewBBACA (`https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html <https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html>`_ or using bioconda)
- Install requirements using the following command:

.. code-block:: bash

    pip install -r requirements.txt

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


Algorithm Explanation
---------------------
This module calls on the PrepExternalSchema module from chewBBACA with a simplified list of arguments.
For more clarification check the PrepExternalSchema documentation (`https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html <https://chewbbaca.readthedocs.io/en/latest/user/modules/PrepExternalSchema.html#>`_)


Outputs
-------
Folder and file structure for the output directory of the `AdaptLoci` module is shown below. The output directory contains the following files and folders:

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

Here are some example commands to use the `AdaptLoci` module:

.. code-block:: bash

    # Adapt loci using default parameters
    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output

    # Adapt loci with custom parameters
    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output -c 4 -b 0.7 -tt 4

Troubleshooting
---------------

If you encounter issues while using the `AdaptLoci` module, consider the following troubleshooting steps:

- Verify that the paths to the input file and output directory are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
