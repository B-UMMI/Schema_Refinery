AdaptLoci - Adapt fastas into a chewBBACA schema
================================================

Description
-----------
The `AdaptLoci` module is designed to adapt loci from FASTA files into a chewBBACA-compatible schema. This module processes input loci, evaluates allele similarity, and selects representative alleles based on configurable parameters. It supports parallel processing to enhance performance and efficiency, making it suitable for large genomic datasets. The module ensures that the adapted loci are ready for downstream analysis and schema refinement tasks.

This module is essential for researchers and bioinformaticians working on genomic schema refinement, providing a robust and flexible tool for adapting loci into a standardized schema format.

Features
--------

- Adapting loci in a schema.
- Configurable parameters for the adaptation process.
- Support for parallel processing using multiple CPUs.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)

Usage
-----

The `AdaptLoci` module can be used as follows:

.. code-block:: bash

    SR AdaptLoci -i /path/to/input_file.tsv -o /path/to/output -c 4 -b 0.6 -tt 11

Command-Line Arguments
----------------------

::

    -i, --input-file
        (Required) TSV file with the loci path to be adapted.

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
