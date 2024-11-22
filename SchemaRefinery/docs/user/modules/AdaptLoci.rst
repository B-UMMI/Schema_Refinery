AdaptLoci - Adapt loci in a schema
==================================

Description
-----------

The `AdaptLoci` module parses command-line arguments and initiates the process to adapt loci. This module sets up an argument parser to handle various command-line options for adapting loci and then calls the main function of the `AdaptLoci` class with the parsed arguments.

The `AdaptLoci` takes as input a set of fastas and returns a chewBBACA schema.

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

    python schema_refinery.py -i /path/to/input_file.tsv -o /path/to/output -c 4 -b 0.6 -tt 11

Command-Line Arguments
----------------------

-i, --input_file
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

.. code-block:: bash

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

Output files and folders description:
-------------------------------------
**OutputFolderName**: The folder where the output files are stored.
    **x.fasta**: The fasta file containing the loci.
    **y.fasta**: The fasta file containing the loci.
    **z.fasta**: The fasta file containing the loci.
    **...**: Other fasta files containing the loci.
    **short**: The folder containing the short loci.
        **x_short.fasta**: The short fasta file containing the loci.
        **y_short.fasta**: The short fasta file containing the loci.
        **z_short.fasta**: The short fasta file containing the loci.
        **...**: Other short fasta files containing the loci.


Examples
--------

Here are some example commands to use the `AdaptLoci` module:

.. code-block:: bash

    # Adapt loci using default parameters
    python schema_refinery.py -i /path/to/input_file.tsv -o /path/to/output

    # Adapt loci with custom parameters
    python schema_refinery.py -i /path/to/input_file.tsv -o /path/to/output -c 4 -b 0.7 -tt 1

Troubleshooting
---------------

If you encounter issues while using the `AdaptLoci` module, consider the following troubleshooting steps:

- Verify that the paths to the input file and output directory are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
