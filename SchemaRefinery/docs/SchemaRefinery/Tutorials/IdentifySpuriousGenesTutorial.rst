IdentifySpuriousGenes - Tutorial
===============================

Objective
---------

This tutorial will guide you through the process of using the `IdentifySpuriousGenes` module to identify spurious genes in genomic data.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to identify spurious genes in genomic data:

.. code-block:: bash

    SR IdentifySpuriousGenes -i /path/to/input_folder -o /path/to/output_folder -a /path/to/allele_call_folder -t 4

- Replace `/path/to/input_folder` with the path to the input folder containing the genomic data.
- Replace `/path/to/output_folder` with the path to the output folder.
- Replace `/path/to/allele_call_folder` with the path to the allele call folder.
