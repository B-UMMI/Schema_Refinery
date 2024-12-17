IdentifyParalogousLoci - Tutorial
=================================

Objective
---------

This tutorial will show you how to identify paralogous loci in a schema using the `IdentifyParalogousLoci` module.

Prerequisites
-------------
- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to identify paralogous loci in a schema:

.. code-block:: bash

    SR IdentifyParalogousLoci -s /path/to/schema_seed -o /path/to/output -c 4 -pm alleles_vs_alleles

- Replace `/path/to/schema_seed` with the path to the schema seed folder.
- Replace `/path/to/output` with the path to the output directory.
