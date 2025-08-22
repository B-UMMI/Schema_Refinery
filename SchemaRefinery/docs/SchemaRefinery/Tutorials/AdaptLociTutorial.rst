AdaptLoci - Tutorial
====================

Objective
---------

This tutorial will guide you through the process of using the :doc:`AdaptLoci </SchemaRefinery/Modules/AdaptLoci>` module to adapt loci from FASTA files into a chewBBACA-compatible schema.

Prerequisites
-------------

- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open a terminal window.
2. Modify the following command and run it to adapt loci from FASTA files into a chewBBACA-compatible schema:

.. code-block:: bash

    SR AdaptLoci -i /path/to/files/tutorial_schema/schema_seed  -o path/to/files/AdaptLoci_Results -tt 11 -c 6

.. important::

	Replace `path/to/files/` with the actual path to the files.

3. Check the output folder containing the adapted schema. The structure of the output folder should match the structure of the input schema, which was already compatible with chewBBACA.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`AdaptLoci documentation </SchemaRefinery/Modules/AdaptLoci>`.
