AdaptLoci - Tutorial
====================

Objective
---------

This tutorial will guide you through the process of using the `AdaptLoci` module to adapt loci from FASTA files into a chewBBACA-compatible schema.

Prerequisites
-------------

- SchemaRefinery installed
- Python between 3.9 and 3.11
- Download the schema file from the `ChewBACCA tutorial https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip`_

Procedure
---------

1. Open a terminal window.

2. Run the following command to adapt loci from FASTA files into a chewBBACA-compatible schema:

.. code-block:: bash

    SR AdaptLoci -i '/path/to/files/tutorial_schema/schema_seed'  -o 'path/to/files/AdaptLoci_Results' -tt 11 -c 6

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the adaptation process to complete.

5. Check the output folder for the adapted schema, it should look like the original schema folder.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `AdaptLoci documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/AdaptLoci.html>`_.

Conclusion
----------

You have successfully adapted loci in genomic schemas using the `AdaptLoci` module.

For more information on the `AdaptLoci` module, refer to the `AdaptLoci documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/AdaptLoci.html>`_.