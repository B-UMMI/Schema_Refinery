CreateSchemaStructure - Tutorial
===============================================

Objective
---------

This tutorial will guide you through the process of using the `CreateSchemaStructure` module for creating a new schema given a schema and a set of recommendations.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- BLAST (`https://www.ncbi.nlm.nih.gov/books/NBK279690/ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_)
- ChewBBACA (https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html or using bioconda)
- Download the schema file from the ChewBACCA tutorial: https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip

Procedure
---------

1. Open the terminal

2. Run the following command to execute

.. code-block:: bash

    SR CreateSchemaStructure -rf '/path/to/recommendations.tsv' -ff '/path/to/tutorial_schema/schema_seed' -o '/path/to/CreateSchemaStructure_output' -c 6 --nocleanup

- Replace `/path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the created schema.
    The format of the schema should be the same as of a chewBBACA schema.
    To check if the changes are correct, compare the fastas of the loci of the original schema with the ones of the new one and confirm if they match with the actions.


Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `CreateSchemaStructure documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/CreateSchemaStructure.html>`_.

Conclusion
----------

You have successfully created a new strcutured schema using the `CreateSchemaStructure` module.

For more information on the `CreateSchemaStructure` module, refer to the `CreateSchemaStructure documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/CreateSchemaStructure.html>`_.
