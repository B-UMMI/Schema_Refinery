MatchSchemas - Tutorial
=======================

Objective
---------

This tutorial will guide you through the process of using the `MatchSchemas` module to match genomic schemas.

Prerequisites
-------------

- SchemaRefinery installed
- Python between 3.9 and 3.11
- Download the schema file from the `ChewBACCA tutorial https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip`_

Procedure
---------

1. Open a terminal window.

2. Run the following command to match two genome schemas:

.. code-block:: bash

    SR MatchSchemas -fs /path/to/tutorial_schema -ss /path/to/tutorial_schema -o path/to/output_folder -c 6 --nocleanup

- Replace `path/to/files/` with the actual path to the files.
- This will run the tutorial schema against itself.
- No clean up will leave all the intermediary files in the output folder.

3. Press Enter to execute the command.

4. Wait for the matching process to complete.

5. Check the output folder for the matched schemas (The output file should have all matches).

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `MatchSchemas documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/MatchSchemas.html>`_.

Conclusion
----------

You have successfully matched genomic schemas using the `MatchSchemas` module.

For more information on the `MatchSchemas` module, refer to the `MatchSchemas documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/MatchSchemas.html>`_.