MatchSchemas - Tutorial
=======================

Objective
---------

This tutorial will guide you through the process of using the `MatchSchemas` module to match genomic schemas.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to cg/wgMSLT genomic schemas:

.. code-block:: bash

    SR MatchSchemas -i /path/to/input_folder -o /path/to/output_folder -t 4

- Replace `/path/to/input_folder` with the path to the folder containing the genomic schemas.
- Replace `/path/to/output_folder` with the path to the output folder.

3. Press Enter to execute the command.

4. Wait for the matching process to complete.

5. Check the output folder for the matched schemas.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`MatchSchemas documentation <SchemaRefinery/Modules/Match_schemas>`.

Conclusion
----------

You have successfully matched genomic schemas using the `MatchSchemas` module.

For more information on the `MatchSchemas` module, refer to the :doc:`MatchSchemas documentation <SchemaRefinery/Modules/Match_schemas>`.