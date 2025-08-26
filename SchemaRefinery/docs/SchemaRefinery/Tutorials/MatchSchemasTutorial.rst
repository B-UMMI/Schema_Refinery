MatchSchemas - Tutorial
=======================

Objective
---------

This tutorial will guide you through the process of using the :doc:`MatchSchemas </SchemaRefinery/Modules/MatchSchemas>` module to find matches between the loci in two schemas.

Prerequisites
-------------

- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open a terminal window.

2. Modify and run the following command to compare the loci in two schemas:

.. code-block:: bash

    SR MatchSchemas -fs /path/to/tutorial_schema -ss /path/to/tutorial_schema -o path/to/output_folder -c 6 --nocleanup

.. important::
	Replace `path/to/files/` with the actual path to the files. This will run the tutorial schema against itself. The `--nocleanup` parameter will leave all the intermediary files in the output folder.

1. Check the outptut folder containing the list of loci in both schemas that are similar.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`MatchSchemas documentation </SchemaRefinery/Modules/MatchSchemas>`.

Conclusion
----------

You have successfully matched two schemas using the :doc:`MatchSchemas </SchemaRefinery/Modules/MatchSchemas>` module.
