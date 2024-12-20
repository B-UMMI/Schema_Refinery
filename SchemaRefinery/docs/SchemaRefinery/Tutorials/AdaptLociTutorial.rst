AdaptLoci - Tutorial
====================

Objective
---------

This tutorial will guide you through the process of using the `AdaptLoci` module to adapt loci from FASTA files into a chewBBACA-compatible schema.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to adapt loci from FASTA files into a chewBBACA-compatible schema:

.. code-block:: bash

    SR AdaptLoci -i /path/to/input_folder -o /path/to/output_folder -t 4

- Replace `/path/to/input_folder` with the path to the folder containing the FASTA files.
- Replace `/path/to/output_folder` with the path to the output folder.

3. Press Enter to execute the command.

4. Wait for the adaptation process to complete.

5. Check the output folder for the adapted loci.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :docs:`AdaptLoci documentation <SchemaRefinery/Modules/AdaptLoci>`.

Conclusion
----------

You have successfully adapted loci from FASTA files into a chewBBACA-compatible schema using the `AdaptLoci` module.

For more information on the `AdaptLoci` module, refer to the :docs:`AdaptLoci documentation <SchemaRefinery/Modules/AdaptLoci>`.