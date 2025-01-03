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
- Download the zenodo file from the link: https://zenodo.org/record/5560007/files/genbank_files.zip

Procedure
---------

1. Open a terminal window.

2. Run the following command to cg/wgMSLT genomic schemas:

.. code-block:: bash

    SR MatchSchema -qs 'path/to/files/zenodo/Data/mpneumoniae_schema/mpneumoniae_schema' -ss 'path/to/files/zenodo/Data/mpneumoniae_schema/mpneumoniae_schema' -o 'path/to/files/output_folder/MatchSchemas_Results' -tt 4 -pm alleles_vs_alleles -c 6

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the matching process to complete.

5. Check the output folder for the matched schemas (The results are also available in the zenodo files).

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `MatchSchemas documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/MatchSchemas.html>`_.

Conclusion
----------

You have successfully matched genomic schemas using the `MatchSchemas` module.

For more information on the `MatchSchemas` module, refer to the `MatchSchemas documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/MatchSchemas.html>`_.