SchemaAnnotation - Tutorial
===========================

Objective
---------

This tutorial will guide you through the process of using the `SchemaAnnotation` module to annotate genomic schemas.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)
- Download the zenodo file from the link: https://zenodo.org/record/5560007/files/genbank_files.zip

Procedure
---------

1. Open a terminal window.

2. Run the following command to annotate genomic schemas:

.. code-block:: bash

    SR SchemaAnnotation -s 'path/to/files/zenodo/Data/mpneumoniae_schema/mpneumoniae_schema' -o 'path/to/files/output_folder/SchemaAnnotation_Results' -ao genbank -gf 'path/to/files/zenodo/Data/genbanks' -c 6 -tt 4

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the annotation process to complete.

5. Check the output folder for the annotated genomic schemas (The results are also available in the zenodo files).

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`SchemaAnnotation documentation <SchemaRefinery/Modules/Schema_annotation>`.

Conclusion
----------

You have successfully annotated genomic schemas using the `SchemaAnnotation` module.

For more information on the `SchemaAnnotation` module, refer to the :doc:`SchemaAnnotation documentation <SchemaRefinery/Modules/Schema_annotation>`.