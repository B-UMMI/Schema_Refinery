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

Procedure
---------

1. Open a terminal window.

2. Run the following command to annotate genomic schemas:

.. code-block:: bash

    SR SchemaAnnotation -i /path/to/input_folder -o /path/to/output_folder  -ao genbank-files -gf path/to/genbank/files

- Replace `/path/to/input_folder` with the path to the input folder containing the genomic schemas.
- Replace `/path/to/output_folder` with the path to the output folder.
- Replace `path/to/genbank/files` with the path to the folder containing the genbank files.

3. Press Enter to execute the command.

4. Wait for the annotation process to complete.

5. Check the output folder for the annotated genomic schemas.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`SchemaAnnotation documentation <SchemaRefinery/Modules/Schema_annotation>`.

Conclusion
----------

You have successfully annotated genomic schemas using the `SchemaAnnotation` module.

For more information on the `SchemaAnnotation` module, refer to the :doc:`SchemaAnnotation documentation <SchemaRefinery/Modules/Schema_annotation>`.