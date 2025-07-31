SchemaAnnotation - Tutorial
===========================

Objective
---------

This tutorial will guide you through the process of using the `SchemaAnnotation` module to annotate genomic schemas.

Prerequisites
-------------

- SchemaRefinery installed
- Python between 3.9 and 3.11
- Download the schema file from the `ChewBACCA tutorial https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip`_
- Download the proteome files from `Uniprot https://www.uniprot.org/proteomes?query=Streptococcus+agalactiae`_

Procedure
---------

1. Open a terminal window.

2. Run the following command to annotate genomic schemas:

.. code-block:: bash

    SR SchemaAnnotation -s 'path/to/tutorial_schema/schema_seed' -o 'path/to/output_folder' -ao uniprot-proteomes -pt 'path/to/unzipped/proteome_file' -c 6 -tt 11 --nocleanup

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the annotation process to complete.

5. Check the output folder for the annotated genomic schemas (it should have all columns of both files except for the Proteome ID column of the annotation file).

.. Note:: For the other options the module follow the same procedure just using the appropriate arguments.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `SchemaAnnotation documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/SchemaAnnotation.html>`_.

Conclusion
----------

You have successfully annotated genomic schemas using the `SchemaAnnotation` module.

For more information on the `SchemaAnnotation` module, refer to the `SchemaAnnotation documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/SchemaAnnotation.html>`_.