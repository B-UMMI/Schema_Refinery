SchemaAnnotation - Tutorial
===========================

Objective
---------

This tutorial will guide you through the process of using the :doc:`SchemaAnnotation </SchemaRefinery/Modules/SchemaAnnotation>` module to annotate a schema.

Prerequisites
-------------

- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.
- Download the list of proteomes for *Streptococcus agalactiae* from `Uniprot <https://www.uniprot.org/proteomes?query=Streptococcus+agalactiae>`_.

Procedure
---------

1. Open a terminal window.

2. Modify and run the following command to annotate the tutorial schema:

.. code-block:: bash

    SR SchemaAnnotation -s 'path/to/tutorial_schema/schema_seed' -o 'path/to/output_folder' -ao uniprot-proteomes -pt 'path/to/unzipped/proteome_file' -c 6 -tt 11 --nocleanup

.. important::
	Replace `path/to/files/` with the actual path to the files.

3. Check the output folder for the loci annotations (it should have all columns of both files except for the Proteome ID column of the annotation file).

.. Note::
	For the other annotation options the module follows the same procedure but uses the appropriate arguments.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`SchemaAnnotation documentation </SchemaRefinery/Modules/SchemaAnnotation>`.

Conclusion
----------

You have successfully annotated a schema using the `SchemaAnnotation` module.
