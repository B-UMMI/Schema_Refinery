IdentifySpuriousGeneSchema - Tutorial
=====================================

Objective
---------

This tutorial will guide you through the process of using the `IdentifySpuriousGene` module for `Schema` to identify spurious genes in a genomic schema.

Prerequisites
-------------

- SchemaRefinery installed
- Python between 3.9 and 3.11
- Download the schema file from the `ChewBACCA tutorial https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip`_

Procedure
---------

1. Open the terminal

2. Run the following command to execute

.. code-block:: bash

    SR IdentifySpuriousGenes -s '/path/to/tutorial_schema/schema_seed' -a '/path/to/Allele_calling'  -o '/path/to/files/output_folder/IdentifySpuriousGenesSchema' -m schema -pm alleles_vs_alleles --t 11 -c 6

- Replace `/path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the identified spurious genes.
    The first lines of the final clusters file should look like:
::
    
   Locus	Action	Class
    GCA-000730255-protein547	Join	1a
    GCA-000427055-protein583	Join	1a
    GCA-000730215-protein2131	Choice	4b
    GCA-000196055-protein1223	Choice	1c
    GCA-000007265-protein1233	Choice	1c
    GCA-000007265-protein534	Choice	1c
    GCA-000012705-protein1877	Choice	4b
    #
    GCA-000730215-protein1962	Join	1a
    GCA-000007265-protein1932	Join	1a
    GCA-000196055-protein485	Choice	1c
    GCA-000196055-protein1146	Choice	1c
    GCA-000196055-protein398	Choice	1c
    GCA-000427075-protein1286	Choice	1c
    #


Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.

Conclusion
----------

You have successfully identified spurious genes in a genomic schema using the `IdentifySpuriousGenes` module for `Schema`.

For more information on the `IdentifySpuriousGenes` module, refer to the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.