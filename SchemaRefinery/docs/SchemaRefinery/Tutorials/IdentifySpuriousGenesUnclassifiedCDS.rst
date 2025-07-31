IdentifySpuriousGenesUnclassifiedCDS - Tutorial
===============================================

Objective
---------

This tutorial will guide you through the process of using the `IdentifySpuriousGenes` module for `Unclassified CDS` to classify unclassified coding sequences (CDS).

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

    SR IdentifySpuriousGenes -s '/path/to/tutorial_schema/schema_seed' -a '/path/to/Allele_calling' -o '/path/to/files/output_folder/IdentifySpuriousGenesUnclassifiedCDS' -m unclassified_cds --t 11 -c 6

- Replace `/path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the identified spurious genes.
    The first lines of the final clusters file should look like:
::
    
    Locus	Action	Class
    GCA-000831145-protein1681	Join	1a
    GCA-000831105-protein622	Join	1a
    GCA-000012705-protein568	Join	1a
    GCA-001275545-protein1163	Join	1a
    GCA-000007265-protein582	Join	1a
    GCA-000730215-protein582	Join	1a
    GCA-000730255-protein607	Join	1a
    GCA-000427075-protein661	Choice	3b
    #
    GCA-000427055-protein712	Join	1a
    GCA-000730255-protein664	Join	1a
    #


Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.

Conclusion
----------

You have successfully classified unclassified coding sequences (CDS) using the `IdentifySpuriousGenes` module for `Unclassified CDS`.

For more information on the `IdentifySpuriousGenes` module, refer to the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.
