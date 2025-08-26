IdentifySpuriousGenesUnclassifiedCDS - Tutorial
===============================================

Objective
---------

This tutorial will guide you through the process of using the :doc:`IdentifySpuriousGenes </SchemaRefinery/Modules/IdentifySpuriousGenes>` module to identify potential new loci or spurious alleles.

Prerequisites
-------------

- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open the terminal

2. Modify and run the following command to evaluate a set of unclassified CDSs:

.. code-block:: bash

    SR IdentifySpuriousGenes -s '/path/to/tutorial_schema/schema_seed' -a '/path/to/Allele_calling' -o '/path/to/files/output_folder/IdentifySpuriousGenesUnclassifiedCDS' -m unclassified_cds --t 11 -c 6

.. important::
	Replace `/path/to/files/` with the actual path to the files.

3. The output directory contains files with the groups of similar alleles that were identfied. The first lines of the final clusters file should look like:

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

To see the expected output structure, refer to the "Outputs" section in the :doc:`IdentifySpuriousGenes documentation </SchemaRefinery/Modules/IdentifySpuriousGenes>`.

Conclusion
----------

You have successfully classified unclassified coding sequences (CDS) using the `IdentifySpuriousGenes` module.
