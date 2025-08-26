IdentifySpuriousGeneSchema - Tutorial
=====================================

Objective
---------

This tutorial will guide you through the process of using the :doc:`IdentifySpuriousGene </SchemaRefinery/Modules/IdentifySpuriousGenes>` module to identify spurious genes in a schema.

Prerequisites
-------------

- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open the terminal

2. Modify and run the following command to identify spurious genes in the schema:

.. code-block:: bash

    SR IdentifySpuriousGenes -s '/path/to/tutorial_schema/schema_seed' -a '/path/to/Allele_calling'  -o '/path/to/files/output_folder/IdentifySpuriousGenesSchema' -m schema -pm alleles_vs_alleles --t 11 -c 6

.. important::
	Replace `/path/to/files/` with the actual path to the files.

1. Check the output directory for the list of identified spurious genes. The first lines of the file containing the clusters of potential spurious loci should look like this:

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

To see the expected output structure, refer to the "Outputs" section in the :doc:`IdentifySpuriousGenes documentation </SchemaRefinery/Modules/IdentifySpuriousGenesOutputDescription>`.

Conclusion
----------

You have successfully identified spurious genes in a schema using the :doc:`IdentifySpuriousGene </SchemaRefinery/Modules/IdentifySpuriousGenes>`.
