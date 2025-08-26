IdentifyParalogousLoci - Tutorial
=================================

Objective
---------

This tutorial will show you how to identify paralogous loci in a schema using the :doc:`IdentifyParalogousLoci </SchemaRefinery/Modules/IdentifyParalogousLoci>` module.

Prerequisites
-------------
- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open a terminal window.

2. Modify the following command and run it to identify paralogous loci in a schema:

.. code-block:: bash

    SR IdentifyParalogousLoci -s 'path/to/tutorial_schema/schema_seed' -o 'path/to/files/output_folder/IdentifyParalogousLoci_Results' -tt 11 -c 6 -pm alleles_vs_alleles

.. important::
	Replace `path/to/files/` with the actual path to the files.

1. Check the output directory for the list of identified paralogous loci. The first lines of the file containing the list of clusters of paralogous loci that were identified should look like:

::
    
    Loci_id	Action
    GCA-000007265-protein1932	Join
    GCA-000730215-protein1962	Join
    #	
    GCA-000427055-protein1391	Join
    GCA-000427035-protein1421	Join
    GCA-000782855-protein1355	Join
    #	
    GCA-000012705-protein2017	Join
    GCA-000427075-protein2050	Join
    #	

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`IdentifyParalogousLoci documentation </SchemaRefinery/Modules/IdentifyParalogousLoci>`.

Conclusion
----------

You have successfully identified paralogous loci in a schema using the :doc:`IdentifyParalogousLoci documentation </SchemaRefinery/Modules/IdentifyParalogousLoci>` module.
