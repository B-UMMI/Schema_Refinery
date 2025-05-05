IdentifyParalogousLoci - Tutorial
=================================

Objective
---------

This tutorial will show you how to identify paralogous loci in a schema using the `IdentifyParalogousLoci` module.

Prerequisites
-------------
- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)
- Download the schema file from the ChewBACCA tutorial: https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip

Procedure
---------

1. Open a terminal window.

2. Run the following command to identify paralogous loci in a schema:

.. code-block:: bash

    SR IdentifyParalogousLoci -s 'path/to/tutorial_schema/schema_seed' -o 'path/to/files/output_folder/IdentifyParalogousLoci_Results' -tt 11 -c 6 -pm alleles_vs_alleles

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the identified paralogous loci.
    The first lines of the final clusters file should look like:
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

To see the expected output structure, refer to the "Outputs" section in the `IdentifyParalogousLoci documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifyParalogousLoci.html>`_.

Conclusion
----------

You have successfully identified paralogous loci in a schema using the `IdentifyParalogousLoci` module.

For more information on the `IdentifyParalogousLoci` module, refer to the `IdentifyParalogousLoci documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifyParalogousLoci.html>`_.