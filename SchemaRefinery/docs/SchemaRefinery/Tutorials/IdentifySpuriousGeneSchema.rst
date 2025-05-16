IdentifySpuriousGeneSchema - Tutorial
=====================================

Objective
---------

This tutorial will guide you through the process of using the `IdentifySpuriousGene` module for `Schema` to identify spurious genes in a genomic schema.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)
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
    
    Loci_id	Action
    GCA-000730255-protein547	Join
    GCA-000427055-protein583	Join
    GCA-000730255-protein547	Choice
    GCA-000730215-protein2131	Choice
    #
    GCA-001448985-protein1872	Join
    GCA-000007265-protein1932	Join
    GCA-000730215-protein1962	Join
    #


Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.

Conclusion
----------

You have successfully identified spurious genes in a genomic schema using the `IdentifySpuriousGenes` module for `Schema`.

For more information on the `IdentifySpuriousGenes` module, refer to the `IdentifySpuriousGenes documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/IdentifySpuriousGenes.html>`_.