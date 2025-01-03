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
- Download the zenodo file from the link: https://zenodo.org/record/5560007/files/genbank_files.zip

Procedure
---------

1. Open the terminal

2. Run the following command to execute

.. code-block:: bash

    SR IdentifySpuriousGenes -s '/path/to/files/zenodo/Data/mpneumoniae_schema/mpneumoniae_schema' -a '/path/to/files/zenodo/Data/NCBI_plus_AllTheBacteria_allelecall_results'  -o '/path/to/files/output_folder/IdentifySpuriousGenesSchema' -m schema -pm alleles_vs_alleles --t 4 -c 6

- Replace `/path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the identified spurious genes (The results are also available in the zenodo files).

Conclusion
----------

You have successfully identified spurious genes in a genomic schema using the `IdentifySpuriousGene` module for `Schema`.

For more information on the `IdentifySpuriousGene` module, refer to the :doc:`IdentifySpuriousGene documentation <SchemaRefinery/Tutorials/IdentifySpuriousGeneSchema>`.