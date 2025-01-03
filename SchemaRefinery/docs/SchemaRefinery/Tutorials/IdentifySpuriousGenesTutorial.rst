IdentifySpuriousGenes - Tutorial
===============================

Objective
---------

This tutorial will guide you through the process of using the `IdentifySpuriousGenes` module to identify spurious genes in genomic data.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to identify spurious genes in genomic data:

.. code-block:: bash

    SR IdentifySpuriousGenes -i /path/to/input_folder -o /path/to/output_folder -a /path/to/allele_call_folder -t 4

- Replace `/path/to/input_folder` with the path to the input folder containing the genomic data.
- Replace `/path/to/output_folder` with the path to the output folder.
- Replace `/path/to/allele_call_folder` with the path to the allele call folder.

3. Press Enter to execute the command.

4. Wait for the identification process to complete.

5. Check the output folder for the identified spurious genes.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :docs:`IdentifySpuriousGenes documentation <SchemaRefinery/Modules/IdentifySpuriousGenes>`.

Conclusion
----------

Based on output results the user can decide what loci to remove and keep based on the recomendations of the `IdentifySpuriousGenes` module for the Schema.

For more information on the `IdentifySpuriousGenes` module, refer to the :docs:`IdentifySpuriousGenes documentation <SchemaRefinery/Modules/IdentifySpuriousGenes>`.
