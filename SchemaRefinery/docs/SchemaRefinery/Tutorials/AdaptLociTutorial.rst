AdaptLoci - Tutorial
====================

Objective
---------

This tutorial will guide you through the process of using the `AdaptLoci` module to adapt loci from FASTA files into a chewBBACA-compatible schema.

Prerequisites
-------------

- SchemaRefinery installed
- Python 3.9 or higher
- Biopython library (`pip install biopython`)

Procedure
---------

1. Open a terminal window.

2. Run the following command to adapt loci from FASTA files into a chewBBACA-compatible schema:

.. code-block:: bash

    SR AdaptLoci -i 'path/to/files/zenodo/Data/fastas_path.txt'  -o 'path/to/files/AdaptLoci_Results' -tt 4 -c 6

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the adaptation process to complete.

5. Check the output folder for the adapted loci.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`AdaptLoci documentation <SchemaRefinery/Modules/AdaptLoci>`.

Conclusion
----------

You have successfully adapted loci from FASTA files into a chewBBACA-compatible schema using the `AdaptLoci` module.

For more information on the `AdaptLoci` module, refer to the :doc:`AdaptLoci documentation <SchemaRefinery/Modules/AdaptLoci>`.