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
- Download the zenodo file from the link: https://zenodo.org/record/5560007/files/genbank_files.zip

Procedure
---------

1. Open a terminal window.

2. Run the following command to identify paralogous loci in a schema:

.. code-block:: bash

    SR IdentifyParalagousLoci -s 'path/to/files/zenodo/Data/mpneumoniae_schema/mpneumoniae_schema' -o 'path/to/files/output_folder/IdentifyParalogousLoci_Results' -tt 4 -c 6 -pm alleles_vs_alleles

- Replace `path/to/files/` with the actual path to the files.

3. Press Enter to execute the command.

4. Wait for the process to complete.

5. Check the output directory for the identified paralogous loci (The results are also available in the zenodo files).

Conclusion
----------

You have successfully identified paralogous loci in a schema using the `IdentifyParalogousLoci` module.

For more information on the `IdentifyParalogousLoci` module, refer to the :doc:`IdentifyParalogousLoci documentation <SchemaRefinery/Tutorials/IdentifyParalogousLociTutorial>`.
