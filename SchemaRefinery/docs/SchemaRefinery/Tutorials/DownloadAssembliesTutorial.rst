DownloadAssemblies - Tutorial
=============================

Objective
---------

This tutorial will show you how to download assemblies from NCBI and ENA661K using the `DownloadAssemblies` module.

Prerequisites
-------------
- SchemaRefinery installed
- Python between 3.9 and 3.11
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Procedure
---------

1. Open a terminal window.

2. Run the following command to download assemblies from NCBI and ENA661K:

.. code-block:: bash

    SR DownloadAssemblies -db NCBI -th 3 --download -fm -o 'path/to/files/DownloadAssemblies_Results' -f '/path/input_table_example.tsv' -e your_email@email.com

- The input table used is in the DownloadAssemblies folder.
- Replace `path/to/files/` with the actual path to the files.
- Replace `your_email@email.com` with your email address.

3. Press Enter to execute the command.

4. Wait for the download process to complete.

5. Check the output directory for the downloaded assemblies.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the `DownloadAssemblies documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/DownloadAssemblies.html>`_.

Conclusion
----------

You have successfully downloaded assemblies from NCBI and ENA661K using the `DownloadAssemblies` module.

For more information on the `DownloadAssemblies` module, refer to the `DownloadAssemblies documentation <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Modules/DownloadAssemblies.html>`_.