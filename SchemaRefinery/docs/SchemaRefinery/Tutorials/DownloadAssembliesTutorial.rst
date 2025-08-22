DownloadAssemblies - Tutorial
=============================

Objective
---------

This tutorial will show you how to download assemblies from NCBI and ENA661K using the :doc:`DownloadAssemblies </SchemaRefinery/Modules/DownloadAssemblies>` module.

Prerequisites
-------------

- NCBI datasets command-line too installed (instructions `here <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/>`_).

Procedure
---------

1. Open a terminal window.

2. Modify the following command and run it to download assemblies from NCBI and ENA661K:

.. code-block:: bash

    SR DownloadAssemblies -db NCBI -th 3 --download -fm -o 'path/to/files/DownloadAssemblies_Results' -f '/path/input_table_example.tsv' -e your_email@email.com

.. important::
	The file passed to `-f` is available `here`. Do not forget to replace the paths in the command with the paths to the files in your system. Provide your email address to `-e`.

3. Check the output directory for the downloaded assemblies.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`DownloadAssemblies documentation </SchemaRefinery/Modules/DownloadAssemblies>`.
