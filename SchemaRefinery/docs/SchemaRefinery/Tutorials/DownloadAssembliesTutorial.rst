DownloadAssemblies - Tutorial
=============================

Objective
---------

This tutorial will show you how to download assemblies from NCBI and ENA661K using the `DownloadAssemblies` module.

Prerequisites
-------------
- SchemaRefinery installed
- Python 3.9 or higher
- Requests library (`pip install requests`)
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Procedure
---------

1. Open a terminal window.

2. Run the following command to download assemblies from NCBI and ENA661K:

.. code-block:: bash

    SR DownloadAssemblies -t "Mycoplasma pneumonia" -db NCBI ENA661K -o /path/to/output_folder -e email@example -th 4 -fm --download

- Replace `email@example` with your email address.
- Replace `/path/to/output_folder` with the path to the output directory.

3. Press Enter to execute the command.

4. Wait for the download process to complete.

5. Check the output directory for the downloaded assemblies.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :ref:`DownloadAssemblies documentation <Download_assemblies>`.

Conclusion
----------

You have successfully downloaded assemblies from NCBI and ENA661K using the `DownloadAssemblies` module.

For more information on the `DownloadAssemblies` module, refer to the :ref:`DownloadAssemblies documentation <Download_assemblies>`.
