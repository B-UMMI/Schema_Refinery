SchemaRefinery - Full Tutorial
==============================

Objective
---------

This tutorial will guide you through the complete workflow of SchemaRefinery, from schema creation to schema refinement using the `SchemaRefinery` modules.

Prerequisites
-------------
- SchemaRefinery installed
- chewBBACA 3.3.10 or higher
- Python 3.9 or higher
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)
- Requests library (`pip install requests`)

Procedure
---------

1. Open a terminal window.

2. Follow the following steps: :ref:`DownloadAssemblies tutorial <DownloadAssembliesTutorial>`.

3. Based on the downloaded assemblies choose those that you want to use as schema seed (e.g best quality, most complete, etc.), create a schema using the `CreateSchema` module from chewBBACA.

.. code-block:: bash

    chewBBACA CreateSchema -i /path/to/input_folder -o /path/to/output_folder -t 4

.. Note:: The `CreateSchema` module will generate a `schema_seed` folder containing the schema seed.

4. Populate the schema seed with the downloaded assemblies using the `AlleleCall` module from chewBBACA.

.. code-block:: bash

    chewBBACA AlleleCall -i /path/to/schema_seed -g /path/to/genome_folder -o /path/to/output_folder -t 4


