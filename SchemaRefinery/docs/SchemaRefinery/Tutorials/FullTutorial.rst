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
- (`NCBI datasets <https://www.ncbi.nlm.nih.gov/datasets/>`_)
- Requests library (`pip install requests`)

Procedure
---------

1. Open a terminal window.

2. Follow the following steps: `DownloadAssemblies tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/DownloadAssembliesTutorial.html>`.

3. Based on the downloaded assemblies choose those that you want to use as schema seed (e.g best quality, most complete, etc.), create a schema using the `CreateSchema` module from chewBBACA.

.. code-block:: bash

    chewBBACA CreateSchema -i /path/to/input_folder -o /path/to/output_folder -t 4

- Replace `/path/to/input_folder` with the path to the folder containing the downloaded assemblies.
- Replace `/path/to/output_folder` with the path to the output folder.
- For more information on the `CreateSchema` module, refer to the `chewBBACA documentation <https://chewbbaca.readthedocs.io/en/latest/user/modules/CreateSchema.html>`_.

.. Note:: The `CreateSchema` module will generate a `schema_seed` folder containing the schema seed.

4. Populate the schema seed with the downloaded assemblies using the `AlleleCall` module from chewBBACA.

.. code-block:: bash

    chewBBACA AlleleCall -i /path/to/schema_seed -g /path/to/genome_folder -o /path/to/output_folder -t 4

- Replace `/path/to/schema_seed` with the path to the `schema_seed` folder.
- Replace `/path/to/genome_folder` with the path to the folder containing the downloaded assemblies.
- Replace `/path/to/output_folder` with the path to the output folder.
- For more information on the `AlleleCall` module, refer to the `chewBBACA documentation <https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCall.html>`_.

5. Follow the following steps: `IdentifySpuriousGenes Unclassified CDS tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/IdentifySpuriousGenesUnclassifiedCDS.html>`.

- In, normal workflow the users would ave to select the best loci to keep based on the recomendations of the `IdentifySpuriousGenes` module. Here we skip this step to show the full workflow.

6. Follow the following steps: `AdaptLoci tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/AdaptLociTutorial.html>`.

- Pass as input the temp_fastas_path.txt, that has the paths for temp_fastas folder generated in the `IdentifySpuriousGenes` module.

7. Follow the following steps: `IdentifySpuriousGenes tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/IdentifySpuriousGenesSchema.html>`.

- In, normal workflow the users would ave to select the best loci to keep based on the recomendations of the `IdentifySpuriousGenes` module. Here we skip this step to show the full workflow.

Optional modules to further refine or create a schema:
------------------------------------------------------

8. Follow the following steps: `MatchSchemas tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/MatchSchemasTutorial.html>`.

- Matches two different schema loci.

9. Follow the following steps: `SchemaAnnotation tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/SchemaAnnotationTutorial.html>`.

- Annotates the schema with additional information from various databases.

10. Follow the following steps: `IdentifyParalogousLoci tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/IdentifyParalogousLociTutorial.html>`.

- Identifies paralogous loci in a schema.

Note
----

The assemblies present in NCBI may change, so the results may vary. For exact results as provided in the Zenodo files, use the provided assemblies.

Conclusion
----------

You have successfully completed the full workflow of SchemaRefinery, from schema creation to schema refinement using the `SchemaRefinery` modules.

For more information on the `SchemaRefinery` modules, refer to the `SchemaRefinery documentation <https://schema-refinery.readthedocs.io/en/latest/index.html>`_.