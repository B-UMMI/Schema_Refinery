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

2. Follow the following steps: :docs:`DownloadAssemblies tutorial <DownloadAssembliesTutorial <SchemaRefinery/Tutorials/DownloadAssembliesTutorial>`.

3. Based on the downloaded assemblies choose those that you want to use as schema seed (e.g best quality, most complete, etc.), create a schema using the `CreateSchema` module from chewBBACA.

.. code-block:: bash

    chewBBACA CreateSchema -i /path/to/input_folder -o /path/to/output_folder -t 4

- Replace `/path/to/input_folder` with the path to the folder containing the downloaded assemblies.
- Replace `/path/to/output_folder` with the path to the output folder.

.. Note:: The `CreateSchema` module will generate a `schema_seed` folder containing the schema seed.

4. Populate the schema seed with the downloaded assemblies using the `AlleleCall` module from chewBBACA.

.. code-block:: bash

    chewBBACA AlleleCall -i /path/to/schema_seed -g /path/to/genome_folder -o /path/to/output_folder -t 4

5. Follow the following steps: :docs:`IdentifySpuriousGenes Unclassified CDS tutorial <SchemaRefinery/Tutorials/IdentifySpuriousGenesUnclassifiedCDS>`.

- In, normal workflow the users would ave to select the best loci to keep based on the recomendations of the `IdentifySpuriousGenes` module. Here we skip this step to show the full workflow.

6. Follow the following steps: :doc:`AdaptLoci tutorial <SchemaRefinery/Tutorials/AdaptLociTutorial>`.

- Pass as input the temp_fastas_path.txt, that has the paths for temp_fastas folder generated in the `IdentifySpuriousGenes` module.

7. Follow the following steps: :doc:`IdentifySpuriousGenes tutorial <SchemaRefinery/Tutorials/IdentifySpuriousGenesTutorialSchema>`.

- In, normal workflow the users would ave to select the best loci to keep based on the recomendations of the `IdentifySpuriousGenes` module. Here we skip this step to show the full workflow.

Optional modules to further refine or create a schema:
------------------------------------------------------

8. Follow the following steps: :doc:`AdaptLoci tutorial <SchemaRefinery/Tutorials/AdaptLociTutorial>`.

9. Follow the following steps: :doc:`MatchSchemas tutorial <SchemaRefinery/Tutorials/MatchSchemasTutorial>`.

10. Follow the following steps: :doc:`SchemaAnnotation tutorial <SchemaRefinery/Tutorials/SchemaAnnotationTutorial>`.

11. Follow the following steps: :doc:`IdentifyParalogousLoci tutorial <SchemaRefinery/Tutorials/IdentifyParalogousLociTutorial>`.