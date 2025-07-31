SchemaRefinery - Full Tutorial
==============================

Objective
---------

This tutorial will guide you through the apossible workflow of SchemaRefinery, from schema creation to schema refinement using the `SchemaRefinery` modules.

Prerequisites
-------------
- SchemaRefinery installed
- chewBBACA 3.3.10 or higher
- Python between 3.9 and 3.11
- `NCBI datasets <https://www.ncbi.nlm.nih.gov/datasets/>`_

Procedure
---------

1. Open a terminal window.

2. Download the assemblies from NCBI needed for creating the schema:
    .. code-block:: bash
        SR DownloadAssemblies -f path/to/input_tsv_file_with_taxon -db NCBI -e youremail@example.com -o path/to/DownloadAssemblies_NCBI_download -fm --download

Based on the downloaded assemblies choose those that you want to use as schema seed (e.g. the best quality, most complete, etc.).

3. Create a schema using the `CreateSchema` module from chewBBACA:
    .. code-block:: bash
        chewBBACA.py CreateSchema -i /path/to/DownloadAssemblies_NCBI_download/assemblies_ncbi_unziped -o /path/to/CreateSchema_chewbbaca_mySchema -t 4

For more information on the `CreateSchema` module, refer to the `chewBBACA documentation <https://chewbbaca.readthedocs.io/en/latest/user/modules/CreateSchema.html>`_.

.. Note:: 
    The `CreateSchema` module will generate a `schema_seed` folder containing the schema seed.

4. Populate the schema seed with the downloaded assemblies using the `AlleleCall` module from chewBBACA.
    .. code-block:: bash
        chewBBACA.py AlleleCall -i /path/to/CreateSchema_chewbbaca_mySchema/schema_seed -g /path/to/DownloadAssemblies_NCBI_download/assemblies_ncbi_unziped -o /path/to/AlleleCall_folder -t 4

For more information on the `AlleleCall` module, refer to the `chewBBACA documentation <https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCall.html>`_.

5. Check the unclassified CDS for spurious genes and create a proto schema from it:
    .. code-block:: bash
        SR IdentifySpuriousGenes -s path/to/CreateSchema_chewbbaca_mySchema/schema_seed -a path/to/AlleleCall_folder -m unclassified_cds -o path/to/IdentifySpuriousGenes_uCDS_mySchema -c 6

In a normal workflow the users would have to select the best loci to keep based on the recomendations of the `IdentifySpuriousGenes` module. Here we skip this step to show the full workflow.

6. Adapt the proto schema created with the unclssified CDS:
    .. code-block:: bash
        SR AdaptLoci -i path/to/IdentifySpuriousGenes_uCDS_mySchema/temp_fastas -o path/to/AdaptLoci_unclassified

Pass as input the ttemp_fastas folder generated in the `IdentifySpuriousGenes` module. Repeat the step 4 with this new schema to create a new AlleleCall folder. 

7. Refine the schema created with the unclassified CDS:
    .. code-block:: bash 
        SR IdentifySpuriousGenes -s path/to/AdaptLoci_unclassified/schema_seed -a path/to/AlleleCall_unclassified -m schema -o path/to/IdentifySpuriousGenes_unclassifiedSchema -c 6

Analyse the clusters formed and change the "Choice" actions into "Join", "Add" or "Drop".

8. Create a final schema using the `altered` **recommendations_annotations.tsv** file from the previous step:
    .. code-block:: bash 
        SR CreateSchemaStructure -s path/to/AdaptLoci_unclassified/schema_seed -rf path/to/IdentifySpuriousGenes_unclassifiedSchema/recommendations_annotations.tsv -o path/to/CreateSchemaStructure_refined_schema -c 6

Use the altered recomendations file in the `CreateSchemaStructure` module folder, as that one has a selection of the loci to be dropped or added.

Optional modules to further refine or create a schema:
------------------------------------------------------

8. Follow the following steps: `MatchSchemas tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/MatchSchemasTutorial.html>`_.

- Matches two different schema loci.

9. Follow the following steps: `SchemaAnnotation tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/SchemaAnnotationTutorial.html>`_.

- Annotates the schema with additional information from various databases.

10. Follow the following steps: `IdentifyParalogousLoci tutorial <https://schema-refinery.readthedocs.io/en/latest/SchemaRefinery/Tutorials/IdentifyParalogousLociTutorial.html>`_.

- Identifies paralogous loci in a schema.

.. Note:: The assemblies present in NCBI may change, so the results may vary.

Conclusion
----------

You have successfully completed a possible workflow of SchemaRefinery, from schema creation to schema refinement using the `SchemaRefinery` modules.

For more information on the `SchemaRefinery` modules, refer to the `SchemaRefinery documentation <https://schema-refinery.readthedocs.io/en/latest/index.html>`_.