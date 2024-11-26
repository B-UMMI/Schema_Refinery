SchemaRefinery - A Tool for Refining Genomic Schemas
====================================================

Description
-----------

The `SchemaRefinery` repository contains tools and modules for refining genomic schemas. These tools help in identifying paralogous loci, spurious genes, and annotating schemas. The repository supports various genomic data processing tasks and provides configurable parameters for different processes.

Modules
-------

The repository includes the following main modules:

1. **:ref:`IdentifyParalogousLoci <identify_paralogous_loci>`**: Identifies paralogous loci in a schema.
2. **:ref:`IdentifySpuriousGenes <identify_spurious_genes>`**: Identifies spurious genes in a schema.
3. **:ref:`SchemaAnnotation <schema_annotation>`**: Annotates schemas with additional information.
4. **:ref:`MatchSchemas <match_schemas>`**: Matches schemas in a directory.
5. **:ref:`DownloadAssemblies <download_assemblies>`**: Downloads genomic assemblies from various databases.
6. **:ref:`AdaptLoci <adapt_loci>`**: Adapts loci in a schema.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Modules
-------

Each module can be used independently by running the corresponding script with the required command-line arguments. Below are examples for each module:

### :ref:`IdentifyParalogousLoci <identify_paralogous_loci>`

.. code-block:: bash

    SR IdentifyParalogousLoci --help

### :ref:`IdentifySpuriousGenes <identify_spurious_genes>`

.. code-block:: bash

    SR IdentifySpuriousGenes --help

### :ref:`SchemaAnnotation <schema_annotation>`

.. code-block:: bash

    SR SchemaAnnotation --help

### :ref:`MatchSchemas <match_schemas>`

.. code-block:: bash

    SR MatchSchemas --help

### :ref:`DownloadAssemblies <download_assemblies>`

.. code-block:: bash

    SR DownloadAssemblies --help

### :ref:`AdaptLoci <adapt_loci>`

.. code-block:: bash

    SR AdaptLoci --help
