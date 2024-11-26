SchemaRefinery - A Tool for Refining Genomic Schemas
====================================================

Description
-----------

The `SchemaRefinery` repository contains tools and modules for refining genomic schemas. These tools help in identifying paralogous loci, spurious genes, and annotating schemas. The repository supports various genomic data processing tasks and provides configurable parameters for different processes.

Modules
-------

The repository includes the following main modules:

1. **IdentifyParalogousLoci**: Identifies paralogous loci in a schema.
2. **IdentifySpuriousGenes**: Identifies spurious genes in a schema.
3. **SchemaAnnotation**: Annotates schemas with additional information.
4. **MatchSchemas**: Matches schemas in a directory.
5. **DownloadAssemblies**: Downloads genomic assemblies from various databases.

Dependencies
------------

- Python 3.6 or higher
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Usage
-----

Each module can be used independently by running the corresponding script with the required command-line arguments. Below are examples for each module:

### IdentifyParalogousLoci

.. code-block:: bash

    SR IdentifyParalogousLoci --help

### IdentifySpuriousGenes

.. code-block:: bash

    SR IdentifySpuriousGenes --help

### SchemaAnnotation

.. code-block:: bash

    SR SchemaAnnotation --help

### MatchSchemas

.. code-block:: bash

    SR MatchSchemas --help

### DownloadAssemblies

.. code-block:: bash

    SR DownloadAssemblies --help

### AdaptLoci

.. code-block:: bash

    SR AdaptLoci --help
