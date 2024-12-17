SchemaRefinery - A Tool for Refining Genomic Schemas
====================================================

Description
-----------

The `SchemaRefinery` is a comprehensive toolkit designed for refining and managing genomic schemas. It provides a suite of modules for various tasks such as identifying paralogous loci, downloading genomic assemblies, and adapting loci into standardized schemas.

Dependencies
------------

- Python 3.9 or higher
- Biopython library (`pip install biopython`)
- NCBI datasets (`https://www.ncbi.nlm.nih.gov/datasets/ <https://www.ncbi.nlm.nih.gov/datasets/>`_)

Installation
------------
Current implementation of `SchemaRefinery` is available as a Python package available at `SchemaRefinery GitHub repository <https://github.com/B-UMMI/Schema_Refinery>`_. To install the package, run the following command:

1. **Install Python**: Ensure you have Python 3.9 or higher installed. You can download Python from the official website: https://www.python.org/downloads/

2. **Install Git**: If you don't have Git installed, you can install it using the following command:

    .. code-block:: bash

        sudo apt-get install git # For Debian-based systems
        sudo yum install git # For Red Hat-based systems
        brew install git # For macOS

3. **Clone the Repository**: Clone the `SchemaRefinery` repository from GitHub:

    .. code-block:: bash

        git clone https://github.com/MForofontov/Schema_Refinery.git # Clone the repository

4. **Change Directory**: Navigate to the cloned repository:

    .. code-block:: bash

        cd Schema_Refinery # Change directory to the repository

5. **Install Dependencies**: Install the required Python packages:

    .. code-block:: bash

        pip install biopython

6. **Install the Package**: Install the `SchemaRefinery` package:
    
    .. code-block:: bash

        python setup.py install

7. **Verify Installation**: Verify the installation by running the following command:
    
        .. code-block:: bash
    
            SR --help

Modules
-------

The repository includes the following main modules:

1. **:ref:`IdentifyParalogousLoci <Identify_paralogous_loci>`**: Identifies paralogous loci in a schema.
2. **:ref:`IdentifySpuriousGenes <Identify_spurious_genes>`**: Identifies spurious genes in a schema.
3. **:ref:`SchemaAnnotation <Schema_annotation>`**: Annotates schemas with additional information.
4. **:ref:`MatchSchemas <Match_schemas>`**: Matches schemas in a directory.
5. **:ref:`DownloadAssemblies <Download_assemblies>`**: Downloads genomic assemblies from various databases.
6. **:ref:`AdaptLoci <Adapt_loci>`**: Adapts loci in fasta format to a schema format.

Modules Usage
-------------

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

Schema Creation Workflow
------------------------

Workflow for creating a schema using the `SchemaRefinery` modules and chewBBACA:

.. image:: source/Schema_creation_workflow.png
   :alt: Schema Creation Workflow
   :width: 80%
   :align: center

Troubleshooting
---------------

If you encounter issues while using the modules, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and other directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
- Ensure that you have a stable internet connection.

if the issue persists, please report it to the development team using github issues.

Contributing
------------

We welcome contributions to the SchemaRefinery project. If you would like to contribute, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them with a clear message.
4. Push your changes to your forked repository.
5. Create a pull request to the main repository.

License
-------

This project is licensed under the MIT License. See the `LICENSE <https://opensource.org/license/mit>`_ file for details.

Contact Information
-------------------

For support or to report issues, please contact the development team at github issues in `SchemaRefinery GitHub repository <https://github.com/B-UMMI/Schema_Refinery>`_.