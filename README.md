# SchemaRefinery - A Tool for Refining Genomic Schemas

## Description

The `SchemaRefinery` repository contains tools and modules for refining genomic schemas. These tools help in identifying paralogous loci, spurious genes, and annotating schemas. The repository supports various genomic data processing tasks and provides configurable parameters for different processes.

# Installation Guide

Follow these steps to install the `SchemaRefinery` package on your system.

1. **Install Git**: Ensure that Git is installed on your system. You can install Git using the following command:

    ```bash
    # For macOS
    brew install git

    # For Ubuntu/Debian
    sudo apt-get install git

    # For Fedora
    sudo dnf install git
    ```

2. **Install Conda**: Ensure that Conda is installed on your system. You can install Miniconda (a minimal Conda installer) using the following commands:

    ```bash
    # For macOS and Linux
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

3. **Clone the Repository**: Clone the [SchemaRefinery](http://_vscodecontentref_/0) repository from GitHub:

    ```bash
    git clone https://github.com/MForofontov/Schema_Refinery.git
    ```

4. **Change Directory**: Navigate to the cloned repository:

    ```bash
    cd Schema_Refinery
    ```

5. **Create a Conda Environment**: It is recommended to create a conda environment to manage dependencies:

    ```bash
    conda create --name schema_refinery python=3.9
    conda activate schema_refinery
    ```

6. **Install Dependencies**: Install the required Python packages:

    ```bash
    conda install blast
    pip install -r requirements.txt
    ```

7. **Install the Package**: Install the [SchemaRefinery](http://_vscodecontentref_/1) package:

    ```bash
    python setup.py install
    ```

8. **Verify Installation**: Verify the installation by running the following command:

    ```bash
    SR --help
    ```

9. **Deactivate the Conda Environment**: Once you are done, you can deactivate the conda environment:

    ```bash
    conda deactivate
    ```

## Modules

The repository includes the following main modules:

1. **IdentifyParalogousLoci**: Identifies paralogous loci in a schema.
2. **IdentifySpuriousGenes**: Identifies spurious genes in a schema.
3. **SchemaAnnotation**: Annotates schemas with additional information.
4. **MatchSchemas**: Matches schemas in a directory.
5. **DownloadAssemblies**: Downloads genomic assemblies from various databases.
6. **AdaptLoci**: Adapts loci in fasta format to a schema format.

## Dependencies

- Python 3.9 or higher
- Biopython library (`pip install biopython`)
- NCBI datasets ([NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/))

## Modules Usage

Each module can be used independently by running the corresponding script with the required command-line arguments. Below are examples for each module:

### IdentifyParalogousLoci

    ```bash
    SR IdentifyParalogousLoci --help
    ```

### IdentifySpuriousGenes

    ```bash
    SR IdentifySpuriousGenes --help
    ```

### SchemaAnnotation

    ```bash
    SR SchemaAnnotation --help
    ```

### MatchSchemas

    ```bash
    SR MatchSchemas --help
    ```

### DownloadAssemblies

    ```bash
    SR DownloadAssemblies --help
    ```

### AdaptLoci

    ```bash
    SR AdaptLoci --help
    ```

## Troubleshooting


If you encounter issues while using the modules, consider the following troubleshooting steps:

- Verify that the paths to the schema, output, and other directories are correct.
- Check the output directory for any error logs or messages.
- Increase the number of CPUs using the `-c` or `--cpu` option if the process is slow.
- Ensure that you have a stable internet connection.

if the issue persists, please report it to the development team using github issues.

## Contributing

We welcome contributions to the SchemaRefinery project. If you would like to contribute, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them with a clear message.
4. Push your changes to your forked repository.
5. Create a pull request to the main repository.

## License

This project is licensed under the MIT License. See the [LICENSE](https://opensource.org/license/mit) file for details.

## Contact Information

For support or to report issues, please contact the development team at GitHub issues in [SchemaRefinery GitHub repository](https://github.com/B-UMMI/Schema_Refinery).
