Installation Guide
==================

Follow these steps to install `SchemaRefinery`. There are three methods to install `SchemaRefinery`: Bioconda, PyPI or using Git. Choose the one that is most convenient based on your system and preferences. Schema Refinery requires Python 3.9 or higher, BLAST+ and NCBI datasets to be installed.

.. important::
	For a better organization of the packages and dependencies, it is advised to create a new conda environment to install Schema Refinery into.

Installation through Bioconda
----------------------------

1. **Install Conda**: Ensure that Conda is installed on your system. You can install Miniconda (a minimal Conda installer) using the following commands:

    .. code-block:: bash

        # For macOS and Linux
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

2. **Add Bioconda and conda-forge channels**: Ensure that conda can search and install packages from the Bioconda and conda-forge channels.

    .. code-block:: bash

        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda config --set channel_priority strict

3. **Create a Conda Environment**: It is recommended to create a conda environment to manage dependencies:

    .. code-block:: bash

        conda create --name schema_refinery python=3.9
        conda activate schema_refinery

4. **Install Schema Refinery**:

    .. code-block:: bash

        conda install bioconda::schemarefinery


Installation through PyPi
--------------------------

1. **Install pip**: For Python 3.4 and higher installed using Windows or macOS `pip` is installed as default. For Linux environments run:

    .. code-block:: bash

        #For Ubuntu and python3
        sudo apt-get install python3 - pip

        #For CentOS and python3
        sudo yum install python3 - pip

2. **Install Conda**: Ensure that Conda is installed on your system. You can install Miniconda (a minimal Conda installer) using the following commands:

    .. code-block:: bash

        # For macOS and Linux
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

3. **Create a Conda Environment**: It is recommended to create a conda environment to install Schema Refinery and manage dependencies.

    .. code-block:: bash

        conda create --name schema_refinery python=3.9
        conda activate schema_refinery

4. **Install SchemaRefinery using pip**:

    .. code-block:: bash

        pip install SchemaRefinery


Installation through GitHub
---------------------------

1. **Install Git**: Ensure that Git is installed on your system. You can install Git using the following command:

    .. code-block:: bash

        # For macOS
        brew install git

        # For Ubuntu/Debian
        sudo apt-get install git

        # For Fedora
        sudo dnf install git

2. **Install Conda**: Ensure that Conda is installed on your system. You can install Miniconda (a minimal Conda installer) using the following commands:

    .. code-block:: bash

        # For macOS and Linux
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

3. **Clone the Repository**: Clone the `SchemaRefinery` repository from GitHub:

    .. code-block:: bash

        git clone https://github.com/MForofontov/Schema_Refinery.git

4. **Change Directory**: Navigate to the cloned repository:

    .. code-block:: bash

        cd Schema_Refinery

5. **Create a Conda Environment**: It is recommended to create a conda environment to manage dependencies:

    .. code-block:: bash

        conda create --name schema_refinery python=3.9
        conda activate schema_refinery

6. **Install Dependencies**: Install BLAST and the required Python packages:

    .. code-block:: bash

        conda install blast
        pip install -r requirements.txt

7. **Install Schema Refinery**:

    .. code-block:: bash

        python setup.py install

8. **Verify Installation**: Verify the installation by running the following command:

    .. code-block:: bash

        SR --help

Additional Information
----------------------

- **Updating the Package**: To update `SchemaRefinery`, navigate to the repository directory and pull the latest changes:

    .. code-block:: bash

        cd Schema_Refinery
        git pull
        python setup.py install

- **Uninstalling the Package**: To uninstall `SchemaRefinery`, use the following command:

    .. code-block:: bash

        pip uninstall SchemaRefinery

.. important::
	If you encounter any issues during installation, ensure that all dependencies are installed and that you are using a compatible version of Python. You can also refer to the `GitHub repository <https://github.com/B-UMMI/Schema_Refinery>`_ for more information and support.
