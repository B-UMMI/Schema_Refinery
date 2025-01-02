Installation Guide
==================

Follow these steps to install the `SchemaRefinery` package on your system.

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

6. **Install Dependencies**: Install the required Python packages:

    .. code-block:: bash

        conda install blast
        pip install -r requirements.txt

7. **Install the Package**: Install the `SchemaRefinery` package:

    .. code-block:: bash

        python setup.py install

8. **Verify Installation**: Verify the installation by running the following command:

    .. code-block:: bash

        SR --help

9. **Deactivate the Conda Environment**: Once you are done, you can deactivate the conda environment:

    .. code-block:: bash

        conda deactivate

Additional Information
----------------------

- **Updating the Package**: To update the `SchemaRefinery` package, navigate to the repository directory and pull the latest changes:

    .. code-block:: bash

        cd Schema_Refinery
        git pull
        python setup.py install

- **Uninstalling the Package**: To uninstall the `SchemaRefinery` package, use the following command:

    .. code-block:: bash

        pip uninstall SchemaRefinery

- **Troubleshooting**: If you encounter any issues during installation, ensure that all dependencies are installed and that you are using a compatible version of Python. You can also refer to the [GitHub repository](https://github.com/MForofontov/Schema_Refinery) for more information and support.