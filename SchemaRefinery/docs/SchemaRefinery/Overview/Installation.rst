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

2. **Clone the Repository**: Clone the `SchemaRefinery` repository from GitHub:

    .. code-block:: bash

        git clone https://github.com/MForofontov/Schema_Refinery.git

3. **Change Directory**: Navigate to the cloned repository:

    .. code-block:: bash

        cd Schema_Refinery

4. **Create a Virtual Environment**: It is recommended to create a virtual environment to manage dependencies:

    .. code-block:: bash

        python -m venv venv
        source venv/bin/activate  # On Windows use `venv\Scripts\activate`

5. **Install Dependencies**: Install the required Python packages:

    .. code-block:: bash

        pip install -r requirements.txt

6. **Install the Package**: Install the `SchemaRefinery` package:

    .. code-block:: bash

        python setup.py install

7. **Verify Installation**: Verify the installation by running the following command:

    .. code-block:: bash

        SR --help

8. **Deactivate the Virtual Environment**: Once you are done, you can deactivate the virtual environment:

    .. code-block:: bash

        deactivate

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