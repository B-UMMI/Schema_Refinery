ApplyRecommendations - Tutorial
===============================

Objective
---------

This tutorial will guide you through the process of using the :doc:`ApplyRecommendations </SchemaRefinery/Modules/ApplyRecommendations>` module to create a new schema given a schema and a set of recommendations.

Prerequisites
-------------

- chewBBACA 3.3.10 or higher (`chewBBACA's installation instructions <https://chewbbaca.readthedocs.io/en/latest/user/getting_started/installation.html>`_).
- Download the schema from `chewBBACA's tutorial <https://github.com/B-UMMI/chewBBACA_tutorial/blob/master/expected_results/Schema_creation/tutorial_schema.zip>`_.

Procedure
---------

1. Open the terminal.
2. Modify the following command and run it to execute the `ApplyRecommendations` module:

.. code-block:: bash

    SR ApplyRecommendations -rf /path/to/recommendations.tsv -ff /path/to/tutorial_schema/schema_seed -o /path/to/ApplyRecommendations_output -c 6 --nocleanup

.. important::
	Do not forget to replace the paths in the command with the paths to the files in your system. The file passed to the `-rf` argument should be the `recommendation.tsv` file available `here`. This file was created using the `IdentifySpuriousGenes` module and the "Choice" actions were altered randomly just for demonstration purposes.

3. Check the output folder containing the adapted schema. Its structure should be compatible with chewBBACA. To check if the changes were applied, compare the FASTA files of the loci of the original schema with the ones in the new schema to verify if the changes were applied correctly.

Example Output Structure
------------------------

To see the expected output structure, refer to the "Outputs" section in the :doc:`ApplyRecommendations documentation</SchemaRefinery/Modules/ApplyRecommendations>`.
