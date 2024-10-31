import os
from typing import Dict, List

try:
    from utils import (sequence_functions as sf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf)

def create_schema_structure(instructions_file: str, 
                            fastas_folder: str, 
                            skip_choices: bool, 
                            output_directory: str) -> None:
    """
    Creates a schema structure based on the instructions provided in the instructions file.

    Parameters
    ----------
    instructions_file : str
        Path to the file containing the instructions.
    fastas_folder : str
        Path to the folder containing the FASTA files.
    skip_choices : bool
        Whether to skip recommendations with 'Choice'.
    output_directory : str
        Path to the directory where the output files will be saved.

    Returns
    -------
    None
        The function writes the output files to the specified directory.
    """
    # Get all FASTA paths in the FASTA folder
    fastas_files: Dict[str, str] = {
        os.path.basename(fasta_file).split('.')[0]: os.path.join(fastas_folder, fasta_file)
        for fasta_file in os.listdir(fastas_folder)
    }
    action_list: Dict[int, Dict[str, List[str]]] = {}
    i: int = 1

    # Read and process the instructions file
    with open(instructions_file, 'r') as f:
        for line in f:
            # Strip any leading/trailing whitespace characters
            line = line.strip()
            # Skip empty lines and lines starting with '#'
            if line == '#':
                i += 1
                continue
            # Split the line into the action and the IDs
            recommendation, ids = line.split('\t')
            # Split the IDs into a list
            ids_list: List[str] = ids.split(',')
            # Save the action and the IDs in the action_list dictionary
            action_list.setdefault(i, {}).update({recommendation: ids_list})

    # For each action in the action_list dictionary
    for action, recommendations in action_list.items():
        # For each recommendation in the action dictionary
        for recommendation, ids_list in recommendations.items():
            # If the recommendation is 'Choice' and skip_choices is True, skip the recommendation
            if skip_choices and recommendation == "Choice":
                continue
            # If the recommendation is 'Joined'
            elif recommendation == "Joined":
                new_file_name: str = ids_list[0]
                output_file: str = os.path.join(output_directory, f'{new_file_name}.fasta')
                # Write the new FASTA file with the desired outcome
                with open(output_file, 'w') as out:
                    allele_id: int = 1  # Initialize the allele_id
                    seen_fastas: List[str] = []  # Initialize the seen_fastas list that stores FASTA hashes
                    # For each ID in the ids_list
                    for id in ids_list:
                        if id in fastas_files:  # If the ID is in the fastas_files dictionary
                            fasta_file: str = fastas_files[id]  # Get the FASTA file path
                            fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_file, out)  # Fetch the FASTA dictionary
                            # For each header and sequence in the FASTA dictionary
                            for header, seq in fasta_dict.items():
                                fasta_hash: str = sf.hash_sequence(seq)  # Get the hash of the sequence
                                # If the FASTA hash is not in the seen_fastas list
                                if fasta_hash not in seen_fastas:
                                    # Write the new header and sequence to the output file
                                    out.write(f'>{new_file_name}_{allele_id}\n{seq}\n')
                                    # Increment the allele_id and add the FASTA hash to the seen_fastas list
                                    allele_id += 1
                                    seen_fastas.append(fasta_hash)
                                else:
                                    continue

                            print(f'File {id} added to {new_file_name} at {output_file}')
                        else:
                            print(f'File {id} not found in the FASTA folder')
            else:
                print(f"The following IDs: {', '.join(ids_list)} have been removed due to drop action")
