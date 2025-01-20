import os
import shutil
from typing import Dict, List


try:
    from utils import (sequence_functions as sf,
                       file_functions as ff,
                       logger_functions as logf,
                       globals as gb)
    from AdaptLoci import AdaptLoci
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf,
                                      file_functions as ff,
                                      logger_functions as logf,
                                      globals as gb)
    from SchemaRefinery.AdaptLoci import AdaptLoci

def create_schema_structure(recommendations_file: str, 
                            fastas_folder: str,
                            output_directory: str,
                            cpu: int,
                            bsr: float,
                            translation_table: int,
                            no_cleanup:bool,) -> None:
    """
    Creates a schema structure based on the recommendations provided in the recommendations file.

    Parameters
    ----------
    recommendations_file : str
        Path to the file containing the recommendations.
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
    action_id: int = 1

    new_fastas_path: List[str] = []
    
    temp_fasta_folder = os.path.join(output_directory, 'temp_fasta')
    ff.create_directory(temp_fasta_folder)
    # Read and process the recommendations file
    with open(recommendations_file, 'r') as f:
        for index, line in enumerate(f):
            # Skip the Header
            if index == 0:
                continue
            # Strip any leading/trailing whitespace characters
            line = line.strip()
            # Skip empty lines and lines starting with '#'
            if line == '#':
                action_id += 1
                continue
            # Split the line into the action and the IDs
            recommendation, ids = line.split('\t')
            # Split the IDs into a list
            ids_list: List[str] = ids.split(',')
            # Save the action and the IDs in the action_list dictionary
            action_list.setdefault(action_id, {}).update({recommendation: ids_list})

    processed_files: List[str] = []
    # For each action in the action_list dictionary
    # Recomendations can be 'Joined', 'Choice' or 'Drop'
    for action_id, recommendations in action_list.items():
        # For each recommendation in the action dictionary
        for recommendation, ids_list in recommendations.items():
            # If the recommendation is 'Joined'
            if "Joined" in recommendation:
                
                new_file_name: str = recommendation.split('_')[1]  # Get the new file name
                output_file: str = os.path.join(temp_fasta_folder, f'{new_file_name}.fasta')
                # Append the new FASTA file path to the new_fastas_path list
                new_fastas_path.append(output_file)
                # Write the new FASTA file with the desired outcome
                with open(output_file, 'w') as out:
                    allele_id: int = 1  # Initialize the allele_id
                    seen_fastas: List[str] = []  # Initialize the seen_fastas list that stores FASTA hashes
                    # For each ID in the ids_list
                    for id_ in ids_list:
                        if id_ in fastas_files:  # If the ID is in the fastas_files dictionary
                            processed_files.append(id_) # Add the ID to the processed_files list
                            fasta_file: str = fastas_files[id_]  # Get the FASTA file path
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
            # If the recommendation is 'Choice'
            elif "Choice" in recommendation:
                for id_ in ids_list:
                    processed_files.append(id_) # Add the ID to the processed_files list
                    output_file = os.path.join(temp_fasta_folder, f'{id_}.fasta')
                    # Append the new FASTA file path to the new_fastas_path list
                    new_fastas_path.append(output_file)
                    fasta_file= fastas_files[id_]  # Get the FASTA file path
                    shutil.copy(fasta_file, output_file)
                    print(f'File {id_} copied to {output_file}')
            else:
                processed_files.extend(ids_list) # Add the IDS to the processed_files list
                print(f"The following IDs: {', '.join(ids_list)} have been removed due to drop action")

    print("Adding all of the remaining files that were not dropped or had action to do.")
    # Remove from fasta _files the processed IDs
    for id_ in processed_files:
        if id_ in fastas_files:
            fastas_files.pop(id_)
    # For each ID in the fastas_files dictionary
    for id_, fasta_file in fastas_files.items():
        output_file = os.path.join(temp_fasta_folder, f'{id_}.fasta') # Create the output file path
        # Append the new FASTA file path to the new_fastas_path list
        new_fastas_path.append(output_file)
        # Copy the FASTA file to the output file
        shutil.copy(fasta_file, output_file) # Copy the FASTA file to the output file
        print(f'File {id_} copied to {output_file}')

    fastas_paths_file: str = os.path.join(output_directory, 'fastas_paths.txt')

    with open(fastas_paths_file, 'w') as out:
        for path in new_fastas_path:
            out.write(f'{path}\n')

    # Create schema structure
    print("Create Schema Structure...")
    # Schema path
    schema_path = os.path.join(output_directory, 'schema')
    ff.create_directory(schema_path)
    AdaptLoci.main(fastas_paths_file, schema_path, cpu, bsr, translation_table)

    if not no_cleanup:
        print("\nCleaning up temporary files...")
        ff.cleanup(output_directory, [schema_path, logf.get_log_file_path(gb.LOGGER)])