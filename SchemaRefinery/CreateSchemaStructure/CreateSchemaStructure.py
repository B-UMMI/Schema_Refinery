import os
import shutil
from typing import Dict, List


try:
    from AdaptLoci import AdaptLoci
    from utils import (sequence_functions as sf,
                                            file_functions as ff,
                                            logger_functions as logf,
                                            print_functions as pf,
                                            globals as gb)
    
except ModuleNotFoundError:
    from SchemaRefinery.AdaptLoci import AdaptLoci
    from SchemaRefinery.utils import (sequence_functions as sf,
                                                            file_functions as ff,
                                                            logger_functions as logf,
                                                            print_functions as pf,
                                                            globals as gb)
    

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

    output_d= os.path.abspath(output_directory)

    # Get all FASTA paths in the FASTA folder
    fastas_files: Dict[str, str] = {
        os.path.basename(fasta_file).split('.')[0]: os.path.join(fastas_folder, fasta_file)
        for fasta_file in os.listdir(fastas_folder)
    }
    action_list: Dict[int, Dict[str, List[str]]] = {}
    action_id: int = 1
    ids_list: List[str] = []
    last_rec: str = None

    new_fastas_path: List[str] = []
    
    temp_fasta_folder = os.path.join(output_d, 'temp_fasta')
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
                ids_list = []
                continue
            # Split the line into the action and the IDs   
            id, recommendation = line.split('\t')
            # Check if the recommendation is different from the preivous one
            # If so, start a new set of IDs
            if recommendation != last_rec:
                # Split the IDs into a list
                ids_list = []
                last_rec = recommendation
            # Save the action and the IDs in the action_list dictionary
            ids_list.append(id)
            action_list.setdefault(action_id, {}).update({recommendation: ids_list})

    processed_files: List[str] = []
    # For each action in the action_list dictionary
    # Recomendations can be 'Joined', 'Choice', 'Drop' or 'Add'
    for action_id, recommendations in action_list.items():
        # For each recommendation in the action dictionary
        for recommendation, ids_list in recommendations.items():
            # If the recommendation is 'Joined'
            if "Join" in recommendation:
                output_file: str = os.path.join(temp_fasta_folder, f'{ids_list[0]}.fasta')
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
                                    out.write(f'>{ids_list[0]}_{allele_id}\n{seq}\n')
                                    # Increment the allele_id and add the FASTA hash to the seen_fastas list
                                    allele_id += 1
                                    seen_fastas.append(fasta_hash)
                                else:
                                    continue

                            pf.print_message(f'File {id_} added to {ids_list[0]} at {output_file}', "info")
                        else:
                            pf.print_message(f'File {id_} not found in the FASTA folder', "info")
            # If the recommendation is 'Choice' or 'Add'
            elif "Choice" in recommendation or "Add" in recommendation:
                for id_ in ids_list:
                    processed_files.append(id_) # Add the ID to the processed_files list
                    output_file = os.path.join(temp_fasta_folder, f'{id_}.fasta')
                    # Append the new FASTA file path to the new_fastas_path list
                    new_fastas_path.append(output_file)
                    fasta_file= fastas_files[id_]  # Get the FASTA file path
                    shutil.copy(fasta_file, output_file)
                    pf.print_message(f'File {id_} copied to {output_file}', "info")
            else:
                processed_files.extend(ids_list) # Add the IDS to the processed_files list
                pf.print_message(f"The following IDs: {', '.join(ids_list)} have been removed due to drop action", "info")

    # Create schema structure
    pf.print_message("Create Schema Structure...", "info")
    # Schema path
    schema_path = os.path.join(output_d, 'schema')
    AdaptLoci.adapt_loci(temp_fasta_folder, schema_path, cpu, bsr, translation_table)

    if not no_cleanup:
        pf.print_message("\nCleaning up temporary files...", "info")
        ff.cleanup(output_d, [schema_path, logf.get_log_file_path(gb.LOGGER)])