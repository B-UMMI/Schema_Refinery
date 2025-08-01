import os
import shutil
import sys
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
                            training_file: str,
                            cpu: int,
                            bsr: float,
                            translation_table: int,
                            no_cleanup:bool) -> None:
    """
    Creates a schema structure based on the recommendations provided in the recommendations file.

    Parameters
    ----------
    recommendations_file : str
        Path to the file containing the recommendations.
    fastas_folder : str
        Path to the folder containing the FASTA files.
    output_directory: str
        Path to the directory where the final schema will be stored.
    training_file: str
        Path to the Prodigal training file that will be included in the directory of the adapted schema.
    cpu : int
        Number of CPU cores that will be used to run the process.
    bsr : float
        The BLAST Score Ratio value that will be used to evaluate
        allele similarity and select representative alleles.
    translation_table : int
        Genetic code used to translate alleles.
    no_cleanup : bool
        Flag to indicate whether to clean up temporary files.

    Returns
    -------
    None
        The function writes the output files to the specified directory.
    """

    output_d= os.path.abspath(output_directory)

    # Get all FASTA paths in the FASTA folder
    fastas_files: Dict[str, str] = {
        os.path.basename(fasta_file).split('.')[0]: os.path.join(fastas_folder, fasta_file)
        for fasta_file in os.listdir(fastas_folder) if fasta_file.endswith('.fasta')
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
            if len(line.split('\t')) == 3:
                locus, recommendation, class_id = line.split('\t', 2)
            else:
                locus, recommendation, class_id, other = line.split('\t', 3)
            # If there still is a action Choice in the recommendation file
            # If yes, then exit the module, that action is not accepted
            if recommendation == "Choice":
                pf.print_message('The input recommendation file still has loci labeled "Choice".', 'warning')
                pf.print_message('Please change these into Add, Join or Drop.', 'warning')
                sys.exit()
            # Check if the actions make sense for the class
            if class_id == '6' and recommendation == "Join":
                pf.print_message('A locus with classification 6 can not have the recommendation Join.', 'warning')
                pf.print_message(f'Please change the recommendation of the locus {locus}.', 'warning')
                sys.exit()
            # Check if the recommendation is different from the preivous one
            # If so, start a new set of IDs
            if recommendation != last_rec:
                # Split the IDs into a list
                ids_list = []
                last_rec = recommendation
            # Save the action and the IDs in the action_list dictionary
            ids_list.append(locus)
            if recommendation in action_list.get(action_id, {}) and len(ids_list)==1 :
                action_id += 1
                action_list.setdefault(action_id, {}).update({recommendation: ids_list})
            else:
                action_list.setdefault(action_id, {}).update({recommendation: ids_list})

    processed_files: List[str] = []
    all_ids: List[str] = []
    total_join: int = 0
    total_groups: int = 0
    total_drop: int = 0
    total_add: int = 0
    # For each action in the action_list dictionary
    # Recomendations can be 'Join', 'Choice', 'Drop' or 'Add'
    for action_id, recommendations in action_list.items():
        # For each recommendation in the action dictionary
        pf.print_message(f'Processing group {action_id}', 'info')
        for recommendation, ids_list in recommendations.items():
            pf.print_message(f'Involves the following loci: {ids_list} with action {recommendation}')
            # If the recommendation is 'Join'
            if "Join" in recommendation:
                if len(ids_list) < 2:
                    pf.print_message('There should be at least 2 locus in a row with the Join recommendation.', 'warning')
                    pf.print_message(f'Reorder the cluster or change the recommendation of the locus {locus}.', 'warning')
                    sys.exit()
                pf.print_message(f'These loci will be joined under the locus name {ids_list[0]}')
                output_file: str = os.path.join(temp_fasta_folder, f'{ids_list[0]}.fasta')
                # Append the new FASTA file path to the new_fastas_path list
                new_fastas_path.append(output_file)
                # Write the new FASTA file with the desired outcome
                allele_id: int = 1  # Initialize the allele_id
                total_alleles: int = 0
                with open(output_file, 'w') as out:
                    seen_fastas: List[str] = []  # Initialize the seen_fastas list that stores FASTA hashes
                    # For each ID in the ids_list
                    for id_ in ids_list:
                        if id_ in fastas_files:  # If the ID is in the fastas_files dictionary
                            processed_files.append(id_) # Add the ID to the processed_files list
                            fasta_file: str = fastas_files[id_]  # Get the FASTA file path
                            fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_file, out)  # Fetch the FASTA dictionary
                            total_join +=1 
                            all_ids.append(id_)
                            # For each header and sequence in the FASTA dictionary
                            for header, seq in fasta_dict.items():
                                total_alleles += 1
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
                total_groups += 1
                pf.print_message(f'In this group there were a total of {total_alleles} alleles, out of which {allele_id-1} were unique.', 'info')
            elif "Add" in recommendation:
                for id_ in ids_list:
                    processed_files.append(id_) # Add the ID to the processed_files list
                    output_file = os.path.join(temp_fasta_folder, f'{id_}.fasta')
                    # Append the new FASTA file path to the new_fastas_path list
                    new_fastas_path.append(output_file)
                    fasta_file= fastas_files[id_]  # Get the FASTA file path
                    shutil.copy(fasta_file, output_file)
                    pf.print_message(f'File {id_} copied to {output_file}', "info")
                    total_add += 1
                    all_ids.append(id_)
            # If the recommendation is 'Drop'
            elif "Drop" in recommendation:
                processed_files.extend(ids_list) # Add the IDS to the processed_files list
                total_drop += len(ids_list)
                for id_ in ids_list:
                    all_ids.append(id_)
                pf.print_message(f"The following IDs: {', '.join(ids_list)} have been removed due to drop action", "info")
            else:
                pf.print_message(f'The action of ids {ids_list} is not recognized. Chose beteen Add, Join and Drop.', 'warning')
                sys.exit()


    # Create schema structure
    pf.print_message("Create Schema Structure...", "info")
    # Schema path
    schema_folder = os.path.join(output_d, 'schema')
    schema_path = os.path.join(schema_folder, 'adapted_schema')
    AdaptLoci.adapt_loci(temp_fasta_folder, schema_folder, training_file, cpu, bsr, translation_table)

    extra_schema_ids: List[str] = []
    for id_ in fastas_files:
        if id_ not in all_ids:
            extra_schema_ids.append(id_)
    pf.print_message(extra_schema_ids)

    # Print final statistics
    pf.print_message('')
    pf.print_message(f'The input schema ({fastas_folder}) has {len(fastas_files)} loci.', 'info')
    pf.print_message(f'\t{total_join} loci were joined into {total_groups} groups.', 'info')
    pf.print_message(f'\t{total_drop} loci were dropped from the final schema.', 'info')
    pf.print_message(f'\t{total_add} loci were directly added into the final schema.', 'info')
    if len(extra_schema_ids) != 0:
        pf.print_message(f'The original schema has loci that were not present in the recommendation file and were then not moved into the finl schema.', 'info')
        pf.print_message(f'These loci were:\n' + "\n".join(extra_schema_ids), 'info')


    final_schema: List[str] = []
    final_schema += [file for file in os.listdir(schema_path) if file.endswith('.fasta')]
    pf.print_message(f'The final schema has a total of {len(final_schema)} loci.', 'info')
    pf.print_message('')

    if not no_cleanup:
        pf.print_message("\nCleaning up temporary files...", "info")
        ff.cleanup(output_d, [schema_folder, logf.get_log_file_path(gb.LOGGER)])