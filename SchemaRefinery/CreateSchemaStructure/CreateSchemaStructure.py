import os

try:
    from utils import (sequence_functions as sf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf)

def create_schmea_structure(instructions_file, fastas_folder, skip_choices, output_directory):
    # Get all fastas paths in the fastas folder
    fastas_files = {os.path.basename(fasta_file).split('.')[0]: os.path.join(fastas_folder, fasta_file)
                    for fasta_file in os.listdir(fastas_folder)}
    action_list = {}
    i = 1
    with open(instructions_file, 'r') as f:
        for line in f:
            # Strip any leading/trailing whitespace characters
            line = line.strip()
            # Skip empty lines and lines starting with '#'
            if line == '#':
                i += 1
                continue
            # Split the line into the action and the ids
            recommendation, ids = line.split('\t')
            # Split the ids into a list
            ids_list = ids.split(',')
            # Save the action and the ids in the action_list dictionary
            action_list.setdefault(i, {}).update({recommendation: ids_list})
    
    # For each action in the action_list dictionary
    for action, ids in action_list.items():
        # For each recommendation in the action dictionary
        for recommendation, ids_list in ids.items():
            # If the recommendation is 'Choice' and skip_choices is True, skip the recommendation
            if skip_choices and recommendation == "Choice":
                continue
            # If the recommendation is 'Joined'.
            elif recommendation == "Joined":
                new_file_name = ids_list[0]
                output_file = os.path.join(output_directory, f'{new_file_name}.fasta')
                # Write the new fasta file with the desired outcome
                with open(output_file, 'w') as out:
                    allele_id = 1 # Initialize the allele_id
                    seen_fastas = [] # Initialize the seen_fastas list that store fasta hashes
                    # For each id in the ids_list
                    for id in ids_list:
                        if id in fastas_files: # If the id is in the fastas_files dictionary
                            fasta_file = fastas_files[id] # Get the fasta file path
                            fasta_dict = sf.fetch_fasta_dict(fasta_file, out) # Fetch the fasta dictionary
                            # For each header and sequence in the fasta dictionary
                            for header, seq in fasta_dict.items():
                                fasta_hash = sf.hash_sequence(seq) # Get the hash of the sequence
                                # If the fasta hash is not in the seen_fastas list
                                if fasta_hash not in seen_fastas:
                                    # Write the new header and sequence to the output file
                                    out.write(f'>{new_file_name}_{allele_id}\n{seq}\n')
                                    # Increment the allele_id and add the fasta hash to the seen_fastas list
                                    allele_id += 1
                                    seen_fastas.append(fasta_hash)
                                else:
                                    continue

                            print(f'File {id} added to {new_file_name} at {output_file}')
                        else:
                            print(f'File {id} not found in the fastas folder')
            else:
                print(f"The following IDS: {''.join(ids_list)} Have been removed due to drop action")