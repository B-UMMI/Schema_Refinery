import os

def create_schmea_structure(instructions_file, fastas_folder, output_directory):
    fastas_files = {fasta_file: os.path.join(fastas_folder, fasta_file) for fasta_file in os.listdir(fastas_folder)}
    with open(instructions_file) as f:
        