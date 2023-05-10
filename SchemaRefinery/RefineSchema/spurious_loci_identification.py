import os
import subprocess

LOCUS_CLASSIFICATIONS_TO_CHECK = ["ASM", 
                                  #"ALM", 
                                #   "NIPH", 
                                #   "NIPHEM"
                                  ]

def check_and_make_directory(dir:str):
    if not os.path.isdir(dir):
        os.mkdir(dir)

def run_blast_for_CDS(locus, classification, subjects_number, blast_results_dir, representative_file, current_cds_file):
    blast_results_file = os.path.join(blast_results_dir, f"blast_results_{locus}_{classification}.txt")
    blast_args = ['blastn', '-query', representative_file, '-subject', current_cds_file, '-out', blast_results_file]

    print(f"Running BLAST for locus: {locus}")
    print(f"Classification: {classification}.")
    print(f"Number of subjects - {subjects_number}")

    blast_proc = subprocess.Popen(blast_args,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()
    if len(stderr) > 0:
        print(stderr)

def run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory):
    blast_results_all_representatives = os.path.join(output_directory, "blast_results_all_representatives")
    check_and_make_directory(blast_results_all_representatives)

    print("Running Blast for all representatives...")

    for locus in loci:
        blast_results_file = os.path.join(blast_results_all_representatives, f"blast_results_all_representatives_{locus}.txt")
        blast_args = ['blastn', '-query', representative_file_dict[locus], '-subject', all_representatives_file, '-out', blast_results_file]

        print(f"Running BLAST for locus: {locus}")

        blast_proc = subprocess.Popen(blast_args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

        stderr = blast_proc.stderr.readlines()
        if len(stderr) > 0:
            print(stderr)

def main(schema, missing_classes_fasta, output_directory):
    check_and_make_directory(output_directory)

    blast_results_dir = os.path.join(output_directory, "blast_results")
    check_and_make_directory(blast_results_dir)

    representatives_dir = os.path.join(output_directory, "representatives")
    check_and_make_directory(representatives_dir)

    schema_files_paths = {f.replace(".fasta", ""): os.path.join(schema, f) for f in os.listdir(schema) if not os.path.isdir(f) and ".fasta" in f}
    loci = schema_files_paths.keys()

    all_representatives_file = os.path.join(representatives_dir, f"All_representatives.fasta")
    representative_file_dict = {}

    # iterate through all representatives
    with open(all_representatives_file, "w") as all_reps_file:
        for locus in loci:
            representative_file = os.path.join(representatives_dir, f"{locus}_representative.fasta")
            representative_file_dict[locus] = representative_file

            with open(schema_files_paths[locus], "r") as locus_file:
                locus_file_lines = locus_file.readlines()
                with open(representative_file, "w") as rep_file:
                    rep_file.writelines([locus_file_lines[0], locus_file_lines[1]])

                all_reps_file.writelines([locus_file_lines[0], locus_file_lines[1]])

    # only received the schema
    if schema and not missing_classes_fasta:
        # Run BLAST for all representatives
        run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory)

    # received both arguments
    if schema and missing_classes_fasta:
        for locus in loci:
            for classification in LOCUS_CLASSIFICATIONS_TO_CHECK:
                subjects_number = 0
                with open(missing_classes_fasta, "r") as missing_classes_fasta_file:
                    lines = missing_classes_fasta_file.readlines()
                    current_line_idx = 0
                    
                    while current_line_idx < len(lines):
                        _, _, info, cds_name = lines[current_line_idx].split("|")
                        if locus in info and classification in info:
                            current_cds_file = os.path.join(output_directory, f"CDS_{locus}_{classification}.fasta")
                            with open(current_cds_file, "a") as out_file:
                                cds = lines[current_line_idx+1] # sequence
                                out_file.writelines([f">{cds_name.split('&')[0]}\n", f"{cds}"])
                            subjects_number += 1

                        # jumping 2 lines of the fasta file to skip the sequence
                        current_line_idx+=2

                # only do blast if any subjects are found
                if subjects_number > 0:
                    run_blast_for_CDS(locus, classification, subjects_number, blast_results_dir, representative_file_dict[locus], current_cds_file)

        # Run BLAST for all representatives
        run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory)
