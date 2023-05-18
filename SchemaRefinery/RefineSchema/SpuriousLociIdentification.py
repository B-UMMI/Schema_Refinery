import os
import copy
import subprocess
from Bio.Seq import Seq

LOCUS_CLASSIFICATIONS_TO_CHECK = ["ASM", 
                                  #"ALM", 
                                  #"NIPH", 
                                  #"NIPHEM"
                                  ]

DNA_BASES = 'AGCT'

def check_and_make_directory(dir:str):
    if not os.path.isdir(dir):
        os.mkdir(dir)

def check_str_alphabet(input_string, alphabet):
    alphabet_chars = set(alphabet)
    string_chars = set(input_string)

    diff = string_chars - alphabet_chars

    return len(diff) == 0

def check_str_multiple(input_string, number):
    return (len(input_string) % number) == 0

def reverse_str(input_string):
    revstr = input_string[::-1]

    return revstr

def reverse_complement(input_string, alphabet):
    translation_table = str.maketrans(alphabet, alphabet[::-1])

    upper_string = input_string.upper()
    complement_string = upper_string.translate(translation_table)
    reverse_complement_string = reverse_str(complement_string)

    return reverse_complement_string

def translate_sequence(dna_str, table_id):
    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq

def translate_dna_aux(dna_sequence, method, table_id):
    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq, DNA_BASES)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq, DNA_BASES)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]

def translate_dna(dna_sequence, table_id, min_len):
    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the sequence has ambiguous bases
    valid_dna = check_str_alphabet(original_seq, DNA_BASES)
    if valid_dna is not True:
        return 'ambiguous or invalid characters'

    # check if sequence size is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return 'sequence length is not a multiple of 3'

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        translated_seq = translate_dna_aux(original_seq, translating_methods[i], table_id)
        if not isinstance(translated_seq, list):
            exception_collector.append('{0}({1})'.format(strands[i],
                                                         translated_seq.args[0]))

        i += 1
        if i == len(strands) or isinstance(translated_seq, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(translated_seq, list):
        return [translated_seq, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str
    
def run_blast(blast_args):
    blast_proc = subprocess.Popen(blast_args,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()
    if len(stderr) > 0:
        print(stderr)


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

def get_alignments(blast_results_file):
    with open(blast_results_file, "r") as f:
        alignments = [line[2:-3] for line in f.readlines() if line[0] == ">"]
    
    return alignments

def run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory, schema):
    blast_results_all_representatives = os.path.join(output_directory, "blast_results_all_representatives")
    check_and_make_directory(blast_results_all_representatives)

    alleles_protein_dir = os.path.join(output_directory, "alleles_protein")
    check_and_make_directory(alleles_protein_dir)

    blast_results_alignments = os.path.join(blast_results_all_representatives, "alignments")
    check_and_make_directory(blast_results_alignments)

    print("Running Blast for all representatives...")

    total_loci = len(loci)

    for idx, locus in enumerate(loci, 1):
        blast_results_file = os.path.join(blast_results_all_representatives, f"blast_results_all_representatives_{locus}.txt")
        blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', all_representatives_file, '-out', blast_results_file]

        print(f"Running BLAST for locus: {locus} - {idx}/{total_loci}")
        run_blast(blast_args)

        # check results file for alignments
        alignments = get_alignments(blast_results_file)
        if locus in alignments: alignments.remove(locus) # remove own locus from alignment hits

        schema_files = {f.replace(".fasta", ""): f for f in os.listdir(schema) if ".fasta" in f}

        total_alignments = len(alignments)
        for i, alignment in enumerate(alignments, 1):
            alignment_file_path = os.path.join(schema, schema_files[alignment])
            alleles_protein_file_path = os.path.join(alleles_protein_dir, f"protein_translation_{alignment}")
            with open(alignment_file_path, "r") as alleles_file:
                lines = alleles_file.readlines()
                with open(alleles_protein_file_path, "w") as alleles_protein_file:
                    for line in lines:
                        protein_translation = translate_dna(line.replace('\n', ''), "Standard", 0)
                        # the protein translation was succesful
                        if isinstance(protein_translation, list):
                            if line[0] == '>':
                                alleles_protein_file.writelines([line])
                            else:
                                alleles_protein_file.writelines([f"{str(protein_translation[0][0])}\n"])
                        else: # protein translation was not succesful
                            # TODO: something to handle dna sequences that couldn't be translated to protein
                            pass
            
            # run blast for alignment locus alleles
            allele_blast_results_file = os.path.join(blast_results_alignments, f"blast_results_alignment_{locus}_-_{alignment}.txt")
            blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', alleles_protein_file_path, '-out', allele_blast_results_file]

            print(f"\tRunning BLAST for alignment: {alignment} - {i}/{total_alignments}")
            run_blast(blast_args)

def main(schema, output_directory, missing_classes_fasta):
    check_and_make_directory(output_directory)

    blast_results_dir = os.path.join(output_directory, "blast_results")
    check_and_make_directory(blast_results_dir)

    representatives_dir = os.path.join(output_directory, "representatives")
    check_and_make_directory(representatives_dir)

    schema_files_paths = {f.replace(".fasta", ""): os.path.join(schema, f) for f in os.listdir(schema) if not os.path.isdir(f) and ".fasta" in f}
    loci = set(schema_files_paths.keys())

    all_representatives_file = os.path.join(representatives_dir, f"All_representatives.fasta")
    representative_file_dict = {}
    not_translated_dna_sequences = {}

    # iterate through all representatives
    filtered_loci = copy.deepcopy(loci)
    with open(all_representatives_file, "w") as all_reps_file:
        for locus in loci:
            representative_file = os.path.join(representatives_dir, f"{locus}_representative.fasta")
            representative_file_dict[locus] = representative_file

            with open(schema_files_paths[locus], "r") as locus_file:
                locus_file_lines = locus_file.readlines()
                protein_translation = translate_dna(locus_file_lines[1].replace('\n', ''), "Standard", 0)

                # the protein translation was succesful
                if isinstance(protein_translation, list):
                    # print(protein_translation)
                    with open(representative_file, "w") as rep_file:
                        rep_file.writelines([locus_file_lines[0], str(protein_translation[0][0])])

                    all_reps_file.writelines([locus_file_lines[0], f"{str(protein_translation[0][0])}\n"])
                else: # protein translation was not succesful
                    # TODO: something to handle dna sequences that couldn't be translated to protein
                    del representative_file_dict[locus]
                    filtered_loci.remove(locus)
                    not_translated_dna_sequences[locus] = protein_translation

    # filtered_loci = ["GCF-022070005-protein1832"]

    # only received the schema
    if schema and not missing_classes_fasta:
        # Run BLAST for all representatives
        run_blast_for_all_representatives(filtered_loci, representative_file_dict, all_representatives_file, output_directory, schema)

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
        run_blast_for_all_representatives(filtered_loci, representative_file_dict, all_representatives_file, output_directory, schema)
