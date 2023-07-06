import os
import copy
import math
import subprocess
from Bio.Seq import Seq
import shutil
import statistics

import plotly.graph_objs as go
import plotly.express.colors as graph_colors
from plotly.offline import plot
from plotly.subplots import make_subplots

LOCUS_CLASSIFICATIONS_TO_CHECK = ["ASM", 
                                  #"ALM", 
                                  #"NIPH", 
                                  #"NIPHEM"
                                  ]

DNA_BASES = 'AGCT'
MAX_GAP_UNITS = 4
ALIGNMENT_RATIO_THRESHOLD = 0.6
OPACITY = 0.2

def hex_to_rgb(hex_color: str) -> tuple:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) == 3:
        hex_color = hex_color * 2
    return int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)

ALIGNMENT_COLORS = [f"rgba{(*hex_to_rgb(color), OPACITY)}" for color in graph_colors.qualitative.Alphabet]
LOCI_COLORS = graph_colors.qualitative.Plotly[:3]

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

def join_intervals(alignments):
    start_stop_list_for_processing = [{"start": i[0], "stop": i[1], "joined_intervals": set()} for i in alignments]
    found_new_interval = True
    
    new_index = 0
    while found_new_interval:
        found_new_interval = False
        for i in range(new_index, len(start_stop_list_for_processing) - 1):
            first = start_stop_list_for_processing[i]
            second = start_stop_list_for_processing[i+1]

            if second["start"] - first["stop"] <= MAX_GAP_UNITS:
                if second["stop"] >= first["stop"]:
                    new_last = second["stop"]
                else:
                    new_last = first["stop"]

                if second["start"] <= first["start"]:
                    new_first = second["start"]
                else:
                    new_first = first["start"]

                new_joined_intervals = first["joined_intervals"].union(second["joined_intervals"]).union(set(((first['start'], first['stop']), (second['start'], second['stop']))))
                new_interval = {"start": new_first, "stop": new_last, "joined_intervals": new_joined_intervals}
                found_new_interval = True
                start_stop_list_for_processing = start_stop_list_for_processing[0:i] + [new_interval] + start_stop_list_for_processing[i+2:]
                new_index = i
                break
            else:
                continue

    for interval in start_stop_list_for_processing:
        interval['joined_intervals'] = list(interval['joined_intervals'])
        interval['joined_intervals'].sort(key=lambda x : x[0])
        interval['joined_intervals'] = [f"{i[0]}-{i[1]}" for i in interval['joined_intervals']]
        if len(interval['joined_intervals']) > 0:
            interval['joined_intervals'] = f"({';'.join(interval['joined_intervals'])})"
        else:
           interval['joined_intervals'] = "" 
    
    return ([f"{interval['start']}-{interval['stop']}{interval['joined_intervals']}" for interval in start_stop_list_for_processing], start_stop_list_for_processing)

def create_subplots(nplots, ncols, titles):
    """Create a Figure objcet with predefined subplots.

    Parameters
    ----------
    nplots : int
        Number of plots that will be displayed.
    ncols : int
        Number of subplots per row.
    titles : list
        Subplot titles.

    Returns
    -------
    subplots_fig : plotly.graph_objs._figure.Figure
        Figure object with predefined subplots.
    """
    # determine number of rows
    nrows = statistics.math.ceil((nplots / ncols))
    subplots_fig = make_subplots(rows=nrows, cols=ncols,
                                 subplot_titles=titles)

    return subplots_fig

def filter_alignments(original:list, inverted:list):
    new_original = copy.deepcopy(original)

    for i in inverted:
        has_alignment = False
        for o in original:
            if i["query_start"] == o["subject_start"] and i["query_end"] == o["subject_end"]:
                has_alignment = True
                break
        
        if not has_alignment:
            new_original.append(i)

    return new_original

def process_alignments_for_graphs(alignments_dict: dict):
    keys_set = set()
    processed_alignments_dict = {}

    print("Processing alignments for graphs...")

    for key in alignments_dict.keys():
        locus_1, locus_2 = key.split(":")
        inverted_key = f"{locus_2}:{locus_1}"

        if inverted_key in alignments_dict.keys():
            if key in keys_set:
                continue
            else:
                original_alignments = alignments_dict[key]
                inverted_alignments = alignments_dict[inverted_key]

                filtered_alignments = filter_alignments(original_alignments, inverted_alignments)

                processed_alignments_dict[key] = filtered_alignments

        else:
            processed_alignments_dict[key] = alignments_dict[key]
            print(f"\tCould not find key: {inverted_key} in alignments.")

        keys_set.add(inverted_key)
    
    return processed_alignments_dict

def printGraphs(alignments_dict: dict, filename:str, output_directory):
    file_template = '{0}/{1}'

    processed_alignments_dict = process_alignments_for_graphs(alignments_dict)
    process_alignments_dict_length = len(processed_alignments_dict)

    if process_alignments_dict_length > 0:
        fig = make_subplots(rows=math.ceil(process_alignments_dict_length/2), cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.15)
        row = 1
        col = 1
        for count, (key, alignments) in enumerate(processed_alignments_dict.items(), 1):
            print(f"Rendering graph for: {key} - {count}/{process_alignments_dict_length}")
            first_alignement_dicts = alignments[0]
            query = first_alignement_dicts["query"]
            subject = first_alignement_dicts["subject"]
            query_length = first_alignement_dicts["query_length"]
            subject_length = first_alignement_dicts["subject_length"]

            x = []
            y = []
            rows = [row, row]
            cols = [col, col]
            for alignment in alignments:
                x.append(alignment["query_start"])
                x.append(alignment["subject_start"])
                y.append(query)
                y.append(subject)

                x.append(alignment["query_end"])
                x.append(alignment["subject_end"])
                
                y.append(query)
                y.append(subject)

            traces = [go.Scatter(x=[0, subject_length], y=[subject, subject], mode="lines", line=dict(color=LOCI_COLORS[0]), marker=dict(color=LOCI_COLORS[0])),
                      go.Scatter(x=[0, query_length], y=[query, query], mode="lines", line=dict(color=LOCI_COLORS[2]), marker=dict(color=LOCI_COLORS[2])),
                    ]
            
            # add dashed lines for alignments
            color_index = 0
            num_traces = 0
            for c in range(len(x) - 1):
                if c%2 != 0:
                    continue
                if num_traces == 2:
                    num_traces = 0
                    color_index += 1
                if color_index == len(ALIGNMENT_COLORS) - 1:
                    color_index = 0

                if num_traces == 0:
                    fill = None
                else:
                    fill = 'tonexty'
                traces.append(go.Scatter(x=[x[c], x[c+1]], y=[y[c], y[c+1]], line=dict(width=4, color=ALIGNMENT_COLORS[color_index], dash='dash'), fill=fill, fillcolor=ALIGNMENT_COLORS[color_index]))
                rows.append(row)
                cols.append(col)
                num_traces += 1

            fig.add_traces(traces, rows=rows, cols=cols)

            if col < 2:
                col += 1
            else:
                col = 1
                row += 1

    fig_height = (statistics.math.ceil(len(fig.data)/2))*100

    fig.update_layout(title='Test',
                                  height=fig_height,
                                  showlegend=False,
                                  bargap=0.5,
                                  template='ggplot2'
                                  )

    # log transform axes for more compact scaling
    fig.update_xaxes(type='linear',
                                 title_text='Length')
    fig.update_yaxes(title_text='Names')

    subplot_fig_plotfile = file_template.format(output_directory,
                                                    f'{filename}.html')
    plot(fig,
             filename=subplot_fig_plotfile,
             auto_open=False)

def process_blast_results(blast_results_file):
    alignments_dict = {}
    with open(blast_results_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            cols = line.replace('\n', '').split("\t")
            query = cols[0]
            subject = cols[1]
            query_length = cols[2]
            subject_length = cols[3]
            query_start = cols[4]
            query_end = cols[5]
            subject_start = cols[6]
            subject_end = cols[7]
            length = cols[8]
            score = cols[9]
            gaps = cols[10]

            key = f"{query}:{subject}"
            value = {
                    "query": query,
                    "subject": subject,
                    "query_length": int(query_length),
                    "subject_length": int(subject_length),
                    "query_start": int(query_start),
                    "query_end": int(query_end),
                    "subject_start": int(subject_start),
                    "subject_end": int(subject_end),
                    "length": int(length),
                    "score": int(score),
                    "gaps": int(gaps)
                    }
            
            if not key in alignments_dict.keys():
                alignments_dict[f"{query}:{subject}"] = [value]
            else:
                alignments_dict[f"{query}:{subject}"].append(value)
    # filter alignments
    alignment_strings = []
    alignment_query = []
    alignment_subject = []
    filtered_alignments_dict = copy.deepcopy(alignments_dict)
    for key, alignments in alignments_dict.items():
        if len(alignments) > 0:
            query = alignments[0]['query']
            subject = alignments[0]['subject']
            gaps = alignments[0]['gaps']

            if query == subject:
                del filtered_alignments_dict[key]
            else:
                query_length = alignments[0]["query_length"]
                subject_length = alignments[0]["subject_length"]

                alignments.sort(key=lambda x : x["query_start"])
                query_start_stops_list = [[entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"]] for entry in alignments]
                alignments.sort(key=lambda x : x["subject_start"])
                subject_start_stops_list = [[entry["subject_start"], entry["subject_end"], entry["query_length"], entry["subject_length"]] for entry in alignments]

                final_query_start_stop_list, alignment_query = join_intervals(query_start_stops_list)
                final_subject_start_stop_list, alignment_subject = join_intervals(subject_start_stops_list)

                #filter for alignment percentage
                alignment_query.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                alignment_subject.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)

                bigger_query_alignment = alignment_query[0]["stop"] - alignment_query[0]["start"]
                bigger_subject_alignment = alignment_subject[0]["stop"] - alignment_subject[0]["start"]

                query_ratio = bigger_query_alignment / query_length
                subject_ratio = bigger_subject_alignment / subject_length

                if query_ratio >= ALIGNMENT_RATIO_THRESHOLD or subject_ratio >= ALIGNMENT_RATIO_THRESHOLD:
                    query_start_stops = ';'.join(final_query_start_stop_list)
                    subject_start_stops = ';'.join(final_subject_start_stop_list)
                    alignment_string = f"{query}\t{subject}\t{query_start_stops}\t{subject_start_stops}\t{query_ratio}\t{subject_ratio}\t{gaps}\n"
                    alignment_strings.append(alignment_string)
                else:
                    del filtered_alignments_dict[key]

    return (alignment_strings, filtered_alignments_dict)

def run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory, schema):
    blast_results_all_representatives = os.path.join(output_directory, "blast_results_all_representatives")
    check_and_make_directory(blast_results_all_representatives)

    alleles_protein_dir = os.path.join(output_directory, "alleles_protein")
    check_and_make_directory(alleles_protein_dir)

    blast_results_alignments = os.path.join(blast_results_all_representatives, "alignments")
    check_and_make_directory(blast_results_alignments)

    report_file_path = os.path.join(output_directory, "report.tsv")

    print("Running Blast for all representatives...")

    all_alignments_dict = {}
    best_allele_alignments_dict = {}
    total_loci = len(loci)
    with open(report_file_path, 'w') as report_file:
        report_file.writelines(["Query\t", "Subject\t", "Query Start-End\t", "Subject Start-End\t", "Query Biggest Alignment Ratio\t", "Subject Biggest Alignment Ratio\t", "Number of Gaps\n"])
        for idx, locus in enumerate(loci, 1):
            blast_results_file = os.path.join(blast_results_all_representatives, f"blast_results_all_representatives_{locus}.tsv")
            blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', all_representatives_file, '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps', '-out', blast_results_file]

            print(f"Running BLAST for locus: {locus} - {idx}/{total_loci}")
            run_blast(blast_args)

            # check results file for alignments
            _, alignments_dict = process_blast_results(blast_results_file)
            schema_files = {f.replace(".fasta", ""): f for f in os.listdir(schema) if ".fasta" in f}

            # report_file.writelines(alignments_string)

            all_alignments_dict.update(alignments_dict)

            total_alignments = len(alignments_dict)
            alignments_list = [key.split(":")[1].split("_")[0] for key in alignments_dict.keys()]

            for i, alignment in enumerate(alignments_list, 1):
                if not alignment in schema_files:
                    #TODO handle loci that appear in short directory but don't appear on main schema directory
                    print(f"\tSkipped alignment: {alignment}")
                    continue
                alignment_file_path = os.path.join(schema, schema_files[alignment])
                alleles_protein_file_path = os.path.join(alleles_protein_dir, f"protein_translation_{alignment}")
                with open(alignment_file_path, "r") as alleles_file:
                    lines = alleles_file.readlines()
                    with open(alleles_protein_file_path, "w") as alleles_protein_file:
                        for j in range(0, len(lines), 2):
                            protein_translation = translate_dna(lines[j+1].replace('\n', ''), "Standard", 0)
                            # the protein translation was succesful
                            if isinstance(protein_translation, list):
                                alleles_protein_file.writelines([lines[j]])
                                alleles_protein_file.writelines([f"{str(protein_translation[0][0])}\n"])
                            else: # protein translation was not succesful
                                # TODO: something to handle dna sequences that couldn't be translated to protein
                                pass
                
                # run blast for alignment locus alleles
                allele_blast_results_file = os.path.join(blast_results_alignments, f"blast_results_alignment_{locus}_-_{alignment}.tsv")
                blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', alleles_protein_file_path, '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps', '-out', allele_blast_results_file]

                print(f"\tRunning BLAST for alignment: {alignment} - {i}/{total_alignments}")
                run_blast(blast_args)

                # process the blast file for this alignment
                allele_alignments_string, allele_alignments_dict = process_blast_results(allele_blast_results_file)

                # Choose alignment with best score to make a graph of it
                best_score = 0
                best_alignment = None
                for key, l in allele_alignments_dict.items():
                    for alig in l:
                        if alig["score"] > best_score:
                            best_score = alig["score"]
                            best_alignment = {key: [alig]}

                best_allele_alignments_dict.update(best_alignment)

                # if not best_allele_alignment:
                #     best_allele_alignment = best_alignment
                # elif best_score > best_allele_alignment[best_key][0]["score"]:
                #     best_allele_alignment = best_alignment

                report_file.writelines(allele_alignments_string)

    print("Rendering all alignement graphs...")
    printGraphs(all_alignments_dict, "Scatter_plots_all_representatives", output_directory)

    print("Rendering best allele alignements graphs...")
    printGraphs(best_allele_alignments_dict, "Scatter_plots_best_allele_alignments", output_directory)
    
    # shutil.rmtree(blast_results_alignments)
    shutil.rmtree(alleles_protein_dir)
    # shutil.rmtree(blast_results_all_representatives)

def main(schema, output_directory, missing_classes_fasta, threshold):
    # use short directory fasta files
    schema_short = os.path.join(schema, "short")

    check_and_make_directory(output_directory)

    blast_results_dir = os.path.join(output_directory, "blast_results")
    check_and_make_directory(blast_results_dir)

    representatives_dir = os.path.join(output_directory, "representatives")
    check_and_make_directory(representatives_dir)

    schema_files_paths = {f.replace("_short.fasta", ""): os.path.join(schema_short, f) for f in os.listdir(schema_short) if not os.path.isdir(f) and ".fasta" in f}
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
                for j in range(0, len(locus_file_lines), 2):
                    protein_translation = translate_dna(locus_file_lines[j+1].replace('\n', ''), "Standard", 0)
                    # the protein translation was succesful
                    if isinstance(protein_translation, list):
                        with open(representative_file, "w") as rep_file:
                            rep_file.writelines([locus_file_lines[j], f"{str(protein_translation[0][0])}\n"])

                        all_reps_file.writelines([locus_file_lines[j], f"{str(protein_translation[0][0])}\n"])

                    else: # protein translation was not succesful
                        # TODO: something to handle dna sequences that couldn't be translated to protein
                        if locus in representative_file_dict.keys():
                            del representative_file_dict[locus]
                            filtered_loci.remove(locus)
                            not_translated_dna_sequences[locus] = protein_translation

    # filtered_loci = ["GCF-000006885-protein699","GCF-000007045-protein1266"]

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

    # shutil.rmtree(blast_results_dir)
    shutil.rmtree(representatives_dir)