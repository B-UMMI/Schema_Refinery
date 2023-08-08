import os
import copy
import math
import subprocess
from Bio.Seq import Seq
import shutil
import statistics
import datapane as dp

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
PIDENT_THRESHOLD = 75
OPACITY = 0.2

def hex_to_rgb(hex_color: str) -> tuple:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) == 3:
        hex_color = hex_color * 2
    return int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)

ALIGNMENT_COLORS = [f"rgba{(*hex_to_rgb(color), OPACITY)}" for color in graph_colors.qualitative.Alphabet]
LOCI_COLORS = graph_colors.qualitative.Plotly[:3]

def check_and_delete_file(file:str):
    if os.path.isfile(file):
        os.remove(file)

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
    protseq = Seq.translate(myseq_obj, table=table_id)

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
    start_stop_list_for_processing = [{"start": alignment[0], "stop": alignment[1], 
                                       "joined_intervals": set(), "pident": alignment[4], 
                                       "gaps": alignment[5], "length": alignment[6]} for alignment in alignments]
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
                new_interval = {"start": new_first, "stop": new_last, "joined_intervals": new_joined_intervals, "pident": f"{first['pident']};{second['pident']}", "gaps": f"{first['gaps']};{second['gaps']}"}
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

def join_and_filter_intervals_for_alleles(alignments):
    start_stop_list_for_processing = [{"start": alignment[0], "stop": alignment[1], 
                                       "joined_intervals": set(), "pident": alignment[4]/100, 
                                       "gaps": alignment[5], "length": alignment[6],
                                       "product_length_pident_list": [alignment[4]/100 * alignment[6]], 
                                       "length_list": [alignment[6]],
                                       "query": alignment[7],
                                       "subject": alignment[8],
                                       "internal_alignments": [alignment[9]]} for alignment in alignments]
    
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
                new_interval = {"start": new_first, 
                                "stop": new_last, 
                                "joined_intervals": new_joined_intervals, 
                                "pident": f"{first['pident']};{second['pident']}", 
                                "gaps": f"{first['gaps']};{second['gaps']}",
                                "product_length_pident_list": first["product_length_pident_list"] + second["product_length_pident_list"], 
                                "length_list": first["length_list"] + second["length_list"],
                                "query": first["query"],
                                "subject": first["subject"],
                                "internal_alignments": first["internal_alignments"] + second["internal_alignments"]}
                found_new_interval = True
                start_stop_list_for_processing = start_stop_list_for_processing[0:i] + [new_interval] + start_stop_list_for_processing[i+2:]
                new_index = i
                break
            else:
                continue

    for interval in start_stop_list_for_processing:
        interval["custom_scoring"] = sum(interval["product_length_pident_list"]) / sum(interval["length_list"])
        del interval["product_length_pident_list"]
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
            # print(f"\tCould not find key: {inverted_key} in alignments.")

        keys_set.add(inverted_key)
    
    return processed_alignments_dict

def build_graph(key: str, alignments: list):
    query, subject = key.split(":")
    first_alignement_dicts = alignments[0]
    query = first_alignement_dicts["query"]
    subject = first_alignement_dicts["subject"]
    query_length = first_alignement_dicts["query_length"]
    subject_length = first_alignement_dicts["subject_length"]

    x = []
    y = []
    for alignment in alignments:
        x.append(alignment["query_start"])
        x.append(alignment["subject_start"])
        y.append(query)
        y.append(subject)

        x.append(alignment["query_end"])
        x.append(alignment["subject_end"])
        
        y.append(query)
        y.append(subject)

    traces = [go.Scatter(x=[1, subject_length], y=[subject, subject], mode="lines", line=dict(color=LOCI_COLORS[0]), marker=dict(color=LOCI_COLORS[0])),
                go.Scatter(x=[1, query_length], y=[query, query], mode="lines", line=dict(color=LOCI_COLORS[2]), marker=dict(color=LOCI_COLORS[2])),
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
        num_traces += 1

    return traces

def printGraphs(representatives_dict: dict, alleles_dict:dict, filename:str, graph_title: str, output_directory):
    all_graphs_structured = []
    processed_representatives_dict = process_alignments_for_graphs(representatives_dict)
    processed_representatives_dict_length = len(processed_representatives_dict)

    # Render graphs to html here
    if processed_representatives_dict_length > 0:
        print(f"Rendering graphs for {graph_title}")
        # fig = make_subplots(rows=math.ceil(process_alignments_dict_length/2), cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.15)

        for (key, alignments) in processed_representatives_dict.items():
            representative_traces = build_graph(key, alignments)
            alleles_graphs = []
            for (allele_key, allele_alignments) in alleles_dict[key].items():
                alignment_traces = build_graph(allele_key, allele_alignments)
                alleles_graphs.append(dp.Plot(go.Figure(alignment_traces).update_layout(
                                            showlegend=False,
                                            bargap=0.5,
                                            template='ggplot2'
                                        )
                                    ))

            all_graphs_structured.append(
                    dp.Select(
                        blocks=[
                            dp.Plot(go.Figure(representative_traces).update_layout(
                                    showlegend=False,
                                    bargap=0.5,
                                    template='ggplot2'
                                    ), 
                                label="Representative"),
                            dp.Group(
                                *alleles_graphs,
                                columns=2,
                                label="Alleles",
                            ),
                        ],
                        type=dp.SelectType.DROPDOWN,
                    ),
                )

        datapane_view = dp.View(
            dp.Text(f'# {graph_title}'),
            dp.Group(
                *all_graphs_structured
            )
        )

        dp.save_report(datapane_view, os.path.join(output_directory, f'{filename}.html'))

    else:
        print(f"Did not print a graph for {filename} - there were no graphs to process.")

def get_alignments_dict(blast_results_file):
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
            pident = cols[11]

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
                    "gaps": int(gaps),
                    "pident": float(pident)
                    }
            
            if not key in alignments_dict.keys():
                alignments_dict[key] = [value]
            else:
                alignments_dict[key].append(value)

    return alignments_dict

def process_blast_results(blast_results_file):
    alignments_dict = get_alignments_dict(blast_results_file)

    # filter alignments
    alignment_strings = []
    alignment_query = []
    alignment_subject = []
    # filter alignments by pident
    alignments_dict = {key: [alignment for alignment in alignments if alignment["pident"] >= PIDENT_THRESHOLD] 
                       for key, alignments in alignments_dict.items()}
    # remove dictionary entries with zero alignments after filtering by pident
    alignments_dict = {key: alignments
                       for key, alignments in alignments_dict.items() if len(alignments) != 0}

    filtered_alignments_dict = copy.deepcopy(alignments_dict)
    for key, alignments in alignments_dict.items():
    
        if len(alignments) > 0:
            query = alignments[0]['query']
            subject = alignments[0]['subject']
            query_before_underscore = query.split("_")[0]
            subject_before_underscore = subject.split("_")[0]

            if query_before_underscore == subject_before_underscore:
                del filtered_alignments_dict[key]
                continue
            else:
                query_length = alignments[0]["query_length"]
                subject_length = alignments[0]["subject_length"]

                alignments.sort(key=lambda x : x["query_start"])
                query_start_stops_list = [[entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"]] for entry in alignments]
                alignments.sort(key=lambda x : x["subject_start"])
                subject_start_stops_list = [[entry["subject_start"], entry["subject_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"]] for entry in alignments]

                final_query_start_stop_list, alignment_query = join_intervals(query_start_stops_list)
                final_subject_start_stop_list, alignment_subject = join_intervals(subject_start_stops_list)
                
                final_gaps = ';'.join([str(alignment["gaps"]) for alignment in alignment_query])
                final_pident = ';'.join([str(alignment["pident"]) for alignment in alignment_query])

                alignment_query.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                alignment_subject.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                bigger_query_alignment = alignment_query[0]["stop"] - alignment_query[0]["start"]
                bigger_subject_alignment = alignment_subject[0]["stop"] - alignment_subject[0]["start"]

                query_ratio = bigger_query_alignment / query_length
                subject_ratio = bigger_subject_alignment / subject_length

                if query_ratio >= ALIGNMENT_RATIO_THRESHOLD or subject_ratio >= ALIGNMENT_RATIO_THRESHOLD:
                    query_start_stops = ';'.join(final_query_start_stop_list)
                    subject_start_stops = ';'.join(final_subject_start_stop_list)
                    alignment_string = f"{query}\t{subject}\t{query_start_stops}\t{subject_start_stops}\t{query_ratio}\t{subject_ratio}\t{query_length}\t{subject_length}\t{final_gaps}\t{final_pident}\n"
                    alignment_strings.append(alignment_string)
                else:
                    del filtered_alignments_dict[key]

    return (alignment_strings, filtered_alignments_dict)

def process_blast_results_for_alleles(blast_results_file):
    alignments_dict = get_alignments_dict(blast_results_file)

    # filter alignments by pident
    alignments_dict = {key: [alignment for alignment in alignments if alignment["pident"] >= PIDENT_THRESHOLD] 
                       for key, alignments in alignments_dict.items()}
    # remove dictionary entries with zero alignments after filtering by pident
    alignments_dict = {key: alignments
                       for key, alignments in alignments_dict.items() if len(alignments) != 0}

    filtered_alignments_dict = {}
    alignment_strings = []
    best_alignment_string = ""
    worst_alignment_string = ""
    best_alignment_dict = {}
    worst_alignment_dict = {}
    best_list_of_alignments = []
    worst_list_of_alignments = []
    best_scoring = 0
    worst_scoring = 1

    for key, alignments in alignments_dict.items():
    
        if len(alignments) > 0:
            query = alignments[0]['query']
            subject = alignments[0]['subject']

            # filter out allele blast with itself result
            if query == subject:
                del filtered_alignments_dict[key]
            else:
                query_length = alignments[0]["query_length"]
                alignments.sort(key=lambda x : x["query_start"])
                query_start_stops_list = [[entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry] for entry in alignments]

                final_query_start_stop_list, alignment_query = join_and_filter_intervals_for_alleles(query_start_stops_list)

                alignment_query.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                bigger_query_alignment = alignment_query[0]["stop"] - alignment_query[0]["start"]

                query_ratio = bigger_query_alignment / query_length

                if query_ratio >= ALIGNMENT_RATIO_THRESHOLD:
                    for idx, alignment in enumerate(alignment_query):
                        custom_scoring = alignment["custom_scoring"]
                        if custom_scoring > best_scoring:
                            best_scoring = custom_scoring
                            best_list_of_alignments = alignment["internal_alignments"]
                            best_alignment_dict = alignment
                            best_alignment_string = final_query_start_stop_list[idx]
                        if custom_scoring < worst_scoring:
                            worst_scoring = custom_scoring
                            worst_alignment_dict = alignment
                            worst_list_of_alignments = alignment["internal_alignments"]
                            worst_alignment_string = final_query_start_stop_list[idx]
    
    if len(alignments_dict) >= 2:
        alignment_strings = [
                            f"{best_alignment_dict['query']}\t{best_alignment_dict['subject']}\t{best_alignment_string}\t{best_alignment_dict['custom_scoring']}\n", 
                            f"{worst_alignment_dict['query']}\t{worst_alignment_dict['subject']}\t{worst_alignment_string}\t{worst_alignment_dict['custom_scoring']}\n", 
                            ]
        best_filtered_key = f"{best_alignment_dict['query']}:{best_alignment_dict['subject']}"
        worst_filtered_key = f"{worst_alignment_dict['query']}:{worst_alignment_dict['subject']}"
        if best_filtered_key != worst_filtered_key:
            filtered_alignments_dict = {
                best_filtered_key: best_list_of_alignments,
                worst_filtered_key: worst_list_of_alignments
            }
        else:
            # keys are equal so we merge all the alignments in the same list
            filtered_alignments_dict = {
                best_filtered_key: best_list_of_alignments + worst_list_of_alignments
            }

    if len(alignments_dict) == 1:
        alignment_strings = [
                            f"{best_alignment_dict['query']}\t{best_alignment_dict['subject']}\t{best_alignment_string}\t{best_alignment_dict['custom_scoring']}\n"
                            ]
        filtered_alignments_dict = {
            f"{best_alignment_dict['query']}:{best_alignment_dict['subject']}": best_list_of_alignments
        }

    return (alignment_strings, filtered_alignments_dict)

def locus_alleles_protein_translation(locus_file_path, translation_file_path):
    successful_translation = False
    with open(locus_file_path, "r") as alleles_file:
        lines = alleles_file.readlines()
        with open(translation_file_path, "w") as alleles_protein_file:
            for j in range(0, len(lines), 2):
                protein_translation = translate_dna(lines[j+1].replace('\n', ''), "Standard", 0)
                # the protein translation was succesful
                if isinstance(protein_translation, list):
                    alleles_protein_file.writelines([lines[j]])
                    alleles_protein_file.writelines([f"{str(protein_translation[0][0])}\n"])
                    successful_translation = True
                else: # protein translation was not succesful
                    # TODO: something to handle dna sequences that couldn't be translated to protein
                    pass
    
    return successful_translation

def run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory, schema):
    blast_results_all_representatives = os.path.join(output_directory, "blast_results_all_representatives")
    check_and_make_directory(blast_results_all_representatives)

    alleles_protein_dir = os.path.join(output_directory, "alleles_protein")
    check_and_make_directory(alleles_protein_dir)

    blast_results_alignments = os.path.join(blast_results_all_representatives, "alignments")
    check_and_make_directory(blast_results_alignments)

    report_file_path = os.path.join(output_directory, "report.tsv")
    alleles_report_file_path = os.path.join(output_directory, "alleles_report.tsv")

    print("Running Blast for all representatives...")

    all_representatives_alignments_dict = {}
    all_allele_alignments_dict = {}
    total_loci = len(loci)
    allele_protein_translation_dict = {}
    # allele_blast_pairs_processed = set()
    with open(report_file_path, 'w') as report_file:
        with open(alleles_report_file_path, 'w') as alleles_report_file:
            report_file.writelines(["Query\t", "Subject\t", "Query Start-End\t", "Subject Start-End\t", "Query Biggest Alignment Ratio\t", "Subject Biggest Alignment Ratio\t", "Query Length\t", "Subject Length\t", "Number of Gaps\t", "Pident - Percentage of identical matches\n"])
            alleles_report_file.writelines(["Query\t", "Subject\t", "Start - End\t", "Custom Score\n"])
            for idx, locus in enumerate(loci, 1):
                blast_results_file = os.path.join(blast_results_all_representatives, f"blast_results_all_representatives_{locus}.tsv")
                blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', all_representatives_file, '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident', '-out', blast_results_file]

                print(f"Running BLAST for locus: {locus} - {idx}/{total_loci}")
                run_blast(blast_args)

                # check results file for alignments
                alignments_string, alignments_dict = process_blast_results(blast_results_file)
                schema_files = {f.replace(".fasta", ""): f for f in os.listdir(schema) if ".fasta" in f}

                report_file.writelines(alignments_string)

                # since we are running the inverse for the alleles
                # we have to filter the alignment if the inverse has ran previously
                # or else we will have repeated results on the alleles report
                filtered_alignments_dict = copy.deepcopy()
                for key in alignments_dict.keys():
                    query, subject = key.split(":")
                    inverse_key = f"{subject}:{query}"
                    if inverse_key in all_representatives_alignments_dict:
                        del filtered_alignments_dict[key]

                all_representatives_alignments_dict.update(filtered_alignments_dict)

                total_alignments = len(filtered_alignments_dict)
                alignments_pair_list = [key.split(":") for key in filtered_alignments_dict.keys()]

                for i, alignment_pair in enumerate(alignments_pair_list, 1):
                    locus_for_key = alignment_pair[0]
                    alignment = alignment_pair[1]
                    alignment_before_underscore = alignment_pair[1].split('_')[0]
                    key_to_process = f"{locus_for_key}:{alignment}"
                    if not key_to_process in all_allele_alignments_dict:
                        all_allele_alignments_dict[key_to_process] = {}

                    # create files for allele protein translation
                    query_locus_file_path = os.path.join(schema, schema_files[locus])
                    alignment_file_path = os.path.join(schema, schema_files[alignment_before_underscore])
                    alleles_query_locus_protein_file_path = os.path.join(alleles_protein_dir, f"protein_translation_{locus}")
                    alleles_alignment_protein_file_path = os.path.join(alleles_protein_dir, f"protein_translation_{alignment_before_underscore}")

                    locus_translation_successful = True
                    alignment_translation_successful = True
                    if locus not in allele_protein_translation_dict:
                        locus_translation_successful = locus_alleles_protein_translation(query_locus_file_path, alleles_query_locus_protein_file_path)
                        allele_protein_translation_dict[locus] = alleles_query_locus_protein_file_path

                    if alignment_before_underscore not in allele_protein_translation_dict:
                        alignment_translation_successful = locus_alleles_protein_translation(alignment_file_path, alleles_alignment_protein_file_path)
                        allele_protein_translation_dict[alignment_before_underscore] = alleles_alignment_protein_file_path

                    # if these entries are still not on the dictionary it means the protein translation was not successful
                    # so we skip them and don't run BLASTp for them
                    if not alignment_translation_successful or not locus_translation_successful:
                        print(f"\tSkipped alignment: {alignment_before_underscore} - {i}/{total_alignments} -> Failure in protein translation.")
                        continue

                    # Run Blast for representative A - Alleles B
                    allele_blast_results_file = os.path.join(blast_results_alignments, f"blast_results_alignment_{locus}_-_{alignment_before_underscore}.tsv")
                    blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', allele_protein_translation_dict[alignment_before_underscore], '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident', '-out', allele_blast_results_file]
                    print(f"\tRunning BLAST for alignment: {alignment_before_underscore} - {i}/{total_alignments}")
                    run_blast(blast_args)
                    allele_alignments_string, allele_alignments_dict = process_blast_results_for_alleles(allele_blast_results_file)
                    
                    alleles_report_file.writelines(allele_alignments_string)
                    all_allele_alignments_dict[key_to_process].update(allele_alignments_dict)

                    # Run Blast for representative B - Alleles A (inverse)
                    allele_blast_results_file = os.path.join(blast_results_alignments, f"blast_results_alignment_{locus}_-_{alignment_before_underscore}.tsv")
                    blast_args = ['blastp', '-query', representative_file_dict[alignment_before_underscore], '-subject', allele_protein_translation_dict[locus], '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident', '-out', allele_blast_results_file]
                    print(f"\tRunning BLAST for the reverse of previous alignment.")
                    run_blast(blast_args)
                    allele_alignments_string, allele_alignments_dict = process_blast_results_for_alleles(allele_blast_results_file)
                    
                    alleles_report_file.writelines(allele_alignments_string)
                    all_allele_alignments_dict[key_to_process].update(allele_alignments_dict)

    # calculate unique loci that had significant alignments
    unique_alignent_ids = set()
    for key in all_representatives_alignments_dict.keys():
        locus_query, _ = key.split(":")
        unique_alignent_ids.add(locus_query.split('_')[0])

    unique_ids_file_path = os.path.join(output_directory, "unique_loci.tsv")
    with open(unique_ids_file_path, 'w') as unique_ids_file:
        unique_ids_file.writelines(["Locus\n"])
        unique_ids_file.writelines([f"{id}\n" for id in unique_alignent_ids])

    with open(info_file_path, 'a') as info_file:
        info_file.writelines([f"There were {len(unique_alignent_ids)} different Loci that aligned with another Locus.\n\n"])

    print("Rendering graphs...")
    printGraphs(all_representatives_alignments_dict, all_allele_alignments_dict, "Scatter_plots_all_representatives", "Scatter plots for all representatives", output_directory)
   
    # shutil.rmtree(blast_results_alignments)
    shutil.rmtree(alleles_protein_dir)
    # shutil.rmtree(blast_results_all_representatives)

def main(schema, output_directory, missing_classes_fasta):
    global info_file_path
    info_file_path = os.path.join(output_directory, "info.txt")

    # delete the old file for the info if it already exists, since we're appending lines to it
    check_and_delete_file(info_file_path)

    # use short directory fasta files
    schema_short = os.path.join(schema, "short")

    check_and_make_directory(output_directory)

    blast_results_dir = os.path.join(output_directory, "blast_results")
    check_and_make_directory(blast_results_dir)

    representatives_dir = os.path.join(output_directory, "representatives")
    check_and_make_directory(representatives_dir)

    schema_files_paths = {f.replace("_short.fasta", ""): os.path.join(schema_short, f) for f in os.listdir(schema_short) if not os.path.isdir(f) and ".fasta" in f}
    loci = set(schema_files_paths.keys())

    print(f"Found {len(loci)} in schema dir (short).")

    all_representatives_file = os.path.join(representatives_dir, f"All_representatives.fasta")
    representative_file_dict = {}
    not_translated_dna_sequences = {}

    # iterate through all representatives files to translate them
    filtered_loci = copy.deepcopy(loci)
    with open(all_representatives_file, "w") as all_reps_file:
        for locus in loci:
            representative_file = os.path.join(representatives_dir, f"{locus}_representative.fasta")
            representative_file_dict[locus] = representative_file

            with open(schema_files_paths[locus], "r") as locus_file:
                locus_file_lines = locus_file.readlines()
                found_translation = False
                for j in range(0, len(locus_file_lines), 2):
                    protein_translation = translate_dna(locus_file_lines[j+1].replace('\n', ''), "Standard", 0)
                    # the protein translation was succesful
                    if isinstance(protein_translation, list):
                        with open(representative_file, "w") as rep_file:
                            rep_file.writelines([locus_file_lines[j], f"{str(protein_translation[0][0])}\n"])

                        all_reps_file.writelines([locus_file_lines[j], f"{str(protein_translation[0][0])}\n"])
                        found_translation = True

                    else: # protein translation was not succesful
                        # TODO: something to handle dna sequences that couldn't be translated to protein
                        print(protein_translation)
                        continue
                if not found_translation:
                    del representative_file_dict[locus]
                    filtered_loci.remove(locus)
                    not_translated_dna_sequences[locus] = protein_translation

    with open(info_file_path, 'w') as info_file:
        info_file.writelines([
            f"Found {len(loci)} in schema dir (short).\n"
            ])

    # only received the schema
    if schema and not missing_classes_fasta:
        # Run BLAST for all representatives
        run_blast_for_all_representatives(filtered_loci, representative_file_dict, all_representatives_file, output_directory, schema)

    # received both arguments
    if schema and missing_classes_fasta:
        # TODO Should run code for CDS
        pass

    # shutil.rmtree(blast_results_dir)
    shutil.rmtree(representatives_dir)
