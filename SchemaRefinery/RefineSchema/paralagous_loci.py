import os
import copy
import datapane as dp
import plotly.graph_objs as go
import plotly.express.colors as graph_colors
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff, 
                       sequence_functions as sf, 
                       clustering_functions as cf, 
                       blast_functions as bf, 
                       aligments_functions as af)
    from RefineSchema import (constants as ct,
                              other as o)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff, 
                                      sequence_functions as sf, 
                                      clustering_functions as cf, 
                                      blast_functions as bf, 
                                      aligments_functions as af)
    from SchemaRefinery.RefineSchema import (constants as ct,
                                             other as o)

ALIGNMENT_COLORS = [f"rgba{(*o.hex_to_rgb(color), ct.OPACITY)}" for color in graph_colors.qualitative.Alphabet]
LOCI_COLORS = graph_colors.qualitative.Plotly[:3]

def build_graph(key: str, alignments: list):
    """
    Build a graph with a list of alignments returns the go.Figure with all the traces built into a graph.

    Parameters
    ----------
    key : str
        alignment_string with the format locus_A;locus_B.
    alignments : list
        List containing alignments, either repesentatives alignment or allele alignment.

    Returns
    -------
    go.Figure(traces) : object
        go.Figure with all of the traces build into a graph.
    """

    query, subject = key.split(";")
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

    return go.Figure(traces)

def renderGraphs(processed_representatives_dict: dict, all_allele_alignments_dict: dict, filename: str, graph_title: str, graph_dir):
    """
    Renders graphs based on input dict.

    Parameters
    ----------
    processed_representatives_dict : dict
        dict that contains all the info necessary to render the graphs.
    all_allele_alignments_dict : dict
        Dict that contains all of the results of blast of representatives vs alleles for both loci that matched during representatives 
        blast step.
    filename : str
        Output file name.
    graph_title : str
        Title displayed in the graph.
    graph_dir : str
        Graph dir path.

    Returns
    -------
    Generates the graph based on input dict in graph dir path.
    """

    all_graphs_structured = []
    processed_representatives_dict_length = len(processed_representatives_dict)

    # Render graphs to html here
    if processed_representatives_dict_length > 0:
        print(f"Rendering graphs for {graph_title}")
        # fig = make_subplots(rows=math.ceil(process_alignments_dict_length/2), cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.15)

        for (key, alignments) in processed_representatives_dict.items():
            representative_graph = build_graph(key, alignments).update_layout(
                                    showlegend=False,
                                    bargap=0.5,
                                    template='ggplot2'
                                    )
            alleles_graphs = []
            for (allele_key, allele_alignments) in all_allele_alignments_dict[key].items():
                alignment_graph = build_graph(allele_key, allele_alignments).update_layout(
                                            showlegend=False,
                                            bargap=0.5,
                                            template='ggplot2'
                                            )
                alleles_graphs.append(dp.Plot(alignment_graph))

            if len(alleles_graphs) > 0:
                all_graphs_structured.append(
                        dp.Select(
                            blocks=[
                                dp.Plot(representative_graph, 
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
            else:
                print("No allele graphs for key: ", key)

        datapane_view = dp.View(
            dp.Text(f'# {graph_title}'),
            dp.Group(
                *all_graphs_structured
            )
        )

        dp.save_report(datapane_view, os.path.join(graph_dir, f'{filename}.html'))

    else:
        print(f"Did not print a graph for {filename} - there were no graphs to process.")



def process_blast_results_for_alleles(blast_results_file, alignment_ratio_threshold, pident_threshold):
    """
    main function to process the received blast results filters the results, organizes the alignments and returns a string with 
    the information of the selected alignments for the report a dictionary with all the selected alignments to build the graphs 
    this one takes into account the selection of best and worst alignment for the alleles.

    Parameters
    ----------
    blast_results_file : str
        Path to the blast result file.
    alignment_ratio_threshold : float
        Minimum alignment threshold value.
    pident_threshold : int
        Minimum pident threshold value.

    Returns
    -------
    return : tuple
        alignment_strings : str
            Used to write inside alleles_report.
        filtered_alignments_dict : dict
            Used downstream to generate graphs.
    """

    alignments_dict = af.get_alignments_dict(blast_results_file)

    # filter alignments by pident
    alignments_dict = {key: [alignment for alignment in alignments if alignment["pident"] >= pident_threshold] 
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
    num_alignments_passing_alignment_ratio_threshold = 0

    for _, alignments in alignments_dict.items():
    
        if len(alignments) > 0:
            query = alignments[0]['query']
            subject = alignments[0]['subject']

            # filter out allele blast with itself result
            if query != subject:
                query_length = alignments[0]["query_length"]
                subject_length = alignments[0]["subject_length"]
                alignments.sort(key=lambda x : x["query_start"])
                
                query_start_stops_list = [[entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"], 
                                           entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry] for entry in alignments]
                
                alignments.sort(key=lambda x : x["subject_start"])
                subject_start_stops_list = [[entry["subject_start"], entry["subject_end"], entry["query_length"], entry["subject_length"], 
                                             entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry] for entry in alignments]

                final_query_start_stop_list, alignment_query = af.join_intervals(query_start_stops_list)
                final_subject_start_stop_list, alignment_subject = af.join_intervals(subject_start_stops_list)

                alignment_query.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                alignment_subject.sort(key=lambda x : (x["stop"] - x["start"]), reverse=True)
                bigger_query_alignment = alignment_query[0]["stop"] - alignment_query[0]["start"]
                bigger_subject_alignment = alignment_subject[0]["stop"] - alignment_subject[0]["start"]

                bigger_query_ratio = bigger_query_alignment / query_length
                bigger_subject_ratio = bigger_subject_alignment / subject_length

                if bigger_query_ratio >= alignment_ratio_threshold or bigger_subject_ratio >= alignment_ratio_threshold:
                    # go through query alignments to find best and worst
                    for idx, alignment in enumerate(alignment_query):
                        custom_scoring = alignment["custom_scoring"]
                        query_alignment_size = alignment["stop"] - alignment["start"]
                        query_ratio = query_alignment_size / query_length

                        if query_ratio >= alignment_ratio_threshold:
                            num_alignments_passing_alignment_ratio_threshold += 1
                            if custom_scoring >= best_scoring:
                                best_scoring = custom_scoring
                                best_list_of_alignments = alignment["internal_alignments"]
                                best_alignment_dict = alignment
                                best_alignment_string = final_query_start_stop_list[idx]
                            if custom_scoring <= worst_scoring:
                                worst_scoring = custom_scoring
                                worst_alignment_dict = alignment
                                worst_list_of_alignments = alignment["internal_alignments"]
                                worst_alignment_string = final_query_start_stop_list[idx]
                    
                    # go through subject alignments to find best and worst
                    for idx, alignment in enumerate(alignment_subject):
                        custom_scoring = alignment["custom_scoring"]
                        subject_alignment_size = alignment["stop"] - alignment["start"]
                        subject_ratio = subject_alignment_size / subject_length

                        if subject_ratio >= alignment_ratio_threshold:
                            num_alignments_passing_alignment_ratio_threshold += 1
                            if custom_scoring >= best_scoring:
                                best_scoring = custom_scoring
                                best_list_of_alignments = alignment["internal_alignments"]
                                best_alignment_dict = alignment
                                best_alignment_string = final_subject_start_stop_list[idx]
                            if custom_scoring <= worst_scoring:
                                worst_scoring = custom_scoring
                                worst_alignment_dict = alignment
                                worst_list_of_alignments = alignment["internal_alignments"]
                                worst_alignment_string = final_subject_start_stop_list[idx]
    
    if num_alignments_passing_alignment_ratio_threshold >= 2:
        alignment_strings = [
                            f"{best_alignment_dict['query']}\t{best_alignment_dict['subject']}\t{best_alignment_string}\t{best_alignment_dict['custom_scoring']}\n", 
                            f"{worst_alignment_dict['query']}\t{worst_alignment_dict['subject']}\t{worst_alignment_string}\t{worst_alignment_dict['custom_scoring']}\n", 
                            ]
        best_filtered_key = f"{best_alignment_dict['query']};{best_alignment_dict['subject']}"
        worst_filtered_key = f"{worst_alignment_dict['query']};{worst_alignment_dict['subject']}"
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

    if num_alignments_passing_alignment_ratio_threshold == 1:
        alignment_strings = [
                            f"{best_alignment_dict['query']}\t{best_alignment_dict['subject']}\t{best_alignment_string}\t{best_alignment_dict['custom_scoring']}\n"
                            ]
        filtered_alignments_dict = {
            f"{best_alignment_dict['query']};{best_alignment_dict['subject']}": best_list_of_alignments
        }

    return (alignment_strings, filtered_alignments_dict)

def locus_alleles_protein_translation(locus_file_path, translation_file_path):
    """
    Translated the allele sequences to protein.

    Parameters
    ----------
    locus_file_path : str
        Path to the locus in the schema
    translation_file_path : str
        Path to the file where to write the translation of the alleles.

    Returns
    -------
    successful_translation : bool
        if translation was successful or not
    Creates files in the folder alleles_protein_dir.
    """

    successful_translation = False
    with open(locus_file_path, "r") as alleles_file:
        lines = alleles_file.readlines()
        with open(translation_file_path, "w") as alleles_protein_file:
            for j in range(0, len(lines), 2):
                protein_translation = sf.translate_dna(lines[j+1].replace('\n', ''), 11, 0, cds=False)
                # the protein translation was succesful
                if isinstance(protein_translation, list):
                    alleles_protein_file.writelines([lines[j]])
                    alleles_protein_file.writelines([f"{str(protein_translation[0][0])}\n"])
                    successful_translation = True
                else: # protein translation was not succesful
                    # TODO: something to handle dna sequences that couldn't be translated to protein
                    pass
    
    return successful_translation


def run_blast_representatives_vs_alleles_multiprocessing(locus_alignment_pairs_list, allele_protein_translation_dict, 
                                                         blast_results_alignments_dir, representative_file_dict, constants_threshold):
    """
    This function, based on representatives ids, runs blast of the chosen locus representative vs all of the other locus alleles and vice versa.

    Parameters
    ----------
    locus_alignment_pairs_list : list
        Contains the locus id and alignments_pair_list.
    allele_protein_translation_dict :
        Dict that contains translation os protein sequences of all of the alleles. (updated in this function)
    blast_results_alignments_dir : str
        Path to the blast alignment directory.
    representative_file_dict : dict
        Dict that contains the path to representative file for each locus(key).
    constants_threshold : list
        List that contains two constants, alignment_ratio_threshold, pident_threshold.

    Returns
    -------
    return : list
        alignments_dict_allele : dict
            Dict used downstream to generate graphs for alleles (this is subdict that updates master dict).
        allele_alignments_string_list : list
            List used to write to allele report file.
    """

    locus = locus_alignment_pairs_list[0]
    alignments_pair_list = locus_alignment_pairs_list[1]
    alignment_ratio_threshold, pident_threshold = constants_threshold


    alignments_dict_allele = {}
    allele_alignments_string_list = []

    for alignment_pair in alignments_pair_list:
        alignment_before_underscore = alignment_pair[1].split('_')[0]
        key_to_process = f"{alignment_pair[0]};{alignment_pair[1]}"

        if not key_to_process in alignments_dict_allele:
            alignments_dict_allele[key_to_process] = {}

        print(f"Blasting {locus} representative vs {alignment_before_underscore} alleles and vice versa.")

        # Run Blast for representative A - Alleles B
        allele_blast_results_file = os.path.join(blast_results_alignments_dir, f"blast_results_alignment_{locus}_-_{alignment_before_underscore}.tsv")
        blast_args = ['blastp', '-query', representative_file_dict[locus], '-subject', allele_protein_translation_dict[alignment_before_underscore], 
                      '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident', '-out', allele_blast_results_file]
        
        bf.run_blast_with_args_only(blast_args)
        allele_alignments_string, allele_alignments_dict = process_blast_results_for_alleles(allele_blast_results_file, alignment_ratio_threshold, pident_threshold)
        
        allele_alignments_string_list.append(allele_alignments_string)
        alignments_dict_allele[key_to_process].update(allele_alignments_dict)

        # Run Blast for representative B - Alleles A (inverse)
        allele_blast_results_file = os.path.join(blast_results_alignments_dir, f"blast_results_alignment_{locus}_-_{alignment_before_underscore}.tsv")
        blast_args = ['blastp', '-query', representative_file_dict[alignment_before_underscore], '-subject', allele_protein_translation_dict[locus], 
                      '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident', '-out', allele_blast_results_file]

        bf.run_blast_with_args_only(blast_args)
        allele_alignments_string, allele_alignments_dict = process_blast_results_for_alleles(allele_blast_results_file, alignment_ratio_threshold, pident_threshold)
        
        allele_alignments_string_list.append(allele_alignments_string)
        alignments_dict_allele[key_to_process].update(allele_alignments_dict)

    return [alignments_dict_allele, allele_alignments_string_list]

def split_dict_into_clusters(clustered_loci, processed_representatives_dict):
    """
    Based on the determined clusters at cf.cluster_based_on_ids function, splits the processed_representatives_dict
    into various dicts and puts them inside a list.

    Parameters
    ----------
    clustered_loci : list
        list containing lists of the ids of the clustered loci.
    processed_representatives_dict : dict
        dict that contains all the info necessary to render the graphs.

    Returns
    -------
    results_dicts : list
        list that contains dicts with the information to render the graphs.
    """

    results_dicts = []

    for loci in clustered_loci:
        loci = list(loci)
        new_dict = {}
        for key_pair in list(processed_representatives_dict.keys()):
            if any(key.split("_")[0] in loci for key in key_pair.split(";")):
                new_dict[key_pair] = processed_representatives_dict[key_pair]
                del processed_representatives_dict[key_pair]

        results_dicts.append(new_dict)
    
    return results_dicts

def run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, output_directory, schema, info_file_path, 
                                      constants_threshold, cpu):
    """
    Main function to run the blast for representatives, alleles and render graphs.

    Parameters
    ----------
    loci : set
        Contains all of the loci ids found in the schema.
    representative_file_dict : dict
        Dict that contains the path to representative file for each locus(key).    
    all_representatives_file : str
        Path to the consolidated file of all of the represenatives sequences.    
    output_directory : str
        Path to the output dir.
    schema : str
        Path to schema dir.
    info_file_path : str
        Path to the info file, where the number of loci found and how many loci aligned with another locus is written.
    constants_threshold : list
        List that contains two constants, alignment_ratio_threshold, pident_threshold.    
    cpu : int
        Number of cpus to use during the run.

    Returns
    -------
    Various files created during the process and graphs representing the potential paralogous loci.
    """
    
    blast_results_all_representatives = os.path.join(output_directory, "blast_results_all_representatives")
    ff.create_directory(blast_results_all_representatives)

    alleles_protein_dir = os.path.join(output_directory, "alleles_protein")
    ff.create_directory(alleles_protein_dir)

    blast_results_alignments_dir = os.path.join(blast_results_all_representatives, "alignments")
    ff.create_directory(blast_results_alignments_dir)

    report_file_path = os.path.join(output_directory, "report.tsv")
    alleles_report_file_path = os.path.join(output_directory, "alleles_report.tsv")

    print("Running Blast for all representatives...")

    all_representatives_alignments_dict = {}
    all_allele_alignments_dict = {}
    allele_protein_translation_dict = {}
    total_loci = len(loci)
    representative_blast_results = []

    i=1
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing, 
                                loci, repeat('blastp'), 
                                repeat(blast_results_all_representatives), 
                                repeat(representative_file_dict), 
                                repeat(all_representatives_file)):
            
            alignment_strings, filtered_alignments_dict = af.process_blast_results(res[1], constants_threshold)
            
            representative_blast_results.append([res[0], [alignment_strings, filtered_alignments_dict]])

            print(f"Running BLAST for locus representatives: {res[0]} - {i}/{total_loci}")
            i+=1
    
    with open(report_file_path, 'w') as report_file:
        report_file.writelines(["Query\t", "Subject\t", "Query Start-End\t", "Subject Start-End\t", "Query Biggest Alignment Ratio\t", 
                            "Subject Biggest Alignment Ratio\t", "Query Length\t", "Subject Length\t", "Number of Gaps\t", 
                            "Pident - Percentage of identical matches\n"])
        
        locus_alignment_pairs_list = []
        for alignment in representative_blast_results:
            locus_alignment_pairs_list.append([alignment[0], af.remove_inverse_alignments(alignment[1][1], all_representatives_alignments_dict)])

            report_file.writelines(alignment[1][0])

    # calculate unique loci that had significant alignments            
    unique_aligment_ids = set()
    for key in all_representatives_alignments_dict.keys():
        for locus in key.split(";"):
            unique_aligment_ids.add(locus.split('_')[0])

    print(f"Total of {len(unique_aligment_ids)} loci had aligments with "
          "other loci with the chosen thresholds.")
    
    unique_ids_file_path = os.path.join(output_directory, "unique_loci.tsv")
    with open(unique_ids_file_path, 'w') as unique_ids_file:
        unique_ids_file.writelines(["Locus\n"])
        unique_ids_file.writelines([f"{id}\n" for id in unique_aligment_ids])
    
    locus_alignment_pairs_list = [[locus, alignment_pair] for [locus, alignment_pair] in locus_alignment_pairs_list if locus in unique_aligment_ids]

    schema_files = {f.replace(".fasta", ""): f for f in os.listdir(schema) if f.endswith(".fasta")}
    number_of_loci = len(unique_aligment_ids)

    for i, locus in enumerate(unique_aligment_ids, 1):

        # create files for allele protein translation
        query_locus_file_path = os.path.join(schema, schema_files[locus])
        alleles_query_locus_protein_file_path = os.path.join(alleles_protein_dir, f"protein_translation_{locus}")
        locus_translation_successful = True

        if locus not in allele_protein_translation_dict:
            print(f"Translating: {locus} {i}/{number_of_loci}")
            locus_translation_successful = locus_alleles_protein_translation(query_locus_file_path, alleles_query_locus_protein_file_path)
            allele_protein_translation_dict[locus] = alleles_query_locus_protein_file_path

        # if these entries are still not on the dictionary it means the protein translation was not successful
        # so we skip them and don't run BLASTp for them
        if not locus_translation_successful:
            print(f"\tSkipped alignment: {locus} -> Failure in protein translation.")
            continue

    with open(alleles_report_file_path, 'w') as alleles_report_file:
        alleles_report_file.writelines(["Query\t", "Subject\t","Start-End\t", "Custom Score\n"])
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
            for res in executor.map(run_blast_representatives_vs_alleles_multiprocessing, locus_alignment_pairs_list,
                                    repeat(allele_protein_translation_dict), repeat(blast_results_alignments_dir), repeat(representative_file_dict)
                                    , repeat(constants_threshold)):
                
                all_allele_alignments_dict.update(res[0])
                for string in res[1]:
                    alleles_report_file.writelines(string)

    with open(info_file_path, 'a') as info_file:
        info_file.writelines([f"There were {len(unique_aligment_ids)} different Loci that aligned with another Locus.\n\n"])

    print("Rendering graphs...")

    graph_dir = os.path.join(output_directory,'graphs')
    ff.create_directory(graph_dir)

    processed_representatives_dict = af.process_alignments_for_graphs(all_representatives_alignments_dict)

    clustered_loci = cf.cluster_based_on_ids(processed_representatives_dict)
    clustered_loci_list = []
    paralagous_path = os.path.join(output_directory,'potential_paralagous_groups.tsv')

    with open(paralagous_path,'w') as paralogous_file:
        for group in clustered_loci:
            paralogous_file.write('\t'.join(map(str, group)) + '\n')

            clustered_loci_list.append(group)

    #Render graphs based on previously clustered groups
    for i, representative_dict in enumerate(split_dict_into_clusters(clustered_loci_list, processed_representatives_dict), 1):
        renderGraphs(representative_dict, all_allele_alignments_dict, f"graphs_{i}", f"Graphs_{i}", graph_dir)

def translate_representatives(schema, output_directory, alignment_ratio_threshold, pident_threshold):
    """
    Translates the representatives and creates the necessary files and dicts for downstream usage.

    Parameters
    ----------
    schema : str
        Path to the schema directory.
    output_directory : str
        Path to the output directory.
    alignment_ratio_threshold : float
        Aligment ratio threshold for BLAST.
    pident_threshold : int
        Pident threshold for BLAST

    Returns
    -------
    returns : list
        filtered_loci : set
            Contains all of the loci ids found in the schema.
        representative_file_dict : dict
            Dict that contains the path to representative file for each locus(key). 
        all_representatives_file : str
            Path to the consolidated file of all of the represenatives sequences.
        info_file_path : str
            Path to the info file, where the number of loci found and how many loci aligned with another locus is written.
        constants_threshold : list
            List that contains two constants, alignment_ratio_threshold, pident_threshold.
    """

    info_file_path = os.path.join(output_directory, "info.txt")

    constants_threshold = [alignment_ratio_threshold, pident_threshold]
    # delete the old file for the info if it already exists, since we're appending lines to it
    ff.check_and_delete_file(info_file_path)

    # use short directory fasta files
    schema_short = os.path.join(schema, "short")

    ff.create_directory(output_directory)

    representatives_dir = os.path.join(output_directory, "representatives")
    ff.create_directory(representatives_dir)

    schema_files_paths = {f.replace("_short.fasta", ""): os.path.join(schema_short, f) for f in os.listdir(schema_short) if not os.path.isdir(f) and f.endswith(".fasta")}
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
                    protein_translation = sf.translate_dna(locus_file_lines[j+1].replace('\n', ''), "Standard", 0, cds=False)
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
        
    return [filtered_loci, representative_file_dict, all_representatives_file, info_file_path, constants_threshold]
        
def main(schema, output_directory, alignment_ratio_threshold, pident_threshold, cpu):

    output_directory = os.path.join(output_directory)
    
    loci, representative_file_dict, all_representatives_file, info_file_path, constants_threshold = translate_representatives(schema, 
                                                                                                                              output_directory, 
                                                                                                                              alignment_ratio_threshold, 
                                                                                                                              pident_threshold)

    run_blast_for_all_representatives(loci, representative_file_dict, all_representatives_file, 
                                      output_directory, schema, info_file_path, constants_threshold, cpu)
