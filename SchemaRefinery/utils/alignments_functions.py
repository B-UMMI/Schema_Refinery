from copy import deepcopy
try:
    from RefineSchema.constants import MAX_GAP_UNITS
except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema.constants import MAX_GAP_UNITS

def join_intervals(alignments):
    """
    function to join alignments that intersect with each other and merge them into one single alignment
    but this one takes into account the custom scoring that needs to be calculated for the alleles.

    Parameters
    ----------
    alignments : list
        List containing the data obtained from blast.

    Returns
    -------
    return : tuple
        final_start_stop_list : list
        start_stop_list_for_processing : list
        Tuple containing final_start_stop_list and start_stop_list_for_processing

    """
    
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

def filter_out_equal_alignments(original:list, inverted:list):
    """
    Loop through a list of alignments and another list of alignments but with the key (locus1:locus2) inverted if it finds
    an interval with different start and stop, it means that its a different alignment and we add it to the original list.

    Parameters
    ----------
    original : list
        List containing alignment related to one key.
    inverted : list
        List containing alignments related to the inverse of previous key.

    Returns
    -------
    new_original : list
        new original alignment.
    """

    new_original = deepcopy(original)

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
    """
    This function tries to join alignments from an inverse of another alignment, if it exists, so they appear in the same graph.

    Parameters
    ----------
    alignments_dict : dict

    Returns
    -------
    processed_alignments_dict : dict
        Processed input
    """

    keys_set = set()
    processed_alignments_dict = {}

    print("Processing alignments for graphs...")

    for key in alignments_dict.keys():
        locus_1, locus_2 = key.split(";")
        inverted_key = f"{locus_2};{locus_1}"

        if inverted_key in alignments_dict.keys():
            if key in keys_set:
                continue
            else:
                original_alignments = alignments_dict[key]
                inverted_alignments = alignments_dict[inverted_key]

                filtered_alignments = filter_out_equal_alignments(original_alignments, inverted_alignments)

                processed_alignments_dict[key] = filtered_alignments

        else:
            processed_alignments_dict[key] = alignments_dict[key]

        keys_set.add(inverted_key)
    
    return processed_alignments_dict

def get_alignments_dict(blast_results_file):
    """
    Organize alignments with the same key "Locus_A:Locus_B" into othe same dictionary builds a dictionary where key = "Locus_A:Locus_B" 
    and value = a list of all alignments for that key correspondence.

    Parameters
    ----------
    blast_results_file : str
        Path to the blast result file.

    Returns
    -------
    alignments_dict : dict
        Dictionary containing the necessary information for graph building and representatives vs alleles blast.
    """

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

            key = f"{query};{subject}"
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

def process_blast_results(blast_results_file, constants_threshold):
    """
    Main function to process the received blast results filters the results, organizes the alignments and return 
    a string with the information of the selected alignments for the report a dictionary with all the selected 
    alignments to build the graphs.

    Parameters
    ----------
    blast_results_file : str
        Path to the blast result file.
    constants_threshold : list
        List that contains two constants, alignment_ratio_threshold, pident_threshold.

    Returns
    -------
    return : tuple
        alignment_strings : str
            Used to write inside alleles_report file
        filtered_alignments_dict : dict
            Used downstream to generate graphs. 
    """

    alignments_dict = get_alignments_dict(blast_results_file)

    alignment_ratio_threshold, pident_threshold = constants_threshold

    # filter alignments
    alignment_strings = []
    alignment_query = []
    alignment_subject = []
    # filter alignments by pident
    alignments_dict = {key: [alignment for alignment in alignments if alignment["pident"] >= pident_threshold] 
                       for key, alignments in alignments_dict.items()}
    # remove dictionary entries with zero alignments after filtering by pident
    alignments_dict = {key: alignments
                       for key, alignments in alignments_dict.items() if len(alignments) != 0}

    filtered_alignments_dict = deepcopy(alignments_dict)
    for key, alignments in alignments_dict.items():
    
        if len(alignments) > 0:
            query = alignments[0]['query']
            subject = alignments[0]['subject']
            
            # Review this code in pralagous finder, in gene fusions it cuts all
            # GCF_ named isolate from analysis!
            query_before_underscore = query.split("_")[0]
            subject_before_underscore = subject.split("_")[0]

            if query_before_underscore == subject_before_underscore:
                del filtered_alignments_dict[key]
                continue
            else:
                query_length = alignments[0]["query_length"]
                subject_length = alignments[0]["subject_length"]

                alignments.sort(key=lambda x : x["query_start"])
                query_start_stops_list = [[entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry] for entry in alignments]
                alignments.sort(key=lambda x : x["subject_start"])
                subject_start_stops_list = [[entry["subject_start"], entry["subject_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry] for entry in alignments]

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

                if query_ratio >= alignment_ratio_threshold or subject_ratio >= alignment_ratio_threshold:
                    query_start_stops = ';'.join(final_query_start_stop_list)
                    subject_start_stops = ';'.join(final_subject_start_stop_list)
                    alignment_string = f"{query}\t{subject}\t{query_start_stops}\t{subject_start_stops}\t{query_ratio}\t{subject_ratio}\t{query_length}\t{subject_length}\t{final_gaps}\t{final_pident}\n"
                    alignment_strings.append(alignment_string)
                else:
                    del filtered_alignments_dict[key]
        
    return (alignment_strings, filtered_alignments_dict)

def get_alignments_dict_from_blast_results(blast_results_file, pident_threshold,
                                           get_coords, get_self_score):
    """
    Reads BLAST results file and extracts the necessary items, based on input
    also fetches the coordinates based on query sequences and self-score contained
    inside the BLAST results file.

    Parameters
    ----------
    blast_results_file : str
        Path to the blast result file.
    pident_threshold : int
        Pident threshold to exclude BLAST results.
    get_coords : bool
        If to fetch coordinates for the BLAST match.
    get_self_score : bool
        If to get self-score from BLAST results (Note: if there are multiple 
        queries from which we can get self-score it returns the largest self-score. 
        Also to note, for self-score to be fecth the query must be also
        in the subjects database).

    Returns
    -------
    alignments_dict : dict
        Dictionary containing the results of the BLAST.
    self_score : int
        Value of self-score inside the BLAST results file.
    alignment_coords : dict
        Contains the coordinates for the query/subject pair, the coordinates are
        in reference to the query.
    """

    alignments_dict = {}
    alignment_coords_pident = {}
    alignment_coords_all = {}
    self_score = 0
    with open(blast_results_file, "r") as f:
        lines = f.readlines()
        i = 1
        for line in lines:
            # Extract the columns into the variables
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
            # Save the dict
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
            # Skip if entry matched itself and get self-score if needed
            if query == subject:
                # Largest self-score is choosen
                if float(pident) == 100 and get_self_score and int(score) > self_score:
                    self_score = int(score)
                continue
            
            if not query in alignments_dict.keys():
                alignments_dict[query] = {}
                if get_coords:
                    alignment_coords_all[query] = {}
                    alignment_coords_pident[query] = {}
            if not subject in alignments_dict[query].keys():
                # Create and save the first entry of BLAST
                alignments_dict[query][subject] = {i: value}
                if get_coords:
                    alignment_coords_all[query][subject] = {}
                    alignment_coords_pident[query][subject] = {}
                    
                    alignment_coords_all[query][subject] = {'query': [[int(query_start),int(query_end)]], 
                                                            'subject': [[int(subject_start),int(subject_end)]],}
                    # palign by pident
                    if float(pident) >= pident_threshold:
                        alignment_coords_pident[query][subject] = {'query': [[int(query_start),int(query_end)]],
                                                                   'subject': [[int(subject_start),int(subject_end)]],}
                    # To still create the dict entries for further values
                    else:
                        alignment_coords_pident[query][subject] = {'query': [],
                                                                   'subject': [],}
            else:
                # Save the other entries based on total number of entries present
                # to get the ID
                k = max(alignments_dict[query][subject].keys()) + 1
                alignments_dict[query][subject].update({k: value})
                if get_coords:
                    alignment_coords_all[query][subject]['query'].append([int(query_start),int(query_end)])
                    alignment_coords_all[query][subject]['subject'].append([int(subject_start),int(subject_end)])
                    # palign by pident
                    if float(pident) >= pident_threshold:
                        alignment_coords_pident[query][subject]['query'].append([int(query_start),int(query_end)])
                        alignment_coords_pident[query][subject]['subject'].append([int(subject_start),int(subject_end)])
            
    return alignments_dict, self_score, alignment_coords_all, alignment_coords_pident

def remove_inverse_alignments(alignments_dict, all_representatives_alignments_dict):
    """
    Since we are running the inverse for the alleles we have to filter the alignment if the inverse so that inverse of this alignment doesn't run again
    or else we will have repeated results on the alleles report.

    Parameters
    ----------
    alignments_dict : dict
        alignment dict containing alignment of the representative loci.
    all_representatives_alignments_dict : dict
        Alignment dict containg all of representative loci without the inverse of themselves.

    Returns
    -------
    alignments_pair_list : list
        list of list containing for each sublist the loci alignment pair ids.
    """

    filtered_alignments_dict = deepcopy(alignments_dict)
    for key in alignments_dict.keys():
        query, subject = key.split(";")
        inverse_key = f"{subject};{query}"
        if inverse_key in all_representatives_alignments_dict:
            del filtered_alignments_dict[key]

    all_representatives_alignments_dict.update(filtered_alignments_dict)

    alignments_pair_list = [key.split(";") for key in filtered_alignments_dict.keys()]

    return alignments_pair_list

def merge_intervals(intervals):
    """
    Merges intersecting intervals.

        Parameters
        ----------
        intervals : list
            List that contains list with coordinates of various BLAST matches

        Returns
        -------
        merged : list
            Dictionary with the result of merging intervals
            that overlapped.
    """

    merged = [deepcopy(intervals[0])]
    for current in intervals[1:]:
        previous = merged[-1]
        # current and previous intervals intersect
        if current[0] <= previous[1]:
            # determine top position
            previous[1] = max(previous[1], current[1])
        # current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged