import os
from typing import Dict, List, Set, Tuple, Union
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    from utils import (
        blast_functions as bf,
        sequence_functions as sf,
        file_functions as ff,
        linux_functions as lf,
        alignments_functions as af,
        clustering_functions as cf,
        iterable_functions as itf,
        pandas_functions as pf,
        Types as tp,
        print_functions as prf
    )
except ModuleNotFoundError:
    from SchemaRefinery.utils import (
        blast_functions as bf,
        sequence_functions as sf,
        file_functions as ff,
        linux_functions as lf,
        alignments_functions as af,
        clustering_functions as cf,
        iterable_functions as itf,
        pandas_functions as pf,
        Types as tp,
        print_functions as prf
    )

def get_protein_annotation_fasta(seqRecord: SeqRecord, genbank_table_columns: List[str]) -> Dict[str, List[str]]:
    """
    Get the translated protein from a Genbank file.

    Parameters
    ----------
    seqRecord : SeqRecord
        BioPython sequence record object.

    Returns
    -------
    Dict[str, List[str]]
        Dict containing the translated protein as key and the values are
        list containing the product, the gene name, and the translated protein.

    Notes
    -----
    Source: https://github.com/LeeBergstrand/Genbank-Downloaders/blob/d904c92788696b02d9521802ebf1fc999a600e1b/SeqExtract.py#L48
    """
    cds_info: Dict[str, List[str]] = {}
    features = seqRecord.features  # Each sequence has a list (called features) that stores seqFeature objects.
    for feature in features:  # For each feature on the sequence
        if feature.type == "CDS":  # CDS means coding sequence (These are the only features we're interested in)
            featQualifiers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
            # data about the sequence feature (for example the translation)

            # Gets the required qualifier. Uses featQualifers.get to return the qualifier or a default value if the quantifier
            # is not found. Calls strip to remove unwanted brackets and ' from qualifier before storing it as a string.
            protein_id: str = str(featQualifiers.get('protein_id', '')).strip('\'[]')

            if protein_id == 'no_protein_id' or cds_info.get(protein_id):
                continue  # Skips the iteration if protein has no id or already exists in the dictionary (duplicate)
            # Get all of the relevant values from the genbank file CDS
            for column_name in genbank_table_columns:
                cds_info.setdefault(protein_id, []).append(str(featQualifiers.get(column_name, '')).strip('\'[]'))

    return cds_info

def genbank_annotations(genbank_files: str, schema_directory: str,
                        output_directory: str, cpu: int,
                        bsr: float, translation_table: int,
                        clustering_sim: float, clustering_cov: float,
                        size_ratio: float, run_mode: str,
                        extra_genbank_table_columns: List[str],
                        genbank_ids_to_add: List[str]) -> str:
    """
    Process GenBank files to extract annotations and perform BLAST searches.

    Parameters
    ----------
    genbank_files : str
        Directory containing GenBank files.
    schema_directory : str
        Directory containing schema files.
    output_directory : str
        Directory to save output files.
    cpu_cores : int
        Number of CPU cores to use for BLAST.
    bsr : float
        BLAST score ratio threshold.

    Returns
    -------
    str
        Path to the annotations file.
    """
    prf.print_message("Loading GenBank files...", "info")
    # List and sort GenBank files
    gbk_files: List[str] = [os.path.join(genbank_files, f) for f in os.listdir(genbank_files)]
    gbk_files.sort()

    prf.print_message("Extracting protein annotations from GenBank files...", "info")
    # Initialize variables
    i: int = 0
    all_cds_info: Dict[str, List[str]] = {}
    all_translation_dict: Dict[str, str] = {}
    processed_proteins: Set[str] = set()
    total_proteins: int = 0
    all_genbank_files_ids: Dict[str, List[str]] = {}
    # Create dictionaries to hashed sequence and its representative ID
    hash_to_rep_id: Dict[str, str] = {}
    # Create dictionaries to store the representatives ID and what it representes
    same_protein_other_annotations: Dict[str, List[Tuple[str, ...]]] = {}

    genbank_table_columns: List[str] = ['gene', 'product', 'translation'] + extra_genbank_table_columns
    # Parse GenBank files and extract protein annotations
    for f in gbk_files:
        file_name: str = os.path.basename(f)
        prf.print_message(f"Extracting protein annotations from GenBank file: {file_name}", "info", end='\r', flush=True)
        recs: List[SeqRecord] = [rec for rec in SeqIO.parse(f, 'genbank')]
        for r in recs:
            cds_info: Dict[str, List[str]] = get_protein_annotation_fasta(r, genbank_table_columns)
            total_proteins += len(cds_info)
            for id_, info in cds_info.items():
                # Skip if protein has no id and no sequence
                if not id_:
                    continue
                # Add all IDs in the Genbank file
                all_genbank_files_ids.setdefault(file_name.removesuffix('.gbff'), []).append(id_)
                translated_protein: str = info.pop(2)
                protein_hash: str = sf.hash_sequence(translated_protein)
                if protein_hash in processed_proteins:
                    # If the protein has already been processed, save the other annotations that it may have
                    same_protein_other_annotations.setdefault(hash_to_rep_id[protein_hash], []).append((id_, *info))
                    continue
                else:
                    processed_proteins.add(protein_hash)
                    all_translation_dict[id_] = translated_protein
                    hash_to_rep_id[protein_hash] = id_

            all_cds_info.update(cds_info)

    prf.print_message(f"Out of {total_proteins} proteins sequences {len(all_translation_dict)} are unique proteins", "info")

    prf.print_message("Clustering protein sequences...", "info")
    all_alleles: Dict[str, List[str]] = {}
    reps_sequences: Dict[str, str] = {}
    reps_groups: Dict[str, List[str]] = {}
    prot_len_dict: Dict[str, int]

    # Cluster protein sequences
    [all_alleles,
     reps_sequences,
     reps_groups,
     prot_len_dict] = cf.minimizer_clustering(
                                            all_translation_dict,
                                            5,
                                            5,
                                            True,
                                            1,
                                            all_alleles,
                                            reps_sequences,
                                            reps_groups,
                                            1,
                                            clustering_sim,
                                            clustering_cov,
                                            True,
                                            size_ratio
    )
    # Print the number of clusters
    prf.print_message(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} clusters.", "info")

    # Save extracted protein sequences to a file
    blast_processing_folder: str = os.path.join(output_directory, 'blast_processing')
    ff.create_directory(blast_processing_folder)
    genbank_protein_file: str = os.path.join(blast_processing_folder, 'selected_genbank_proteins.fasta')
    with open(genbank_protein_file, 'w') as outfile:
        for protein_id, values in reps_sequences.items():
            outfile.write(f">{protein_id}\n{values}\n")

    translation_dict: Dict[str, str]
    reps_ids: Dict[str, List[str]]
    translations_paths: Dict[str, str]
    translation_dict, reps_ids, translations_paths = sf.translate_schema_loci(
        schema_directory,
        output_directory,
        translation_table,
        run_mode
    )

    # Create BLASTdb
    blastdb_path: str = os.path.join(blast_processing_folder, 'blastdb')
    ff.create_directory(blastdb_path)
    blast_db_files: str = os.path.join(blastdb_path, 'genbank_protein_db')
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    bf.make_blast_db(makeblastdb_exec, genbank_protein_file, blast_db_files, 'prot')

    max_id_length: int = len(max(reps_ids))
    # Get Path to the blastp executable
    get_blastp_exec: str = lf.get_tool_path('blastp')
    # Calculate self-score
    self_score_dict: Dict[str, float] = bf.calculate_self_score(
        translations_paths,
        get_blastp_exec,
        blast_processing_folder,
        max_id_length,
        cpu
    )

    prf.print_message("Running BLASTp...", "info")
    blastp_results_folder: str = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, List[Union[str, float]]] = {}
    best_bsr_values_per_genbank_file: Dict[str, Dict[str, List[Union[str, float]]]] = {k: {} for k in all_genbank_files_ids.keys()}
    total_blasts: int = len(reps_ids)
    # Run BLASTp in parallel
    blastp_results_files = bf.run_blastp_operations(cpu,
                                                    get_blastp_exec,
                                                    blast_db_files,
                                                    translations_paths,
                                                    blastp_results_folder,
                                                    total_blasts,
                                                    max_id_length)

    for blast_result_file in blastp_results_files:
        # Get the alignments
        filtered_alignments_dict: tp.BlastDict 
        filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(blast_result_file, 0, True, False, True, True, False)

        # Since BLAST may find several local alignments, choose the largest one to calculate BSR.
        for query, subjects_dict in filtered_alignments_dict.items():
            # Get the loci name
            loci: str = query.split('_')[0]
            # Create the dict of the query
            bsr_values.setdefault(query, {})
            for subject_id, results in subjects_dict.items():
                if '|' in subject_id:
                    subject_id = subject_id.split('|')[1]  # Get the original ID and not the modified Blast version
                # Highest score (First one)
                subject_score: float = next(iter(results.values()))['score']
                # Calculate BSR value
                bsr_value: float = bf.compute_bsr(subject_score, self_score_dict[query])
                # Check if the BSR value is higher than the threshold
                if bsr_value >= bsr:
                    # Round BSR values if they are superior to 1.0 to 1 decimal place
                    if bsr_value > 1.0:
                        bsr_value = round(bsr_value, 1)
                    # Save all of the different matches that this query had and their BSR values
                    bsr_values[query].update({subject_id: bsr_value})
                else:
                    continue

                # Extract extra information
                extra_info: List[str] = ['NA' if element == '' else element for element in all_cds_info[subject_id]]
                # Check if the BSR value is the best for the locus
                current_best_bsr: List[Union[str, float]] = best_bsr_values.get(loci, [])
                # If there is a previous BSR value for the locus, check if the current BSR value is higher
                # We are interested in the best match only
                if not current_best_bsr:
                    best_bsr_values[loci] = [subject_id, bsr_value, *extra_info]
                elif bsr_value > current_best_bsr[1]:
                    best_bsr_values[loci] = [subject_id, bsr_value, *extra_info]

                # Get best value for genbank file
                genbank_file: str = itf.identify_string_in_dict_get_key(subject_id, all_genbank_files_ids)
                current_best_in_genbank_file: List[Union[str, float]] = best_bsr_values_per_genbank_file[genbank_file].get(loci,[])
                if current_best_in_genbank_file and bsr_value > current_best_in_genbank_file[1]:
                    best_bsr_values_per_genbank_file[genbank_file][loci] = [subject_id, bsr_value, *extra_info]
                else:
                    best_bsr_values_per_genbank_file[genbank_file][loci] = [subject_id, bsr_value, *extra_info]

    prf.print_message("\nExtracting best annotations for genbank files for proteins that were deduplicated and clustered...", "info")
    # Add all the other best matches that files may have but are being represented by other sequence
    for gbk_file, loci_values in list(best_bsr_values_per_genbank_file.items()):
        for loci_id, values in list(loci_values.items()):
            # Find all of the proteins that reps representes
            proteinid = values[0]
            bsr_value = values[1]
            same_protein = same_protein_other_annotations.get(proteinid)
            # Check if the protein is being represented by another sequence
            if same_protein:
                # For all of the elements that the representative represents add them to the dict
                for id_, *info in same_protein:
                    # Get genbank file for that ID
                    genbank_file = itf.identify_string_in_dict_get_key(id_, all_genbank_files_ids)
                    # If the genbank file is the same as the one that the representative is in, skip (may be copy protein)
                    if gbk_file == genbank_file:
                        continue
                    # Add 'NA' if element is empty
                    extra_info = ['NA' if element == '' else element for element in all_cds_info[id_]]
                    # Verify if genbank file is in the dict
                    best_bsr_values_per_genbank_file[genbank_file].setdefault(loci_id, [id_, bsr_value, extra_info])
            # For clustered elements
            rep_cluster = all_alleles.get(proteinid)
            if rep_cluster:
                for values in rep_cluster:
                    id_ = values[0]
                    # Get genbank file for that ID
                    genbank_file = itf.identify_string_in_dict_get_key(id_, all_genbank_files_ids)
                    # If the genbank file is the same as the one that the representative is in, skip (may be paralogous protein)
                    if gbk_file == genbank_file:
                        continue
                    # Add 'NA' if element is empty
                    extra_info = ['NA' if element == '' else element for element in all_cds_info[id_]]
                    # Verify if genbank file is in the dict
                    best_bsr_values_per_genbank_file[genbank_file].setdefault(loci_id, [id_, bsr_value, extra_info])

    merge_files: List[str] = []
    # Save annotations
    tab = '\t'
    header: str = f"Locus\tGenbank_ID\tGenbank_product\tGenbank_name\tGenbank_BSR{tab if extra_genbank_table_columns else ''}{tab.join(extra_genbank_table_columns)}\n"
    # Create the best annotations file
    best_annotations_all_genbank_files: str = os.path.join(output_directory, "best_annotations_all_genbank_files")
    ff.create_directory(best_annotations_all_genbank_files)
    best_annotations_file: str = os.path.join(best_annotations_all_genbank_files, 'best_genbank_annotations.tsv')
    # Add the best annotations file to the list of files to merge
    merge_files.append(best_annotations_file)
    # Get the loci that did not match or the BSR value was lower than the threshold
    not_matched_or_bsr_failed_loci: set = set(translations_paths.keys()) - set(best_bsr_values.keys())
    # Write the annotations to the file
    with open(best_annotations_file, 'w') as at:
        at.write(header)
        for loci, subject_info in best_bsr_values.items():
            subject_id = subject_info[0]  # Get the original ID and not the modified Blast version
            bsr_value = subject_info[1]  # Get the BSR value
            # Check if any element is empty
            update_cds_info: List[str] = ['NA' if element == '' else element for element in all_cds_info[subject_info[0]]]
            # Write the annotations to the file
            if len(subject_info[0]) == 2: # If no extra columns are present
                at.write(f"{loci}\t{subject_id}\t{update_cds_info[:2]}\t{bsr_value}\n")
            else: # If extra columns are present
                at.write(f"{loci}\t{subject_id}\t{tab.join(update_cds_info[:2])}\t{bsr_value}\t{tab.join(update_cds_info[2:])}\n")
    
        # Write the loci that did not match or the BSR value was lower than the threshold
        for loci in not_matched_or_bsr_failed_loci:
            at.write(f"{loci}{(tab+'NA')*(len(genbank_table_columns) + 1)}\n")

    best_annotations_per_genbank_file: str = os.path.join(output_directory, "best_annotations_per_genbank_file")
    ff.create_directory(best_annotations_per_genbank_file)
    for file, loci_results in best_bsr_values_per_genbank_file.items():
        annotations_file_genbank: str = os.path.join(best_annotations_per_genbank_file, f"{file}.tsv")
        if file in genbank_ids_to_add:
            merge_files.append(annotations_file_genbank)
        with open(annotations_file_genbank, 'w') as at:
            at.write(header)
            for loci, subject_info in loci_results.items():
                subject_id: str = subject_info[0]  # Get the original ID and not the modified Blast version
                bsr_value: float = subject_info[1]  # Get the BSR value
                # Check if any element is empty and replace with NA
                update_cds_info = ['NA' if element == '' else element for element in all_cds_info[subject_id]]
                # Write the annotations to the file
                if len(subject_info[0]) == 2: # If no extra columns are present
                    at.write(f"{loci}\t{subject_id}\t{update_cds_info[:2]}\t{bsr_value}\n")
                else: # If extra columns are present
                    at.write(f"{loci}\t{subject_id}\t{tab.join(update_cds_info[:2])}\t{bsr_value}\t{tab.join(update_cds_info[2:])}\n")
    
    # Merge all annotations files that user wants
    annotations_file: str = os.path.join(output_directory, 'best_genbank_annotations.tsv')
    pf.merge_files_into_same_file_by_key(merge_files, 'Locus', annotations_file)

    return annotations_file
            