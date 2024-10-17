import os
import concurrent.futures
from itertools import repeat 
from typing import Dict, List, Set, Tuple, Union
from Bio import SeqIO, SeqRecord

try:
    from utils import (
        blast_functions as bf,
        sequence_functions as sf,
        file_functions as ff,
        linux_functions as lf,
        alignments_functions as af,
        clustering_functions as cf,
        iterable_functions as itf,
    )
except ModuleNotFoundError:
    from SchemaRefinery.utils import (
        blast_functions as bf,
        sequence_functions as sf,
        file_functions as ff,
        linux_functions as lf,
        alignments_functions as af,
        clustering_functions as cf,
        iterable_functions as itf
    )

def get_protein_annotation_fasta(seqRecord: SeqRecord) -> Dict[str, List[str]]:
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

            if protein_id == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.

            gene: str = str(featQualifiers.get('gene', '')).strip('\'[]')
            product: str = str(featQualifiers.get('product', 'no_product_name')).strip('\'[]')
            translated_protein: str = str(featQualifiers.get('translation', 'no_translation')).strip('\'[]')

            cds_info[protein_id] = [product, gene, translated_protein]

    return cds_info

def genbank_annotations(genbank_files: str, schema_directory: str,
                        output_directory: str, cpu: int,
                        bsr: float, translation_table: int,
                        clustering_sim: float, clustering_cov: float,
                        size_ratio: float, run_mode: str) -> str:
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
    
    # Create output directory for genbank annotations
    output_directory: str = os.path.join(output_directory, 'genbank_annotations')
    ff.create_directory(output_directory)

    print("\nLoading GenBank files...")
    # List and sort GenBank files
    gbk_files: List[str] = [os.path.join(genbank_files, f) for f in os.listdir(genbank_files)]
    gbk_files.sort()

    print('\nExtracting protein annotations from GenBank files...')
    # Initialize variables
    i: int = 0
    all_cds_info: Dict[str, List[str]] = {}
    all_translation_dict: Dict[str, str] = {}
    processed_proteins: Set[str] = set()
    same_protein_other_annotations: Dict[str, List[Tuple[str, str, str]]] = {}
    total_proteins: int = 0
    all_genbank_files_ids: Dict[str, List[str]] = {}

    # Parse GenBank files and extract protein annotations
    for f in gbk_files:
        file_name: str = os.path.basename(f)
        print(f"\rExtracting protein annotations from GenBank file: {file_name}", end='', flush=True)
        recs: List[SeqRecord] = [rec for rec in SeqIO.parse(f, 'genbank')]
        for r in recs:
            cds_info: Dict[str, List[str]] = get_protein_annotation_fasta(r)
            total_proteins += len(cds_info)
            for id_, info in cds_info.items():
                # Skip if protein has no id and no sequence
                if not id_:
                    continue
                # Add all IDs in the Genbank file
                all_genbank_files_ids.setdefault(file_name.removesuffix('.gbff'), []).append(id_)
                product: str
                gene: str
                translated_protein: str
                product, gene, translated_protein = info
                protein_hash: str = sf.hash_sequence(translated_protein)
                if protein_hash in processed_proteins:
                    # If the protein has already been processed, save the other annotations that it may have
                    same_protein_other_annotations.setdefault(protein_hash, []).append((id_, product, gene))
                    continue
                else:
                    processed_proteins.add(protein_hash)
                    all_translation_dict[id_] = translated_protein

            all_cds_info.update(cds_info)

    print(f"\nOut of {total_proteins} proteins sequences {len(all_translation_dict)} are unique proteins\n")
    
    print("Clustering protein sequences...")
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
    print(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} clusters.\n")

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

    print("\nRunning BLASTp...")
    blastp_results_folder: str = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, List[Union[str, float]]] = {}
    best_bsr_values_per_genbank_file: Dict[str, Dict[str, List[Union[str, float]]]] = {k: {} for k in all_genbank_files_ids.keys()}
    total_blasts: int = len(reps_ids)
    i: int = 1

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(
            bf.run_blastdb_multiprocessing,
            repeat(get_blastp_exec),
            repeat(blast_db_files),
            translations_paths.values(),
            translations_paths.keys(),
            repeat(blastp_results_folder)
        ):
            # Get the alignments
            filtered_alignments_dict: Dict[str, Dict[str, Dict[str, Union[int, float]]]]
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

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
                    # Check if the BSR value is the best for the locus
                    current_best_bsr: List[Union[str, float]] = best_bsr_values.get(loci)
                    # If there is a previous BSR value for the locus, check if the current BSR value is higher
                    if current_best_bsr and bsr_value > current_best_bsr[1]:
                        best_bsr_values[loci] = [subject_id, bsr_value]
                    else:
                        best_bsr_values[loci] = [subject_id, bsr_value]

                    # Get best value for genbank file
                    genbank_file: str = itf.identify_string_in_dict_get_key(subject_id, all_genbank_files_ids)
                    current_best_in_genbank_file: List[Union[str, float]] = best_bsr_values_per_genbank_file[genbank_file].get(loci)
                    if current_best_in_genbank_file and bsr_value > current_best_in_genbank_file[1]:
                        best_bsr_values_per_genbank_file[genbank_file][loci] = [subject_id, bsr_value]
                    else:
                        best_bsr_values_per_genbank_file[genbank_file][loci] = [subject_id, bsr_value]

            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    # Save annotations
    header: str = 'Locus\tgenebank_origin_id\tgenebank_origin_product\tgenebank_origin_name\tBSR'
    annotations_file: str = os.path.join(output_directory, 'genbank_annotations.tsv')
    not_matched_or_bsr_failed_loci: set = set(translations_paths.keys()) - set(best_bsr_values.keys())
    with open(annotations_file, 'w') as at:
        at.write(header + '\n')
        for loci, subject_info in best_bsr_values.items():
            subject_id: str = subject_info[0]  # Get the original ID and not the modified Blast version
            bsr_value: float = subject_info[1]  # Get the BSR value
            product: str = all_cds_info[subject_info[0]][0]  # Get the product name
            gene: str = all_cds_info[subject_info[0]][1]  # Get the gene name
            # Check if the gene name is empty
            if gene == '':
                gene = 'NA'
            # Write the annotations to the file
            at.write(f"{loci}\t{subject_id}\t{product}\t{gene}\t{bsr_value}\n")
        # Write the loci that did not match or the BSR value was lower than the threshold
        for loci in not_matched_or_bsr_failed_loci:
            at.write(f"{loci}\tNA\tNA\tNA\tNA\n")

    best_annotations_per_genbank_file: str = os.path.join(output_directory, "best_annotations_per_genbank_file")
    ff.create_directory(best_annotations_per_genbank_file)
    for file, loci_results in best_bsr_values_per_genbank_file.items():
        annotations_file_genbank: str = os.path.join(best_annotations_per_genbank_file, f"{file}.tsv")
        with open(annotations_file_genbank, 'w') as at:
            for loci, subject_info in loci_results.items():
                subject_id: str = subject_info[0]  # Get the original ID and not the modified Blast version
                bsr_value: float = subject_info[1]  # Get the BSR value
                product: str = all_cds_info[subject_info[0]][0]  # Get the product name
                gene: str = all_cds_info[subject_info[0]][1]  # Get the gene name
                # Check if the gene name is empty
                if gene == '':
                    gene = 'NA'
                # Write the annotations to the file
                at.write(f"{loci}\t{subject_id}\t{product}\t{gene}\t{bsr_value}\n")
        
    return annotations_file
            