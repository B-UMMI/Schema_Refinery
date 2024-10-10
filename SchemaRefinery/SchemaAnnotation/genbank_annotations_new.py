import os
import concurrent.futures
from itertools import repeat 
from typing import Dict, List
from Bio import SeqIO

try:
    from utils import blast_functions as bf
    from utils.sequence_functions import translate_sequence
    from utils import file_functions as ff
    from utils import linux_functions as lf
    from utils import alignments_functions as af
    from utils import sequence_functions as sf
    from utils import clustering_functions as cf
except ModuleNotFoundError:
    from SchemaRefinery.utils import blast_functions as bf
    from SchemaRefinery.utils.sequence_functions import translate_sequence
    from SchemaRefinery.utils import file_functions as ff
    from SchemaRefinery.utils import linux_functions as lf
    from SchemaRefinery.utils import alignments_functions as af
    from SchemaRefinery.utils import sequence_functions as sf
    from SchemaRefinery.utils import clustering_functions as cf

def get_protein_annotation_fasta(seqRecord):
    """Get the translated protein from a Genbank file.

    Parameters
    ----------
    seqRecord : Biopython SeqRecord
        BioPython sequence record object.

    Returns
    -------
    fasta : list
        List containing the protein in fasta format.
    fasta_dict : dict
        Dict containing the translated protein as key and the values are
        list containing the protein_id, the product and the gene name.

    Notes
    -----
    Source: https://github.com/LeeBergstrand/Genbank-Downloaders/blob/d904c92788696b02d9521802ebf1fc999a600e1b/SeqExtract.py#L48
    """
    cds_info = {}
    features = seqRecord.features  # Each sequence has a list (called features) that stores seqFeature objects.
    for feature in features:  # For each feature on the sequence
        if feature.type == "CDS":  # CDS means coding sequence (These are the only features we're interested in)
            featQualifiers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
            # data about the sequence feature (for example the translation)

            # Gets the required qualifier. Uses featQualifers.get to return the qualifier or a default value if the quantifier
            # is not found. Calls strip to remove unwanted brackets and ' from qualifier before storing it as a string.
            protein_id = str(featQualifiers.get('protein_id', '')).strip('\'[]')

            if protein_id == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.

            gene = str(featQualifiers.get('gene', '')).strip('\'[]')
            product = str(featQualifiers.get('product', 'no_product_name')).strip('\'[]')
            translated_protein = str(featQualifiers.get('translation', 'no_translation')).strip('\'[]')

            cds_info[protein_id] = [product, gene, translated_protein]

    return cds_info

def genbank_annotations(genbank_files: str, schema_directory: str,
                        output_directory: str, cpu: int,
                        bsr: float, translation_table: int,
                        clustering_sim: float, clustering_cov: float,
                        size_ratio: float) -> str:
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
    output_directory = os.path.join(output_directory, 'genbank_annotations')
    ff.create_directory(output_directory)

    print("Loading GenBank files...")
    # List and sort GenBank files
    gbk_files: List[str] = [os.path.join(genbank_files, f) for f in os.listdir(genbank_files)]
    gbk_files.sort()
    
    print('\n')
    print('Extracting protein annotations from GenBank files...')
    # Parse GenBank files and extract protein annotations
    i = 0
    total_length = len(gbk_files)
    all_cds_info = {}
    all_translation_dict = {}
    processed_proteins = set()
    total_proteins = 0
    for f in gbk_files:
        print(f"\rExtracting protein annotations from GenBank file: {os.path.basename(f)}", end='', flush=True)
        recs = [rec for rec in SeqIO.parse(f, 'genbank')]
        for r in recs:
            cds_info = get_protein_annotation_fasta(r)
            total_proteins += len(cds_info)
            for id_, info in cds_info.items():
                # Skip if protein has no id and no sequence
                if not id_:
                    continue
                product, gene, translated_protein = info
                protein_hash = sf.hash_sequence(translated_protein)
                if protein_hash in processed_proteins:
                    continue
                else:
                    processed_proteins.add(protein_hash)
                    all_translation_dict[id_] = translated_protein

            all_cds_info.update(cds_info)

    print(f"\nOut of {total_proteins} proteins sequences {len(all_translation_dict)} are unique proteins\n")
    
    print("Clustering protein sequences...")
    all_alleles = {}
    reps_sequences = {}
    reps_groups = {}

    all_alleles, reps_sequences, reps_groups, prot_len_dict = cf.minimizer_clustering(
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
    print(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} custers.\n")
    # Save extracted protein sequences to a file
    blast_processing_folder = os.path.join(output_directory, 'blast_processing')
    ff.create_directory(blast_processing_folder)
    genbank_protein_file = os.path.join(blast_processing_folder, 'selected_genbank_proteins.fasta')
    with open(genbank_protein_file, 'w') as outfile:
        for protein_id, values in reps_sequences.items():
            outfile.write(f">{protein_id}\n{values}\n")
        
    # BLAST alleles for each locus against file with all CDSs from origin genomes
    short_folder = os.path.join(schema_directory, 'short')
    fasta_files_short_dict = {
        loci.split('.')[0].split('_')[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }
    print("Translating sequences...")
    reps_translations_folder = os.path.join(output_directory, 'reps_translations')
    ff.create_directory(reps_translations_folder)
    translation_dict = {}
    reps_ids = {}
    translations_paths = {}
    for loci, loci_path in fasta_files_short_dict.items():
        fasta_dict = sf.fetch_fasta_dict(loci_path, False)
        print(f"\rTranslating schema reps: {loci}", end='', flush=True)
        for allele_id, sequence in fasta_dict.items():
            reps_ids.setdefault(loci, []).append(allele_id)
            
            # Translate sequences and update translation dictionary
            trans_path_file = os.path.join(reps_translations_folder, f"{loci}.fasta")
            translations_paths[loci] = trans_path_file
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                            trans_path_file,
                                                            None,
                                                            0,
                                                            False,
                                                            translation_table,
                                                            False)
            translation_dict.update(trans_dict)
    print('\n')

    # Create BLASTdb
    blastdb_path = os.path.join(blast_processing_folder, 'blastdb')
    ff.create_directory(blastdb_path)
    blast_db_files = os.path.join(blastdb_path, 'genbank_protein_db')
    makeblastdb_exec = lf.get_tool_path('makeblastdb')
    bf.make_blast_db(makeblastdb_exec, genbank_protein_file, blast_db_files, 'prot')
    
    print("\nCalculate self-score for the CDSs...")
    self_score_reps_folder = os.path.join(blast_processing_folder, 'self_score_reps')
    ff.create_directory(self_score_reps_folder)
    self_score_dict = {}
    max_id_length = len(max(reps_ids))
    # Get Path to the blastp executable
    get_blastp_exec = lf.get_tool_path('blastp')
    i = 1
    # Calculate self-score
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                fasta_files_short_dict.keys(),
                                repeat(get_blastp_exec),
                                fasta_files_short_dict.values(),
                                repeat(self_score_reps_folder)):
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True, True, True, True)
    
            # Save self-score
            self_score_dict.update(self_score)
                            
            print(f"\rRunning BLASTp to calculate self-score for {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1
    # Print newline
    print('\n')
    
    print("Running BLASTp...")
    blastp_results_folder = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values = {}
    best_bsr_values = {}
    total_blasts = len(reps_ids)
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing,
                                repeat(get_blastp_exec),
                                repeat(blast_db_files),
                                translations_paths.values(),
                                translations_paths.keys(),
                                repeat(blastp_results_folder)):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

            # Since BLAST may find several local aligments choose the largest one to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                loci = query.split('_')[0]
                bsr_values.setdefault(query, {})
                for subject_id, results in subjects_dict.items():
                    if '|' in subject_id:
                        subject_id = subject_id.split('|')[1] # Get the original ID and not the modified Blast version
                    #Highest score (First one)
                    subject_score = next(iter(results.values()))['score']
                    bsr_value = bf.compute_bsr(subject_score, self_score_dict[query])
                    if bsr_value >= bsr:
                        bsr_values[query].update({subject_id: bsr})
                    else:
                        continue

                    current_best_bsr = best_bsr_values.get(loci)
                    if current_best_bsr and bsr > current_best_bsr[1]:
                        best_bsr_values[loci] = [subject_id, bsr]
                    else:
                        best_bsr_values[loci] = [subject_id, bsr]
                    
        
            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    # Save annotations
    header = 'Locus\tgenebank_origin_id\tgenebank_origin_product\tgenebank_origin_name\tgenebank_origin_bsr'
    annotations_file = os.path.join(output_directory, 'genbank_annotations.tsv')
    with open(annotations_file, 'w') as at:
        at.write(header + '\n')
        for loci, subject_info in best_bsr_values.items():
            subject_id = subject_info[0]
            bsr = subject_id[1]
            product = all_cds_info[subject_info[0]][0]
            gene = all_cds_info[subject_info[0]][1]
            
            at.write(f"{loci}\t{subject_id}\t{product}\t{gene}\t{bsr}\n")
            