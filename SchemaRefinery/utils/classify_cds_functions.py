import os
import concurrent.futures
import copy
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       blast_functions as bf,
                       alignments_functions as af,
                       iterable_functions as itf,
                       linux_functions as lf,)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      iterable_functions as itf,
                                      linux_functions as lf,)

def identify_problematic_new_loci(clusters_to_keep, clusters, cds_present, not_included_cds, constants, output_path):
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)
    niphems_in_possible_new_loci = {}
    niphs_in_possible_new_loci = {}
    temp_niphs_in_possible_new_loci = {}
    total_possible_new_loci_genome_presence = {}
    problematic_loci = {}
    drop_possible_loci = set()
    for class_, cluster_keep in clusters_to_keep.items():
        for cluster_id in cluster_keep:
            niphems_in_possible_new_loci.setdefault(cluster_id, [])
            temp_niphs_in_possible_new_loci.setdefault(cluster_id, {})
            if class_ == '1a':
                cluster = itf.flatten_list([clusters[i] for i in clusters_to_keep['1a'][cluster_id]])
            else:
                cluster = clusters[cluster_id]
            for allele_id in cluster:
                sequence = not_included_cds[allele_id]
                hashed_seq = sf.seq_to_hash(str(sequence))
                allele_presence_in_genomes = decoded_sequences_ids[hashed_seq][1:]
                # NIPHs
                temp_niphs_in_possible_new_loci[cluster_id].setdefault(allele_id, set(allele_presence_in_genomes))
                # NIPHEMS
                # Convert list to set
                unique_elements = set(allele_presence_in_genomes)
                if len(unique_elements) == len(allele_presence_in_genomes):
                    continue
                else:
                    niphems_genomes = itf.get_duplicates(allele_presence_in_genomes)
                    niphems_in_possible_new_loci.setdefault(cluster_id, []).extend(niphems_genomes)

            niphs_in_genomes = set(itf.get_shared_elements(temp_niphs_in_possible_new_loci))
            niphs_in_possible_new_loci.setdefault(cluster_id, niphs_in_genomes)
            
            niphems_in_genomes = set(niphems_in_possible_new_loci[cluster_id])

            problematic_genomes_in_possible_new_loci = niphs_in_genomes | niphems_in_genomes

            total_possible_new_loci_genome_presence = len(set(itf.flatten_list(temp_niphs_in_possible_new_loci[cluster_id].values())))
            
            problematic_proportion = len(problematic_genomes_in_possible_new_loci)/total_possible_new_loci_genome_presence
            problematic_loci.setdefault(cluster_id, problematic_proportion)
            if problematic_proportion >= constants[8]:
                drop_possible_loci.add(cluster_id)

    
    # Write the groups that were removed due to the presence of NIPHs or NIPHEMs.
    niphems_and_niphs_file = os.path.join(output_path, 'niphems_and_niphs_groups.tsv')
    with open(niphems_and_niphs_file, 'w') as niphems_and_niphs:
        niphems_and_niphs.write('Group_ID\tProportion_of_NIPHs_and_NIPHEMs\tOutcome\n')
        for group, proportion in problematic_loci.items():
            niphems_and_niphs.write(f"{group}\t{proportion}\t{'Dropped' if group in drop_possible_loci else 'Kept'}\n")
            
    return problematic_loci, drop_possible_loci

def write_cluster_members_to_file(output_path, clusters_to_keep, clusters, frequency_in_genomes, drop_possible_loci):
    """
    Write cluster members to file.

    Parameters
    ----------
    output_path : str
        The path where the output will be written.
    clusters_to_keep : dict
        The dictionary containing the CDSs to keep.
    clusters : dict
        The dictionary containing the clusters.
    frequency_in_genomes : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.

    Returns
    -------
    None, writes to file.
    """
    write_cds = clusters_to_keep
    write_cds.setdefault('Dropped', drop_possible_loci)
    cluster_members_output = os.path.join(output_path, 'cluster_members.tsv')
    with open(cluster_members_output, 'w') as cluster_members_file:
        cluster_members_file.write('Cluster_ID\tRepresentatives_IDs\tRep_cluster_members\tFrequency_of_rep'
                                   '\tClassification\n')
        for class_, cds_list in clusters_to_keep.items():
            for cds in cds_list:
                classification = class_
                if class_ == '1a':
                    cluster_members_file.write(str(cds))
                    cds = clusters_to_keep[class_][cds]
                else:
                    cluster_members_file.write(cds)
                    cds = [cds]
                for rep_id in cds:
                    cluster_members_file.write('\t' + str(rep_id))
                    cds_ids = [cds_id for cds_id in clusters[rep_id]]
                    for count, cds_id in enumerate(cds_ids):
                        if count == 0:
                            cluster_members_file.write('\t' + cds_id + '\t' + str(frequency_in_genomes[rep_id])
                                                       + '\t' + classification + '\n')
                        else:
                            cluster_members_file.write('\t\t' + cds_id + '\n')

def update_ids_and_save_changes(clusters_to_keep, clusters, cds_original_ids, dropped_cds,
                                not_included_cds, results_output):
    """
    Update the IDs based on clustering and joining operations and save the changes.

    This function iterates through each class and its corresponding group of CDS (Coding DNA Sequences) to keep,
    updates the IDs based on the provided clusters and the original to new ID mappings, and saves the final ID changes
    to a TSV (Tab-Separated Values) file in the specified output directory.

    Parameters
    ----------
    clusters_to_keep : dict
        A dictionary where each key is a class and each value is a group of CDS to keep.
    clusters : dict
        A dictionary mapping representative IDs to their cluster members.
    cds_original_ids : dict
        A dictionary mapping original IDs to their new IDs after processing.
    dropped_cds : dict
        A dictionary mapping all of the dropped CDSs to the cause of drop.
    not_included_cds : dict
        Dict that contains DNA sequences for each CDS.
    results_output : str
        The directory path where the ID changes file will be saved.

    Notes
    -----
    The function iterates through the `clusters_to_keep` dictionary, updating IDs for each CDS based on their membership
    in the provided `clusters`. It generates a new ID for each CDS, updates `cds_original_ids` with these new IDs,
    and writes the original and new IDs to a TSV file named 'cds_id_changes.tsv' in the `results_output` directory.

    The ID updating process involves generating a new ID by appending an index to the main representative ID for each
    CDS in a cluster. This index is incremented for each CDS in the cluster.

    Examples
    --------
    Assuming the existence of appropriate dictionaries for `clusters_to_keep`, `clusters`, `cds_original_ids`, and a valid
    path for `results_output`, the function can be called as follows:

    >>> update_ids_and_save_changes(clusters_to_keep, clusters, cds_original_ids, '/path/to/output')
    
    This would process the IDs as described and save the changes to '/path/to/output/cds_id_changes.tsv'.
    """

    # Iterate through each class and its CDS group
    for class_, cds_group in clusters_to_keep.items():
        for cds in cds_group:
            main_rep = cds # The main representative ID for the CDS group
            
            # If the class is not '1a', treat the CDS as a single-element list
            if class_ != '1a':
                cds = [cds]
            else:
                # For class '1a', get the CDS group from clusters_to_keep
                cds = clusters_to_keep[class_][cds]
            
            index = 1  # Initialize an index for creating new IDs
            
            # Iterate through each representative ID in the CDS group
            for rep_id in list(cds):
                # Get all CDS IDs in the cluster for the representative ID
                cds_ids = clusters[rep_id]
                # Delete clusters with old IDs
                del clusters[rep_id]
                # Create new rep ID
                # Iterate through each CDS ID in the cluster
                for cds_id in list(cds_ids):
                    # Skip cases
                    if cds_id in dropped_cds:
                        continue
                    if not clusters.get(rep_id):
                        clusters[rep_id] = []
                    # Create a new ID using the main representative ID and the index
                    new_id = f"{main_rep}_{index}"
                    # Update the original ID with the new ID in cds_original_ids
                    cds_id_first = itf.identify_string_in_dict_get_key(cds_id, cds_original_ids)
                    cds_id_second = itf.identify_string_in_dict_get_value(cds_id, cds_original_ids)[-1]
                    # Replace in FASTA dict
                    not_included_cds[new_id] = not_included_cds.pop(cds_id_second)
                    # Add new cluster ID
                    clusters[rep_id].append(new_id)
                    cds_original_ids[cds_id_first].append(new_id)
                    index += 1  # Increment the index for the next ID

    # Add why CDS was dropped
    for cds_member, cause in dropped_cds.items():
        # Get the original ID for the CDS member
        cds_id = itf.identify_string_in_dict_get_key(cds_member, cds_original_ids) or cds_member
        # Append the cause of dropping to the original ID
        if cds_original_ids.get(cds_id):
            if len(cds_original_ids[cds_id]) == 2:
                cds_original_ids[cds_id].extend(['', cause])
            else:
                cds_original_ids[cds_id].append(cause)
        else:
            cds_original_ids.setdefault(cds_id, ['','',cause])
            
    # Prepare to write the ID changes to a file
    tab = "\t"
    id_changes_file = os.path.join(results_output, 'cds_id_changes.tsv')
    
    # Open the file and write the header and ID changes
    with open(id_changes_file, 'w') as id_changes:
        id_changes.write('Original_ID\tID_after_clustering\tID_after_joining\tDrop_reason\n')
        for original_ids, changed_ids in cds_original_ids.items():
            # Write each original ID and its changed IDs to the file
            id_changes.write(f"{original_ids}\t{tab.join(changed_ids)}\n")

def find_new_representatives(groups_trans_reps_paths, groups_trans, groups_paths_reps,
                             cpu, not_included_cds, constants, results_output):

    def run_blast_for_bsr(groups_trans, groups_trans_reps_paths, iterations_folder, cpu):
        print('\n')
        # Create a folder to store the results of the BLASTp self-score calculations.
        blastp_folder = os.path.join(iterations_folder, 'BLASTp')
        ff.create_directory(blastp_folder)
        save_bsr_score = {}
        total_blasts = len(groups_trans)
        i = 1
        # Run Blastp and calculate BSR.
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
            for res in executor.map(bf.run_blast_fastas_multiprocessing,
                                    groups_trans, 
                                    repeat(get_blastp_exec),
                                    repeat(blastp_folder),
                                    repeat(groups_trans_reps_paths),
                                    groups_trans.values()):
                
                filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, False, False)


                # Since BLAST may find several local aligments choose the first one (highest one) to calculate BSR.
                for query, subjects_dict in filtered_alignments_dict.items():
                    for subject_id, results in subjects_dict.items():
                        #Highest score (First one)
                        subject_score = next(iter(results.values()))['score']
                        save_bsr_score.setdefault(query, {}).update({subject_id: bf.compute_bsr(subject_score, self_score_dict_reps[query])})

                print(f"\rRunning BLASTp to confirm identified NIPHs: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
                i += 1

        return save_bsr_score
    
    def self_score_calc(temp_groups_trans_reps_paths, iterations_folder, cpu):
        # Create a folder to store the results of the BLASTp self-score calculations.
        self_score_folder = os.path.join(iterations_folder, 'self_score')
        ff.create_directory(self_score_folder)

        i = 1
        self_score_dict_reps = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
            for res in executor.map(bf.run_self_score_multiprocessing,
                                    temp_groups_trans_reps_paths.keys(),
                                    repeat(get_blastp_exec),
                                    temp_groups_trans_reps_paths.values(),
                                    repeat(self_score_folder)):
                
                _, self_scores, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, True, False, False, True)
        
                # Save self-score.
                self_score_dict_reps.update(self_scores)
                                
                print(f"\rRunning BLASTp to calculate self-score to identify new representatives {res[0]: <{max_id_length}}", end='', flush=True)
                i += 1
        return self_score_dict_reps


    bsr_value = constants[7]
    bsr_value_upper_rep = bsr_value + 0.1 if bsr_value <= 0.9 else 1
    blast_dir = os.path.join(results_output, 'BLAST_find_new_representatives')
    ff.create_directory(blast_dir)

    # Get the path to the BLASTp executable.
    get_blastp_exec = lf.get_tool_path('blastp')
    # Get the max length of the IDs.
    max_id_length = len(max(groups_trans))

    iteration = 1
    continue_to_run_blasts = True
    temp_group_paths = copy.deepcopy(groups_trans)
    temp_groups_trans_reps_paths = copy.deepcopy(groups_trans_reps_paths)
    # Loop till there are no more possible new reps cases.
    while continue_to_run_blasts:
        print(f"\nRunning BLASTp to identify new representatives: Iteration {iteration}...")
        # Create a folder to store the results of the BLASTp self-score calculations and BLASTp.
        iterations_folder = os.path.join(blast_dir, f"iteration_{iteration}")
        ff.create_directory(iterations_folder)
        # Run BLASTp to calculate self-score for the representatives.
        self_score_dict_reps = self_score_calc(temp_groups_trans_reps_paths, iterations_folder, cpu)
        # Run BLASTp to confirm possible new representatives.
        save_bsr_score = run_blast_for_bsr(temp_group_paths, temp_groups_trans_reps_paths, iterations_folder, cpu)

        # Filter the cases by BSR value in inverse order (lowest at the top).
        flattened = [
        (query, subject_id, bsr)
        for query, subjects_ids in save_bsr_score.items()
        for subject_id, bsr in subjects_ids.items()
        ]
        # Sort the cases by BSR value.
        sorted_flattened = sorted(flattened, key=lambda x: x[2], reverse=True)
        # Save the cases in a dict.
        sorted_save_bsr_score = {}
        for query, subject_id, bsr in sorted_flattened:
            if query not in sorted_save_bsr_score:
                sorted_save_bsr_score[query] = {}
            sorted_save_bsr_score[query][subject_id] = bsr

        # Identify cases that are in new reps threshold and cases that are in 
        # normal allele threshold
        new_reps_ids = {}
        del_matched_ids = {}
        for rep, alleles in sorted_save_bsr_score.items():
            loci = rep.split('_')[0]

            new_reps_ids.setdefault(loci, {})
            del_matched_ids.setdefault(loci, {})

            for allele, bsr in alleles.items():
                # Identify possible new reps base on bsr value
                if bsr_value_upper_rep >= bsr:
                    if not new_reps_ids[loci].get(allele):
                        new_reps_ids[loci].setdefault(allele, bsr)
                # We need only the cases where all of the matches are normal.
                else:
                    del_matched_ids[loci].setdefault(allele, bsr)

        # Remove the cases that are not in the range of bsr value to add as representative for other representatives.
        for loci, alleles in del_matched_ids.items():
            for allele in alleles:
                if new_reps_ids[loci].get(allele):
                    del new_reps_ids[loci][allele]

        # Clean the dict
        itf.remove_empty_dicts_recursive(new_reps_ids)   
        # Add one possible reps (not adding all because some may be for the same protein)
        for loci, alleles in list(new_reps_ids.items()):
            # Fetch translations of the alleles.
            translation = sf.read_fasta_file_dict(groups_trans[loci])
            # Fetch translation for the representatives.
            rep_translation = sf.read_fasta_file_dict(groups_trans_reps_paths[loci])
            for allele in list(alleles):
                # If the allele is already a representative, remove it.
                if allele in rep_translation:
                    del new_reps_ids[loci][allele]
                    continue
                # If the allele is not in the representatives, add it.
                elif allele not in rep_translation:
                    del new_reps_ids[loci][allele]
                    print(f"Representative for {loci} is {allele}")
                    with open(groups_paths_reps[loci], 'a') as fasta_file:
                        fasta_file.write(f">{allele}\n{str(not_included_cds[allele])}\n")
                    with open(groups_trans_reps_paths[loci], 'a') as fasta_file:
                        fasta_file.write(f">{allele}\n{str(translation[allele].seq)}\n")
                    # Since we want only one representative for each iteration
                    # because if we add all at the same time we may add several representatives
                    # that represent the same proteins.
                    break

        # Clean the dict
        itf.remove_empty_dicts_recursive(new_reps_ids)
        # What BLASTp to do again to confirm other possible reps
        temp_group_paths = {}
        temp_groups_trans_reps_paths = {}
        # Add new temp paths for loci to run again.
        for loci, alleles in new_reps_ids.items():
            temp_group_paths.update({loci: groups_trans[loci]})
            temp_groups_trans_reps_paths.update({loci: groups_trans_reps_paths[loci]})
        # If there are no more possible new reps to add.
        if not temp_group_paths:
            continue_to_run_blasts = False
        else:
            iteration += 1

def remove_dropped_cds_from_analysis(dropped_cds, not_included_cds,
                                     cds_translation_dict, protein_hashes):
    """
    Removes dropped CDS from the analysis based on the provided parameters.

    Parameters
    ----------
    dropped_cds : dict
        Dictionary containing dropped CDS with their reasons.
    not_included_cds : dict
        Dictionary of CDS that are not included in the analysis.
    niphems_presence_in_genome : dict
        Dictionary indicating the presence of NIPHEMs in the genome.
    cds_translation_dict : dict
        Dictionary mapping CDS to their translations.
    niphs_presence_in_genomes : dict
        Dictionary indicating the presence of NIPHs in genomes.
    protein_hashes : dict
        Dictionary mapping protein hashes to their associated CDS.
    niphs_in_genomes : dict
        Dictionary indicating the presence of NIPHs in genomes.

    Returns
    -------
    None
    """
    for dropped_id, reason in list(dropped_cds.items()):
        if cds_translation_dict.get(dropped_id):
            for similiar_protein_id in itf.identify_string_in_dict_get_value(dropped_id, protein_hashes):
                dropped_cds[similiar_protein_id] = reason

    for dropped_id, reason in list(dropped_cds.items()):
        #Remove from associated translation hashes dict.
        translation_hash = itf.identify_string_in_dict_get_key(dropped_id, protein_hashes)
        if translation_hash is not None:
            dropped_cds[dropped_id] = reason
            if not_included_cds.get(dropped_id):
                del not_included_cds[dropped_id]
            # If this CDSs is the representative in translation dict
            if cds_translation_dict.get(dropped_id):
                del cds_translation_dict[dropped_id]
            protein_hashes[itf.identify_string_in_dict_get_key(dropped_id, protein_hashes)].remove(dropped_id)
            # Remove all hash of the protein if it has no more CDSs associated with it.
            if len(protein_hashes[translation_hash]) == 0:
                del protein_hashes[translation_hash]

def add_cds_to_dropped_cds(drop_possible_loci, dropped_cds, clusters_to_keep,
                           clusters, reason, processed_drop):
    """
    Adds CDS to the dropped CDS list based on the provided parameters.

    Parameters
    ----------
    drop_possible_loci : list
        List of possible loci dropped.
    dropped_cds : dict
        Dictionary to store dropped CDS with their reasons.
    clusters_to_keep : dict
        Dictionary containing CDS to keep, classified by their class type.
    clusters : dict
        Dictionary containing clusters of CDS.
    cds_original_ids : dict
        Dictionary containing original IDs of CDS.
    reason : str
        Reason for dropping the CDS.
    processed_drop : list
        List of already processed drop IDs.

    Returns
    -------
    None
    """

    for drop_id in drop_possible_loci:
        if drop_id in processed_drop:
            continue
        else:
            processed_drop.append(drop_id)
            
        if itf.identify_string_in_dict_get_key(drop_id, clusters_to_keep['1a']):
            del clusters_to_keep['1a'][drop_id]
        else:
            class_ = itf.identify_string_in_dict_get_key(drop_id, {key: value for key, value in clusters_to_keep.items() if key != '1a'})
            if class_:
                clusters_to_keep[class_].remove(drop_id)

        dropped_1a = itf.identify_string_in_dict_get_value(drop_id, clusters_to_keep['1a'])
        if dropped_1a is not None:
            for rep_id in dropped_1a:
                for cds_id in clusters[rep_id]:
                    dropped_cds[cds_id] = reason
        else:
            for cds_id in clusters[drop_id]:
                dropped_cds[cds_id] = reason

def replace_ids_in_clusters(clusters, frequency_cds, dropped_cds, not_included_cds, prot_len_dict,
                            cds_translation_dict, protein_hashes, cds_presence_in_genomes,
                            reps_kmers_sim):
    """
    Replace the IDs of cluster alleles with new IDs in the format 'cluster_x' and update all relevant dictionaries.

    Parameters
    ----------
    clusters : dict
        Dictionary of clusters with their members.
    frequency_cds : dict
        Dictionary mapping CDS IDs to their frequencies.
    dropped_cds : dict
        Dictionary of dropped CDS IDs.
    not_included_cds : dict
        Dictionary of CDS IDs not included in the analysis.
    prot_len_dict : dict
        Dictionary mapping CDS IDs to their protein lengths.
    cds_translation_dict : dict
        Dictionary mapping CDS IDs to their translation sequences.
    protein_hashes : dict
        Dictionary mapping protein hashes to lists of CDS IDs.
    cds_presence_in_genomes : dict
        Dictionary mapping CDS IDs to their presence in genomes.
    niphems_presence_in_genome : dict
        Dictionary mapping NIPHEM IDs to their presence in genomes.
    niphs_presence_in_genomes : dict
        Dictionary mapping NIPH IDs to their presence in genomes.
    niphs_in_genomes : dict
        Dictionary mapping NIPH group IDs to lists of CDS IDs.
    reps_kmers_sim : dict
        Dictionary mapping representative cluster IDs to their k-mer similarities.

    Returns
    -------
    cds_original_ids : dict
        Dictionary mapping original CDS IDs to their new IDs.
    """
    cds_original_ids = {}
    # Replace the IDS of cluster alleles to x_1 and replace all of the alleles in
    # the variables.
    for cluster, members in list(clusters.items()):
        i = 1
        new_members_ids = []
        for member in list(members):
            # Get the new ID.
            new_id = f"{cluster}_{i}"
            # Add the new ID to the dict.
            cds_original_ids[member] = [new_id]
            # Replace the old ID with the new ID for frequency_cds.
            frequency_cds[new_id] = frequency_cds.pop(member)
            # Dropped cds
            dropped_id = itf.identify_string_in_dict_get_key(member, dropped_cds)
            if dropped_id is not None:
                dropped_cds[new_id] = dropped_cds.pop(dropped_id)
            # Save the new members IDs.
            new_members_ids.append(new_id)
            # Replace the old ID with the new ID for the DNA sequences.
            not_included_cds[new_id] = not_included_cds.pop(member)
            # Replace in hashes dict
            translation_hash = itf.identify_string_in_dict_get_key(member, protein_hashes)
            # Replace in prot_len_dict
            if prot_len_dict.get(member):
                prot_len_dict[new_id] = prot_len_dict.pop(member)

            index = protein_hashes[translation_hash].index(member)
            # Replace the old ID with the new ID for the translation sequences.
            # Since only representatives are in the dict we first check if it is present
            if cds_translation_dict.get(member):
                cds_translation_dict[new_id] = cds_translation_dict.pop(member)
            else: # Add the sequences previousy deduplicated
                rep_id = protein_hashes[translation_hash][0]
                cds_translation_dict[new_id] = cds_translation_dict[rep_id]

            # Replace the value at the found index
            protein_hashes[translation_hash][index] = new_id  # Replace `new_id` with the actual value you want to set
            # Replace the old ID with the new ID for the protein hashes.
            cds_presence_in_genomes[new_id] = cds_presence_in_genomes.pop(member)
            i += 1
        clusters[cluster] = new_members_ids
    # Change IDs in reps_kmers_sim
    for cluster_id, elements_id in list(reps_kmers_sim.items()):
        cluster_rep_id = f"{cluster_id}_1"
        reps_kmers_sim.setdefault(cluster_rep_id, {})
        for element_id, kmers in elements_id.items():
            if cds_original_ids.get(element_id):
                new_id = cds_original_ids[element_id][0]
            reps_kmers_sim[cluster_rep_id].setdefault(new_id, kmers)
        del reps_kmers_sim[cluster_id]

    return cds_original_ids

    