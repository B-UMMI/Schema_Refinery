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

def identify_problematic_new_loci(drop_possible_loci, clusters_to_keep, clusters, cds_present,
                                  not_included_cds, constants, output_path):
    """
    Identify problematic new loci based on the presence of NIPHs and NIPHEMs.

    Parameters
    ----------
    drop_possible_loci : set
        A set to store the IDs of clusters that should be dropped.
    clusters_to_keep : dict
        A dictionary where keys are class labels and values are dictionaries with cluster IDs as keys 
        and lists of cluster members as values, for key 1a the values are dicts with joined cluster
        ID as key and list as value.
    clusters : dict
        A dictionary where keys are cluster IDs and values are lists of allele IDs.
    cds_present : str
        Path to the distinct.hashtable file.
    not_included_cds : dict
        A dictionary where keys are allele IDs and values are sequences not included in the the schema.
    constants : list
        A list of constants used in the function. The 9th element (index 8) is the threshold for 
        problematic proportion.
    output_path : str
        The path to the directory where the output file will be saved.

    Returns
    -------
    drop_possible_loci : set
        The updated set of cluster IDs that should be dropped.
    """
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present) # Dict to store the decoded sequences.
    niphems_in_possible_new_loci = {} # Dict to store the NIPHEMs in each possible new loci.
    niphs_in_possible_new_loci = {} # Dict to store the NIPHs in each possible new loci.
    temp_niphs_in_possible_new_loci = {} # Temp dict to store the NIPHs in each possible new loci.
    total_possible_new_loci_genome_presence = {} # Dict to store the total number of genomes in each possible new loci.
    problematic_loci = {} # Dict to store the proportion of problematic genomes in each possible new loci.
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
            
    return drop_possible_loci

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

def write_dropped_cds_to_file(dropped_cds, results_output):
    """
    Write dropped CDS to file.

    Parameters
    ----------
    dropped_cds : dict
        The dictionary containing the dropped CDSs.
    results_output : str
        The path where the output will be written.

    Returns
    -------
    None, writes to file.
    """
    dropped_cds_output = os.path.join(results_output, 'dropped_cds.tsv')
    with open(dropped_cds_output, 'w') as dropped_cds_file:
        dropped_cds_file.write('CDS_ID\tReason_for_dropping\n')
        for cds_id, reason in dropped_cds.items():
            dropped_cds_file.write(f"{cds_id}\t{reason}\n")

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
            
    # Prepare to write the ID changes to a file
    tab = "\t"
    id_changes_file = os.path.join(results_output, 'cds_id_changes.tsv')
    
    # Open the file and write the header and ID changes
    with open(id_changes_file, 'w') as id_changes:
        id_changes.write('Original_ID\tID_after_clustering\tID_after_joining\n')
        for original_ids, changed_ids in cds_original_ids.items():
            # Write each original ID and its changed IDs to the file
            id_changes.write(f"{original_ids}\t{tab.join(changed_ids)}\n")

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
    cds_translation_dict : dict
        Dictionary mapping CDS to their translations.
    protein_hashes : dict
        Dictionary mapping protein hashes to their associated CDS.

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
