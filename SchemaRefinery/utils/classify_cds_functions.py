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

def identify_problematic_cds(cds_presence_in_genomes, cds_translation_dict, protein_hashes, not_included_cds, cds_output,
                             bsr_value, dropped_cds, size_threshold, cpu):
    """
    Identifies problematic CDS (Coding DNA Sequences) based on specified criteria and outputs the results
    and Remove the instace of CDS from all the dicts.

    Parameters
    ----------
    cds_presence_in_genomes : dict
        A dictionary mapping each genome to the presence data of CDS.
    cds_translation_dict : dict
        A dictionary mapping CDS identifiers to their translated protein sequences.
    protein_hashes : set
        A set of unique hashes representing protein sequences, used to identify duplicates or problematic sequences.
    not_included_cds : set
        A set of CDS identifiers that are not included in the analysis, potentially due to previous filtering.
    cds_output : str
        The file path where the results of the analysis will be saved.
    bsr_value : float
        The BLAST Score Ratio (BSR) threshold used to determine problematic sequences. Sequences with a BSR below this value may be considered problematic.
    dropped_cds : dict
        A dictionary mapping all of the dropped CDSs to the cause of drop.
    cpu : int
        The number of CPU cores to be used for parallel processing tasks within the function.

    Returns
    -------
    None
        This function does not return a value but writes the results of the analysis to the specified output file.

    Notes
    -----
    This function is part of a larger pipeline for analyzing genomic data, specifically focusing on the identification
    of problematic CDS based on duplication, absence in certain genomes, or low BLAST Score Ratios. The results are used
    to refine the dataset for further analysis.
    """

    print("\nIdentifying possible NIPHEMs...")
    # Identify NIPHEMs.
    same_origin_genome = {}
    niphems_presence_in_genome = {}
    only_niphems_in_genomes = {}
    # Iterate over each CDS and check for NIPHEMs.
    for id_, cds_in_genomes in cds_presence_in_genomes.items():
        # Remove the protein number from the ID.
        genome_id = itf.remove_by_regex(id_, r'-protein\d+')
        same_origin_genome.setdefault(genome_id, set()).add(id_)
        # If there are duplicates in genomes.
        if len(cds_in_genomes) != len(set(cds_in_genomes)):
            # If all of the genomes contain only NIPHEMs in genomes.
            if itf.check_if_all_elements_are_duplicates(cds_in_genomes):
                # Save the CDS that are only NIPHEMs in genomes.
                only_niphems_in_genomes.setdefault(id_, set(cds_in_genomes))
                dropped_cds.setdefault(id_, 'Dropped_due_to_being_only_NIPHEM_in_genomes')
            else:
                # Add to the dict for further processing when to calculate if to exclude
                # possible new loci
                niphems_presence_in_genome.setdefault(id_, cds_in_genomes)
    # Write the identified NIPHEMs to a file.
    niphems_file = os.path.join(cds_output, 'identified_NIPHEMs_CDSs.tsv')
    tab = "\t"
    with open(niphems_file, 'w') as niphems:
        niphems.write('CDS_ID\tGenome_presence:\n')
        for cds, genomes_id in only_niphems_in_genomes.items():
            niphems.write(f"{cds}{tab}{tab.join([str(i) for i in genomes_id])}\n")
    # Print the results.
    print(f"There were identified {len(niphems_presence_in_genome)} CDSs containing NIPHEMs in genomes "
          f"and {len(only_niphems_in_genomes)} were removed for being present in genomes that only contain NIPHEMs.")

    # Identify CDSs present in the same genome.
    same_origin_genome = {genome_id: [cds for cds in cds_ids if cds_translation_dict.get(cds)] for genome_id, cds_ids in same_origin_genome.items()}
    # Filter out genomes with only one CDS.
    same_origin_genome = {genome_id: cds_ids for genome_id, cds_ids in same_origin_genome.items() if len(cds_ids) > 1}
    # Identify which cases to run.
    sequences_to_run = {}
    sequences_to_run_against = {}
    for genome_id, cds_ids in same_origin_genome.items():
        for cds in cds_ids:
            if cds_translation_dict.get(cds):
                sequences_to_run.setdefault(cds, cds_translation_dict[cds])
                sequences_to_run_against.setdefault(cds, [cds_id for cds_id in cds_ids if (cds_id != cds and cds_translation_dict.get(cds_id))])
    # Filter out the sequences that are from the same genome however don't have similiar protein size.
    for cds, cds_ids_to_run_against in list(sequences_to_run_against.items()):
        cds_protein_len = len(cds_translation_dict[cds])
        min_len = cds_protein_len - size_threshold * cds_protein_len
        max_len = cds_protein_len + size_threshold * cds_protein_len
        for member_id in list(cds_ids_to_run_against):
            member_protein_len = len(cds_translation_dict[member_id])
            if not max_len >= member_protein_len >= min_len:
                sequences_to_run_against[cds].remove(member_id)
    # Remove empty dicts.
    itf.remove_empty_dicts_recursive(sequences_to_run_against)

    #Create folders.
    niphs_folder = os.path.join(cds_output, 'NIPHs_and_NIPHEMs_processing')
    ff.create_directory(niphs_folder)
    translation_sequences_folder = os.path.join(niphs_folder, 'translation_sequences')
    ff.create_directory(translation_sequences_folder)
    translation_sequences_to_run_against_folder = os.path.join(niphs_folder, 'translation_sequences_to_run_against')
    ff.create_directory(translation_sequences_to_run_against_folder)

    #Write all of the FASTAs.
    sequences_fasta_path = {}
    to_run_against_paths = {}
    # Write the FASTAs for the CDSs and the CDSs to run against.
    for cds, cds_ids_to_run_against in sequences_to_run_against.items():
        member_file = os.path.join(translation_sequences_folder, f"{cds}.fasta")
        sequences_fasta_path[cds] = member_file
        #FASTAs to run.
        with open(member_file, 'w') as m_file:
            m_file.write(f">{cds}\n{sequences_to_run[cds]}\n")
        member_file = os.path.join(translation_sequences_to_run_against_folder, f"{cds}.fasta")
        to_run_against_paths[cds] = member_file
        #FASTAs proteins to run against
        for member_id in cds_ids_to_run_against:
            write_type = 'w' if not os.path.exists(member_file) else 'a'
            with open(member_file, write_type) as m_file:
                m_file.write(f">{member_id}\n{sequences_to_run[member_id]}\n")

    # Run BLASTp to identify possible NIPHs.
    print("\nIdentifying possible NIPHs...")
    self_score_folder = os.path.join(niphs_folder, 'self_score')
    ff.create_directory(self_score_folder)
    # Get the path to the BLASTp executable.
    get_blastp_exec = lf.get_tool_path('blastp')
    i = 1
    # Create a dictionary to store the self-score of each CDS.
    self_score_dict_niphs = {}
    # Get the max length of the IDs.
    max_id_length = len(max(sequences_fasta_path))
    # Calculate self-score.
    print("Calculating self-score for possible NIPHs...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                sequences_fasta_path.keys(),
                                repeat(get_blastp_exec),
                                sequences_fasta_path.values(),
                                repeat(self_score_folder)):
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, True, True, True, False)
    
            # Save self-score.
            self_score_dict_niphs[res[0]] = self_score
                            
            print(f"\rRunning BLASTp to calculate self-score for possible NIPHs {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1
    # Run BLASTp to confirm possible NIPHs.
    niphs_blastp_results_folder = os.path.join(niphs_folder, 'niphs_blastp_results')
    ff.create_directory(niphs_blastp_results_folder)
    save_bsr_score = {}
    total_blasts = len(sequences_fasta_path)
    i = 1
    # Run Blastp and calculate BSR.
    print("\nRunning BLASTp to confirm possible NIPHs...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blast_fastas_multiprocessing,
                                sequences_fasta_path, 
                                repeat(get_blastp_exec),
                                repeat(niphs_blastp_results_folder),
                                repeat(sequences_fasta_path),
                                to_run_against_paths.values()):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, False, True, False, False)


            # Since BLAST may find several local aligments choose the first one (highest one) to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    #Highest score (First one)
                    subject_score = next(iter(results.values()))['score']
                    save_bsr_score.setdefault(query, {}).update({subject_id: bf.compute_bsr(subject_score, self_score_dict_niphs[query])})

            print(f"\rRunning BLASTp to confirm identified NIPHs: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    #Identify NIPHS
    niphs_in_genomes = {}
    #Filter the BSR score.
    filtered_save_bsr_score = {query: {subject_id: bsr for subject_id, bsr in subjects_ids.items() if bsr >= bsr_value} for query, subjects_ids in save_bsr_score.items()}
    #Remove empty dicts.
    itf.remove_empty_dicts_recursive(filtered_save_bsr_score)

    # When some IDs didnt get in the same group
    to_merge_lists = [[query] + [subject for subject in subjects_ids.keys()] for query, subjects_ids in filtered_save_bsr_score.items()]
    niphs_in_genomes = {index: set(value) for index, value in enumerate(cf.cluster_by_ids_bigger_sublists(to_merge_lists))}

    for index, niphs in list(niphs_in_genomes.items()):
        for cds in list(niphs):
            same_protein_ids = itf.identify_string_in_dict_get_value(cds, protein_hashes)
            niphs_in_genomes[index].update(same_protein_ids)

    #Write the identified NIPHs to a file.
    niphs_presence_in_genomes = {}
    count_niphs_groups = 0
    count_niphs_cds = 0
    total_niphs = len(niphs_in_genomes)
    niphs_file = os.path.join(cds_output, 'identified_NIPHs_CDSs.tsv')
    # Iterate over possible NIPHs and write to file.
    for niph_id, cds_ids in list(niphs_in_genomes.items()):
        temp_niph_holder = []
        # Get the presence of the CDS in the genomes.
        for cds_id in cds_ids:
            niphs_presence_in_genomes[cds_id] = cds_presence_in_genomes[cds_id]
            temp_niph_holder.append(set(cds_presence_in_genomes[cds_id]))
        # Check if all of the sets are the same.
        if itf.check_if_all_sets_are_same(temp_niph_holder):
            # Write to file.
            write_type = 'w' if not os.path.exists(niphs_file) else 'a'
            with open(niphs_file, write_type) as niphs:
                count_niphs_groups += 1
                # Remove from niphs_in_genomes
                del niphs_in_genomes[niph_id]
                if write_type == 'w':
                    niphs.write('CDS_ID\tGenome_presence:\n')
                for cds_id in cds_ids:
                    count_niphs_cds += 1
                    niphs.write(f"{cds_id}{tab}{tab.join([str(i) for i in niphs_presence_in_genomes[cds_id]])}\n")
                    # Remove from cds not included in the schema dict.
                    dropped_cds.setdefault(cds_id, 'Dropped_to_being_only_NIPH_in_genomes')
                niphs.write("\n")
    # Convert the niphs_in_genomes dict to a list of sets.
    niphs_in_genomes = {key: list(value) for key, value in niphs_in_genomes.items()}
    print(f"There were Identified {total_niphs} groups of CDSs containing NIPHs and {count_niphs_groups}"
          f" groups ({count_niphs_cds} CDSs) were removed for being present in genomes that only contain NIPHs.")

    return niphems_presence_in_genome, niphs_in_genomes, niphs_presence_in_genomes

def identify_problematic_loci(niphems_presence_in_genome, niphs_in_genomes, niphs_presence_in_genomes,
                            cds_presence_in_genomes, clusters_to_keep, clusters, problematic_proportion,
                            dropped_cds, drop_possible_loci, results_output):
    """
    Removes loci deemed problematic based on a specified proportion of NIPHS and NIPHEMS present in the genomes.

    Parameters
    ----------
    niphems_presence_in_genome : dict
        A dictionary mapping each genome to its Niphems presence data.
    niphs_in_genomes : dict
        A dictionary mapping each genome to its Niphs data.
    niphs_presence_in_genomes : dict
        A dictionary mapping each genome to the presence data of Niphs.
    cds_presence_in_genomes : dict
        A dictionary mapping each genome to the presence data of CDS (Coding Sequences).
    clusters_to_keep : list
        A list of CDS identifiers that should be retained.
    clusters : dict
        A dictionary mapping cluster identifiers to their respective genomic data.
    problematic_proportion : float
        The proportion threshold above which a locus is considered problematic.
    dropped_cds : dict
        A dictionary mapping all of the dropped CDSs to the cause of drop.
    drop_possible_loci : set
        A set of loci identifiers that have been dropped from the analysis.
    results_output : str
        The file path where the results of the analysis will be saved.

    Returns
    -------
    proportion_of_niph_genomes : dict
        A dictionary mapping each genome to the proportion of NIPHS and NIPHEMS present in it.
    dropped_due_to_niphs_or_niphems : set
        A set of loci identifiers that were dropped due to the presence of NIPHS or NIPHEMS.
    clusters_to_keep_all_members : dict
        A dictionary mapping each group to its member CDS identifiers.
    clusters_to_keep_all_genomes : dict
        A dictionary mapping each group to the genomes in which it is present.

    Notes
    -----
    This function is designed to work with genomic data, specifically focusing on the presence and absence of certain
    NIPHEMs, NIPHs. It filters out loci based on a defined problematic proportion and updates the genomic data structures
    accordingly.
    """
    # Create list to store the potential paralogous new loci.
    potential_paralagous_new_loci = []
    # Pre process NIPHs
    for key, niphs in list(niphs_in_genomes.items()):
        # Get ids without the allele identifier of CDSs NIPHs that were not removed.
        niphs_ids = [niph.split('_')[0] for niph in niphs if niph not in dropped_cds]
        # If there are no pairs of NIPHs
        if len(niphs_ids) < 2 :
            continue
        possible_loci_ids = niphs_ids
        id_class_1a = [itf.identify_string_in_dict_get_key(niph_id, clusters_to_keep['1a']) for niph_id in niphs_ids]
        # If all IDs are the same (mantain them).
        if all(id_class_1a) and len(set(id_class_1a)) == 1:
            continue
        elif len(set(niphs_ids)) == 1:
            continue
        # If IDs are different (meaning that they were not clustered together or joined).
        # we remove them and write to file the potential paralogous.
        else:
            if any(id_class_1a):
                # Get indices of all values that are not None
                indices_not_none = [index for index, value in enumerate(id_class_1a) if value is not None]
                # Replace the IDs with the joined IDs.
                for index in indices_not_none:
                    possible_loci_ids[index] = id_class_1a[index]
            
            potential_paralagous_new_loci.append(set(possible_loci_ids))

    # Get all of the genomes that one groups is present in.
    clusters_to_keep_all_members = {}
    clusters_to_keep_all_genomes = {}
    for class_, cds_group in list(clusters_to_keep.items()):
    # Iterate over each group in class.
        for group in list(cds_group):
            clusters_to_keep_all_members.setdefault(group, set())
            clusters_to_keep_all_genomes.setdefault(group, set())
            # If the group is a joined group.
            if class_ == '1a':
                for cds in cds_group[group]:
                    for cds_allele in clusters[cds]:
                        clusters_to_keep_all_members[group].add(cds_allele)
                        clusters_to_keep_all_genomes[group].update(cds_presence_in_genomes[cds_allele])
            # If the group is not a joined group.
            else:
                for cds_allele in clusters[group]:
                    clusters_to_keep_all_members[group].add(cds_allele)
                    clusters_to_keep_all_genomes[group].update(cds_presence_in_genomes[cds_allele])

    # Get all of the NIPHs and NIPHEMs in the genomes to consider.
    get_niphems_in_genomes = {}
    proportion_of_niph_genomes = {}
    genomes_that_are_niphs_and_niphems = {}
    # Process NIPHS.
    for niphs in niphs_in_genomes.values():
        intersection_set = None
        if len([niph for niph in niphs if niph not in dropped_cds]) < 2:
            continue
        for niph in niphs:
            if niph in dropped_cds:
                continue
            # Get the ID without the allele identifier.
            niph_genome_id = niph.split('_')[0]
            # Get the ID for the joined IDs.
            id_class_1a = itf.identify_string_in_dict_get_key(niph_genome_id, clusters_to_keep['1a'])
            # Get all of the IDs of the genomes that intersect only in the NIPHs (two similiar alleles present in the same genomes).
            if not intersection_set:
                # Get the first set.
                intersection_set = set(niphs_presence_in_genomes[niph])
            else:
                # Get the intersection of the sets.
                intersection_set.intersection_update(niphs_presence_in_genomes[niph])
        # joined group or are in the same cluster.
        if intersection_set:
            genomes_that_are_niphs_and_niphems.setdefault(id_class_1a or niph_genome_id, set()).update(intersection_set)

    # Process NIPHEMs.
    for niphem in list(niphems_presence_in_genome):
        if niphem in dropped_cds:
            continue
        niphem_genome_id = niphem.split('_')[0]
        id_class_1a = itf.identify_string_in_dict_get_key(niphem_genome_id, clusters_to_keep['1a'])
        # Add which genomes are present in duplicate for that allele (two or more of the same genome ID).
        ids_of_genomes = itf.get_duplicates(niphems_presence_in_genome[niphem])
        get_niphems_in_genomes.setdefault(id_class_1a or niphem_genome_id, []).append(ids_of_genomes)
        
        # Add identified NIPHEMs to the dict that contains the NIPHs and NIPHEMs.
        # Add the genomes that are NIPHEMs.
        genomes_that_are_niphs_and_niphems.setdefault(id_class_1a or niphem_genome_id, set()).update(set(ids_of_genomes))

    # Get the proportion of NIPHs and NIPHEMs in the genomes for each group were they are present.
    dropped_due_to_niphs_or_niphems = set()
    for key, genomes in list(genomes_that_are_niphs_and_niphems.items()):
        # Get the proportion of NIPHs and NIPHEMs in the genomes.
        proportion = len(genomes) / len(clusters_to_keep_all_genomes[key])
        proportion_of_niph_genomes.setdefault(key, proportion)
        # If the proportion is greater than the threshold, remove the group.
        if proportion >= problematic_proportion:
            # Save cases to drop
            dropped_due_to_niphs_or_niphems.add(key)
            drop_possible_loci.add(key)
    
    # Create file Path.
    potential_paralagous = os.path.join(results_output, 'potential_paralagous.tsv')
    # Write paralogous loci
    write_type = 'w'
    with open(potential_paralagous, write_type) as potential_paralagous_report:
        for paralagous_group in potential_paralagous_new_loci:
            paralagous_group = [group for group in paralagous_group if group not in drop_possible_loci]
            if len(paralagous_group) >= 2:
                potential_paralagous_report.write('\t'.join(paralagous_group) + '\n')

    # Write the groups that were removed due to the presence of NIPHs or NIPHEMs.
    niphems_and_niphs_file = os.path.join(results_output, 'niphems_and_niphs_groups.tsv')
    with open(niphems_and_niphs_file, 'w') as niphems_and_niphs:
        niphems_and_niphs.write('Group_ID\tProportion_of_NIPHs_and_NIPHEMs\tOutcome\n')
        for group, proportion in proportion_of_niph_genomes.items():
            niphems_and_niphs.write(f"{group}\t{proportion}\t{'Dropped' if group in dropped_due_to_niphs_or_niphems else 'Kept'}\n")

    return proportion_of_niph_genomes, dropped_due_to_niphs_or_niphems, clusters_to_keep_all_members, clusters_to_keep_all_genomes

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

def remove_dropped_cds_from_analysis(dropped_cds, not_included_cds, niphems_presence_in_genome,
                                     cds_translation_dict, niphs_presence_in_genomes, protein_hashes,
                                     niphs_in_genomes):
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
        if translation_hash:
            dropped_cds[dropped_id] = reason
            niph_group_id = itf.identify_string_in_dict_get_key(dropped_id, niphs_in_genomes)
            if niph_group_id:
                niphs_in_genomes[niph_group_id].remove(dropped_id)
                if len(niphs_in_genomes[niph_group_id]) < 2:
                    del niphs_in_genomes[niph_group_id]
            if not_included_cds.get(dropped_id):
                del not_included_cds[dropped_id]
            # Remove from NIPHEMs dict
            if niphems_presence_in_genome.get(dropped_id):
                del niphems_presence_in_genome[dropped_id]
            # If this CDSs is the representative in translation dict
            if cds_translation_dict.get(dropped_id):
                del cds_translation_dict[dropped_id]
            # Remove from NIPHs
            if niphs_presence_in_genomes.get(dropped_id):
                del niphs_presence_in_genomes[dropped_id]
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
        if dropped_1a:
            for rep_id in dropped_1a:
                for cds_id in clusters[rep_id]:
                    dropped_cds[cds_id] = reason
        else:
            for cds_id in clusters[drop_id]:
                dropped_cds[cds_id] = reason

def replace_ids_in_clusters(clusters, frequency_cds, dropped_cds, not_included_cds, prot_len_dict,
                            cds_translation_dict, protein_hashes, cds_presence_in_genomes,
                            niphems_presence_in_genome, niphs_presence_in_genomes, niphs_in_genomes,
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
            if dropped_id:
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
            # Replace the old ID with the new ID for the NIPHEMs dict.
            if niphems_presence_in_genome.get(member):
                niphems_presence_in_genome[new_id] = niphems_presence_in_genome.pop(member)
            # Replace the old ID with the new ID for the NIPHs genome dict.
            if niphs_presence_in_genomes.get(member):
                niphs_presence_in_genomes[new_id] = niphs_presence_in_genomes.pop(member)
            # Replace the old ID with the new ID for the CDSs that matched as NIPHs.
            niphs_group_id = itf.identify_string_in_dict_get_key(member, niphs_in_genomes)
            if niphs_group_id or niphs_group_id == 0:
                niphs_group_member_index = niphs_in_genomes[niphs_group_id].index(member)
                niphs_in_genomes[niphs_group_id][niphs_group_member_index] = new_id
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

def process_new_loci(fastas_folder, constants):
    new_loci_folder = os.path.join(fastas_folder, 'new_possible_loci_fastas')
    possible_new_loci = {fastafile: os.path.join(new_loci_folder, fastafile) for fastafile in os.listdir(new_loci_folder) if fastafile.endswith('.fasta')}
    possible_new_loci_short_dir = os.path.join(new_loci_folder, 'short')
    possible_new_loci_short = {fastafile: os.path.join(possible_new_loci_short_dir, fastafile) for fastafile in os.listdir(possible_new_loci_short_dir) if fastafile.endswith('.fasta')}
    master_file_path = os.path.join(fastas_folder, 'master.fasta')
    possible_new_loci_translation_folder = os.path.join(fastas_folder, 'possible_new_loci_translation_folder')
    ff.create_directory(possible_new_loci_translation_folder)

    alleles = {}
    translation_dict_possible_new_loci = {}
    for new_loci in possible_new_loci.values():
        loci_id = ff.file_basename(new_loci).split('.')[0]
        alleles.setdefault(loci_id, {})
        fasta_dict = sf.fetch_fasta_dict(new_loci, False)
        os.remove(new_loci)
        for allele_id, sequence in fasta_dict.items():
            new_allele_id = f"{loci_id}_{allele_id}"
            alleles.setdefault(loci_id, {}).update({new_allele_id: str(sequence)})
            write_type = 'a' if os.path.exists(new_loci) else 'w'
            with open(new_loci, write_type) as new_loci_file:
                new_loci_file.write(f">{new_allele_id}\n{str(sequence)}\n")

            write_type = 'a' if os.path.exists(master_file_path) else 'w'
            with open(master_file_path, write_type) as master_file:
                master_file.write(f">{new_allele_id}\n{str(sequence)}\n")
        
        fasta_dict = sf.fetch_fasta_dict(new_loci, False)

        trans_path_file = os.path.join(possible_new_loci_translation_folder, f"{loci_id}.fasta")

        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        trans_path_file,
                                                        None,
                                                        constants[5],
                                                        False,
                                                        constants[6],
                                                        False)
        
        translation_dict_possible_new_loci.update(trans_dict)

    for new_loci_reps in possible_new_loci_short.values():
        loci_id = ff.file_basename(new_loci_reps).split('_')[0]
        fasta_dict = sf.fetch_fasta_dict(new_loci_reps, False)
        os.remove(new_loci_reps)
        for allele_id, sequence in fasta_dict.items():
            new_allele_id = f"{loci_id}_{allele_id}"
            alleles.setdefault(loci_id, {}).update({new_allele_id: sequence})
            write_type = 'a' if os.path.exists(new_loci_reps) else 'w'
            with open(new_loci_reps, write_type) as new_loci_reps_file:
                new_loci_reps_file.write(f">{new_allele_id}\n{str(sequence)}\n")
                
    return alleles, master_file_path, possible_new_loci, translation_dict_possible_new_loci
    