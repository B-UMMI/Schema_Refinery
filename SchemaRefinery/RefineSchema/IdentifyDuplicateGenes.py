import os
from typing import List, Dict, Any, Set

try:
    from utils import (file_functions as ff,
                       iterable_functions as itf,
                       sequence_functions as sf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      iterable_functions as itf,
                                      sequence_functions as sf)

def identify_duplicate_gene(distinct_hashtable: str, 
                              schema_directory: str, 
                              output_directory: str, 
                              problematic_threshold: float,
                              no_cleanup: bool) -> None:
    """
    Identify problematic loci based on the presence of NIPHs and NIPHEMs.

    This function identifies problematic loci by analyzing the presence of NIPHs (Non-Identical Protein Hits)
    and NIPHEMs (Non-Identical Protein Hits with Multiple occurrences) in the provided CDS sequences. It calculates
    the proportion of problematic genomes for each locus and determines if it should be dropped based on a threshold.

    Parameters
    ----------
    distinct_hashtable : str
        Path to the distinct.hashtable file containing CDS sequences.
    schema_directory : str
        Path to the directory containing schema loci.
    output_directory : str
        Path to the directory where the output file will be saved.
    problematic_threshold : float
        Threshold for determining if a locus is problematic.

    Returns
    -------
    None
        The function writes the results to a TSV file and does not return any value.

    Notes
    -----
    - The function first decodes the CDS sequences from the provided file.
    - It then iterates through the schema loci, identifying NIPHs and NIPHEMs for each locus.
    - It calculates the proportion of problematic genomes for each locus and determines if it should be dropped.
    - Finally, it writes the results to a TSV file in the specified output directory.
    """
    print('Identifying problematic loci based on the presence of NIPHs and NIPHEMs...')
    print("Importing DNA hashes and the genomes where they are identified...")
    # Decode the CDS sequences from the provided file
    decoded_sequences_ids: Dict[str, List[Any]] = itf.decode_CDS_sequences_ids(distinct_hashtable)

    # Initialize dictionaries to store NIPHEMs and NIPHs in possible new loci
    niphems_in_loci: Dict[int, List[int]] = {}
    niphs_in_loci: Dict[int, Set[int]] = {}
    temp_niphs_in_loci: Dict[int, Dict[int, Set[int]]] = {}
    total_loci_genome_presence: Dict[int, int] = {}
    problematic_loci: Dict[int, float] = {}
    total_problematic_genomes_in_loci = {}

    # Fetch schema loci to a dictionary
    fastas_dict: Dict[int, Dict[int, str]] = sf.fetch_loci_to_dict(schema_directory)
    print("Calculating NIPHs and NIPHEMs proportion in loci...")
    # Iterate through each locus in the schema loci
    for loci_id, fastas in fastas_dict.items():
        niphems_in_loci.setdefault(loci_id, [])
        temp_niphs_in_loci.setdefault(loci_id, {})

        # Iterate through each allele in the locus
        for allele_id, sequence in fastas.items():
            hashed_seq: str = sf.seq_to_hash(str(sequence))
            allele_presence_in_genomes: List[int] = decoded_sequences_ids[hashed_seq][1:]
            
            # NIPHs
            temp_niphs_in_loci[loci_id].setdefault(allele_id, set(allele_presence_in_genomes))
            
            # NIPHEMs
            unique_elements: Set[int] = set(allele_presence_in_genomes)
            if len(unique_elements) != len(allele_presence_in_genomes):
                niphems_genomes: List[int] = itf.get_duplicates(allele_presence_in_genomes)
                niphems_in_loci.setdefault(loci_id, []).extend(niphems_genomes)

        # Get shared NIPHs in genomes
        niphs_in_genomes: Set[int] = set(itf.get_shared_elements(temp_niphs_in_loci[loci_id]))

        niphs_in_loci.setdefault(loci_id, niphs_in_genomes)
        
        # Get NIPHEMs in genomes
        niphems_in_genomes: Set[int] = set(niphems_in_loci[loci_id])
        
        # Calculate problematic genomes in possible new loci
        problematic_genomes_in_loci: Set[int] = niphs_in_genomes | niphems_in_genomes

        total_problematic_genomes_in_loci[loci_id] = problematic_genomes_in_loci
        
        # Calculate total genome presence
        total_loci_genome_presence[loci_id] = len(set(itf.flatten_list(temp_niphs_in_loci[loci_id].values())))
        
        # Calculate problematic proportion
        problematic_proportion: float = len(problematic_genomes_in_loci) / total_loci_genome_presence[loci_id]
        problematic_loci.setdefault(loci_id, problematic_proportion)

    # Write the groups that were removed due to the presence of NIPHs or NIPHEMs
    niphems_and_niphs_file: str = os.path.join(output_directory, 'niphems_and_niphs_groups.tsv')
    with open(niphems_and_niphs_file, 'w') as niphems_and_niphs:
        niphems_and_niphs.write('Loci\tcount_of_NIPHEMs\tcount_of_NIPHs\tTotal_unique_problematic_genomes\tTotal_genomes\tProportion_of_NIPHs_and_NIPHEMs\tOutcome\n')
        for loci_id, proportion in problematic_loci.items():
            outcome: str = 'Drop' if proportion >= problematic_threshold else 'Keep'
            niphems_and_niphs.write(f"{loci_id}\t{len(niphems_in_loci[loci_id])}\t{len(set(itf.get_shared_elements(temp_niphs_in_loci[loci_id])))}\t{len(total_problematic_genomes_in_loci[loci_id])}\t{total_loci_genome_presence[loci_id]}\t{proportion}\t{outcome}\n")
    # Clean up temporary files
    if not no_cleanup:
        print("\nCleaning up temporary files...")
        # Remove temporary files
        ff.cleanup(output_directory, [niphems_and_niphs_file])