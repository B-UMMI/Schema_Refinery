import os

try:
    from utils import (file_functions as ff,
                       iterable_functions as itf,
                       sequence_functions as sf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      iterable_functions as itf,
                                      sequence_functions as sf)

def dropped_loci_to_file(schema_loci, dropped, results_output):
    dropped_file = os.path.join(results_output, "dropped.tsv")
    
    with open(dropped_file, 'w') as d:
        d.write("ID\tReason\tWhere from\n")
        for drop in dropped:
            if drop in schema_loci:
                dropped_from = 'from_schema'
            else:
                dropped_from = 'from_possible_new_loci'
            d.write(f"{drop}\t{'Dropped_due_to_cluster_frequency_filtering'}\t{dropped_from}\n")

def process_new_loci(schema_folder, allelecall_directory, constants, processing_mode, results_output):
    schema = {fastafile: os.path.join(schema_folder, fastafile) for fastafile in os.listdir(schema_folder) if fastafile.endswith('.fasta')}
    schema_short_dir = os.path.join(schema_folder, 'short')
    schema_short = {fastafile: os.path.join(schema_short_dir, fastafile) for fastafile in os.listdir(schema_short_dir) if fastafile.endswith('.fasta')}
    master_file_path = os.path.join(results_output, 'master.fasta')
    
    possible_new_loci_translation_folder = os.path.join(results_output, 'schema_translation_folder')
    ff.create_directory(possible_new_loci_translation_folder)

    to_blast_paths = schema if processing_mode.split('_')[0] == 'alleles' else schema_short
    to_run_against = schema_short if processing_mode.split('_')[-1] == 'rep' else schema

    all_alleles = {}
    alleles = {}
    translation_dict = {}
    frequency_in_genomes = {}
    temp_frequency_in_genomes = {}
    cds_present = os.path.join(allelecall_directory, "temp", "2_cds_preprocess/cds_deduplication/distinct.hashtable")
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)
    # Alleles to run
    for loci in to_blast_paths.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        alleles.setdefault(loci_id, {})
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():
            alleles.setdefault(loci_id, {}).update({allele_id: str(sequence)})

    # Write master file to run against
    for loci in to_run_against.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        alleles.setdefault(loci_id, {})
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():
            alleles.setdefault(loci_id, {}).update({allele_id: str(sequence)})
            # Write to master file
            write_type = 'a' if os.path.exists(master_file_path) else 'w'
            with open(master_file_path, write_type) as master_file:
                master_file.write(f">{allele_id}\n{str(sequence)}\n")

    # Count loci presence and translate all of the alleles.
    for loci in schema.values():
        loci_id = ff.file_basename(loci).split('.')[0]
        all_alleles.setdefault(loci_id, [])
        fasta_dict = sf.fetch_fasta_dict(loci, False)
        for allele_id, sequence in fasta_dict.items():  
            all_alleles[loci_id].append(allele_id)
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus the first (subtract the first CDS genome).
            if hashed_seq in decoded_sequences_ids:
                #Count frequency of only presence, do not include the total cds in the genomes.
                temp_frequency_in_genomes.setdefault(loci_id, []).append(len(set(decoded_sequences_ids[hashed_seq][1:])))

        frequency_in_genomes.setdefault(loci_id, sum(temp_frequency_in_genomes[loci_id]))

        trans_path_file = os.path.join(possible_new_loci_translation_folder, f"{loci_id}.fasta")

        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        trans_path_file,
                                                        None,
                                                        constants[5],
                                                        False,
                                                        constants[6],
                                                        False)
        
        translation_dict.update(trans_dict)
                
    return alleles, master_file_path, translation_dict, frequency_in_genomes, to_blast_paths, all_alleles