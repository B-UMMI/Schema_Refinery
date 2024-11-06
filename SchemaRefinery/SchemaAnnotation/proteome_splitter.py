import os
import gzip
import pickle
from typing import Dict, Tuple, List
from Bio import SeqIO

try:
    from utils.file_functions import create_directory
except ModuleNotFoundError:
    from SchemaRefinery.utils.file_functions import create_directory

def proteome_splitter(proteomes_directory: str, output_directory: str) -> Tuple[str, str, str, Dict[str, List[str]]]:
    """
    Split proteome records into Swiss-Prot and TrEMBL records and save them to separate files.

    Parameters
    ----------
    proteomes_directory : str
        Path to the directory containing proteome files.
    output_directory : str
        Path to the directory where output files will be saved.

    Returns
    -------
    Tuple[str, str, str]
        Paths to the TrEMBL records file, Swiss-Prot records file, and descriptions file.
    """
    
    split_directory: str = os.path.join(output_directory, 'split_proteomes')
    create_directory(split_directory)
    # List proteomes in proteome directory
    proteomes: List[str] = [os.path.join(proteomes_directory, f)
                            for f in os.listdir(proteomes_directory)]

    # Initialize dictionaries to store Swiss-Prot, TrEMBL records, and descriptions
    swiss_records: Dict[str, str] = {}
    trembl_records: Dict[str, str] = {}
    descriptions: Dict[str, str] = {}
    proteome_file_ids: Dict[str, List[str]] = {}
    # Process each proteome file
    for file in proteomes:
        file_name: str = os.path.basename(file).split('.')[0]
        with gzip.open(file, 'rt') as gzfasta:
            for rec in SeqIO.parse(gzfasta, 'fasta'):
                recid: str = rec.id
                prot: str = str(rec.seq)
                desc: str = rec.description
                # Add record ID to appropriate proteome file
                proteome_file_ids.setdefault(file_name, []).append(recid)
                
                if recid.startswith('tr'):
                    trembl_records[recid] = prot
                elif recid.startswith('sp'):
                    swiss_records[recid] = prot
                descriptions[recid] = desc

    # Save TrEMBL records to FASTA file
    tr_filename: str = 'trembl_prots.fasta'
    tr_file_path: str = os.path.join(split_directory, tr_filename)
    with open(tr_file_path, 'w') as tout:
        records: List[str] = ['>{0}\n{1}'.format(k, v) for k, v in trembl_records.items()]
        rectext: str = '\n'.join(records)
        tout.write(rectext + '\n')

    # Save Swiss-Prot records to FASTA file
    sp_filename: str = 'swiss_prots.fasta'
    sp_file_path: str = os.path.join(split_directory, sp_filename)
    with open(sp_file_path, 'w') as tout:
        records = ['>{0}\n{1}'.format(k, v) for k, v in swiss_records.items()]
        rectext = '\n'.join(records)
        tout.write(rectext + '\n')

    # Save descriptions with Pickle
    descriptions_filename: str = 'prots_descriptions'
    descriptions_file_path: str = os.path.join(split_directory, descriptions_filename)
    with open(descriptions_file_path, 'wb') as dout:
        pickle.dump(descriptions, dout)

    return tr_file_path, sp_file_path, descriptions_file_path, proteome_file_ids