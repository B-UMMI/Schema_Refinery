import hashlib
from Bio.Seq import Seq
from Bio import SeqIO

try:
    from RefineSchema.constants import DNA_BASES
except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema.constants import DNA_BASES

def check_str_alphabet(input_string, alphabet):
    """
    Verifies if all characters in the input string are present in the specified alphabet.

    Parameters
    ----------
    input_string : str
        The input sequence to be checked.
    alphabet : list or set
        The set of valid characters.

    Returns
    -------
    return : bool
        True if all characters in the input string are present in the alphabet, False otherwise.
    """

    alphabet_chars = set(alphabet)
    string_chars = set(input_string)

    diff = string_chars - alphabet_chars

    return len(diff) == 0

def check_str_multiple(input_string, number):
    """
    Verifies if the length of the input string is a multiple of the specified number.

    Parameters
    ----------
    input_string : str
        The input sequence to be checked.
    number : int
        The number to check for being a multiple of.

    Returns
    -------
    return : bool
        True if the length of the input string is a multiple of the specified number, False otherwise.
    """

    return (len(input_string) % number) == 0

def reverse_str(string):
    """
    Reverse the character order in the input string.

    Parameters
    ----------
    string : str
        The string to be reversed.

    Returns
    -------
    revstr : str
        The reverse of the input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """
    Determines the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.

    Returns
    -------
    reverse_complement : str
        The reverse complement of the input DNA sequence.
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence.upper())

    # determine complement strand
    complement_bases = []
    for base in bases:
        complement_bases.append(base_complement.get(base, base))

    complement_strand = ''.join(complement_bases)

    # reverse strand
    reverse_complement = reverse_str(complement_strand)

    return reverse_complement


def translate_sequence(dna_str, table_id, cds=True):
    """
    Translate a DNA sequence into a protein sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        DNA sequence to be translated.
    table_id : int
        Identifier of the translation table to use.
    cds : bool, optional
        Indicates whether the input sequence is a complete and valid CDS (default is True).

    Returns
    -------
    protseq : str
        Protein sequence obtained by translating the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    # sequences must be a complete and valid CDS
    protseq = Seq.translate(myseq_obj, table=table_id, cds=cds)

    return protseq


def determine_duplicated_seqs(sequences):
	"""Create mapping between sequences and sequence identifiers.

	Parameters
	----------
	sequences : dict
		Dictionary with sequence identifiers as keys and
		sequences as values.

	Returns
	-------
	equal_seqs : dict
		Dictionary with sequences as keys and sequence
		identifiers that are associated with each
		sequence as values.
	"""
	equal_seqs = {}
	for seqid, seq in sequences.items():
		# if protein sequence was already added as key
		if seq in equal_seqs:
			# append new protid
			equal_seqs[seq].append(seqid)
		# else add new protein sequence as key and protid
		# as value
		else:
			equal_seqs[seq] = [seqid]

	return equal_seqs

def determine_longest(seqids, sequences):
	"""Find the longest sequence in a set of sequences.

	Parameters
	----------
	seqids : list
		List with sequence identifiers.
	sequences : dict
		Dictionary with sequence identifiers as keys
		and sequences as values.

	Returns
	-------
	chosen : str
		Sequence identifier of the longest sequence.
	"""
	seqids_tups = [(seqid, sequences[seqid]) for seqid in seqids]
	sorted_tups = sorted(seqids_tups, key=lambda x: len(x[1]), reverse=True)
	chosen = sorted_tups[0][0]

	return chosen

def translate_dna_aux(dna_sequence, method, table_id, cds=True):
    """
    Translate a DNA sequence to a protein sequence using specified method.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence to be translated.
    method : str
        Translation method to use:
            - 'original': Translate the original sequence.
            - 'revcomp': Translate the reverse complement of the sequence.
            - 'rev': Translate the reverse of the sequence.
            - 'revrevcomp': Translate the reverse reverse complement of the sequence.
    cds : bool, optional
        Indicates whether the input sequence is a complete and valid CDS (default is True).

    Returns
    -------
    return : tuple or str
        If translation is successful, returns a tuple containing the translated protein sequence
        and the processed DNA sequence. If an error occurs, returns an error message.
    """

    myseq = dna_sequence

    try:
        if method == 'original':
            protseq = translate_sequence(myseq, table_id, cds=cds)
        elif method == 'revcomp':
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        elif method == 'rev':
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        elif method == 'revrevcomp':
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        else:
            return "Invalid method. Available methods: 'original', 'revcomp', 'rev', 'revrevcomp'"
        
        return [protseq, myseq]
    
    except Exception as e:
        return e

def translate_dna(dna_sequence, table_id, min_len, cds=True):
    """
    Translates DNA sequence to protein sequence.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence.
    table_id : str
        translation method.
    min_len : int
        minimum length for the sequence.
    cds : bool, optional
        if has cds

    Returns
    -------
    returns : str
        If failed verification steps, return string explaining what has failed.
    translated_seq : str
        Translated sequence.
    coding_strand : str
        DNA strand in correct orientation.
    """

    original_seq = dna_sequence.upper()
    exceptions = []
    coding_strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # Check if the sequence has ambiguous bases
    if not check_str_alphabet(original_seq, DNA_BASES):
        exceptions.append('ambiguous or invalid characters')

    # Check if sequence size is multiple of three
    if not check_str_multiple(original_seq, 3):
        exceptions.append('sequence length is not a multiple of 3')

    # Check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        exceptions.append(f'sequence shorter than {min_len} nucleotides')

    # Try to translate in 4 different orientations
    for method, coding_strand in zip(translating_methods, coding_strands):
        translated_seq = translate_dna_aux(original_seq, method, table_id, cds=cds)
        if isinstance(translated_seq, list):
            return [translated_seq, coding_strand]

        exceptions.append(f'{coding_strand}({translated_seq})')

    return ','.join(exceptions)


def read_fasta_file_iterator(file):
    """
    Reads fasta files and puts it into a dict.

    Parameters
    ----------
    file : str
        Path to the file

    Returns
    -------
    return : dict
        Returns iterator to read fasta file.
    """
    return SeqIO.parse(file, "fasta")

def read_fasta_file_dict(file):
    """
    Reads a FASTA file and returns a dictionary where the keys are sequence identifiers and the values are sequence records.

    Parameters
    ----------
    file : str
        Path to the FASTA file.

    Returns
    -------
    retun : dict
        A dictionary where the keys are sequence identifiers and the values are Bio.SeqRecord objects.
    """
    return SeqIO.to_dict(SeqIO.parse(file, "fasta"))

def seq_to_hash(seq):
    """
    Converts a sequence into hash.

    Parameters
    ----------
    seq : str
        Sequence

    Returns
    -------
    return : str
        Returns a hash string.
    """

    return hashlib.sha256(seq.encode('utf-8')).hexdigest()

def hash_sequences(file_path):
    """
    Hashes sequences in fasta file based on input file_path.

    Parameters
    ----------
    file_paths : str
        Contains file path to the fasta files.

    Returns
    -------
    hash_set : set
        Returns a list containing all of the sequences hashes present in the input files.
    """

    hash_set = set()
    for rec in read_fasta_file_iterator(file_path):
        hash_set.add(seq_to_hash(str(rec.seq)))

    return hash_set

def hash_sequence(input_string, hash_type='sha256'):
	"""Compute hash of an input string.

	Parameters
	----------
	input_string : str
		Input string to hash.
	hash_type : str
		Hash type/function that will be used to compute the
		hash (any of the hash functions available in the
		hashlib module).

	Returns
	-------
	hashed_string : str
		String representation of the HASH object
		in hexadecimal digits.
	"""
	# get hash function object from hashlib
	hashing_function = getattr(hashlib, hash_type)

	# default encoding is UTF-8
	hashed_string = hashing_function(input_string.encode()).hexdigest()

	return hashed_string

def sequence_lengths(fasta_file, hashed=False):
	"""Determine length of sequences in a FASTA file.

	Read Fasta file and create dictionary with mapping
	between sequence identifiers and sequence lengths.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.
	hashed : bool
		If False, sequence headers are used as
		keys. If True, sequence hashes will be
		used as keys.

	Returns
	-------
	lengths : dict
		Dictionary with sequence identifiers as keys and
		sequence lengths as values.
	"""
	records = sequence_generator(fasta_file)
	if hashed is False:
		lengths = {rec.id: len(rec.seq) for rec in records}
	else:
		lengths = {hash_sequence(str(rec.seq)): len(rec.seq) for rec in records}

	return lengths

def fasta_stats(fasta_file):
	"""Determine the number of sequences in a FASTA file and length stats.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.

	Returns
	-------
	fasta_file : str
		Path to the FASTA file.
	total_seqs: int
		Total number of records in the FASTA file.
	mean_length: float
		Mean sequence length.
	"""
	seq_lengths = sequence_lengths(fasta_file)
	min_length = min(seq_lengths.values())
	max_length = max(seq_lengths.values())
	mean_length = sum(seq_lengths.values())/len(seq_lengths)
	total_seqs = len(seq_lengths)

	return [fasta_file, total_seqs, min_length, max_length, mean_length]

def import_sequences(input_file):
	"""Import sequences from a FASTA file.

	Parameters
	----------
	input_file : str
		Path to a FASTA file.

	Returns
	-------
	records_dict : dict
		Dictionary with sequence identifiers as keys and
		sequences as values.
	"""
	records = sequence_generator(input_file)
	# Only want record identifier and sequence, no need to use SeqIO.to_dict
	records_dict = {rec.id: str(rec.seq.upper()) for rec in records}

	return records_dict

def translate_seq_deduplicate(seq_dict, path_to_write, untras_path, min_len, count_seq,
                              translation_table, deduplicate = True):
    """
    Translates the DNA sequence to protein and verifies if that protein is alredy
    present in the dict, thus ensuring that the dict contains deduplicated sequences,
    it writes the sequences to a FASTA files and return the dict.
    
    Parameters
    ----------
    seq_dict : dict
        Dict that contains sequence ID as key and the sequence as value.
    path_to_write : str
        Path to the file to create and write.
    untras_path : str or None
        Path to write untranslated sequences, if there is no need to write use None
    min_len : int
        minimum length for the sequence.
    count_seq : bool
        If there is need to print into stdout the number of processed sequences.
    deduplicate : bool, optional
        If the process of sequence deduplication is needed.
    
    Returns
    -------
    translation_dict : dict
        Dict that contais sequence ID as key and translated sequence as value.
    protein_hashes : dict
        Dict that contais sequence hash as key and sequences IDs as values.
    untras_seq : dict
        Dict that contains untranslated id as key and error as value.
    """
    translation_dict = {}
    protein_hashes = {}
    untras_seq = {}
    if count_seq:
        total = len(seq_dict)
        
    with open(path_to_write, 'w+') as translation:
        for i, (id_s, sequence) in enumerate(seq_dict.items(),1):
                
            # Translate
            protein_translation = translate_dna(str(sequence),
                                                    translation_table,
                                                    min_len,
                                                    True)
            # Verify if was translated
            if type(protein_translation) == list:
                protein_translation = str(protein_translation[0][0])
                if count_seq:
                    print(f"\rTranslated {i}/{total} CDS", end='', flush=True)
            else:
                if count_seq:
                    print(f"\rFailed to translate {id_s}", end='', flush=True)
                untras_seq.setdefault(id_s, protein_translation)
                continue
            
            if deduplicate:
                # Hash the sequence
                prot_hash = seq_to_hash(protein_translation)
                # Find unique proteins
                if prot_hash not in protein_hashes:
                    protein_hashes[prot_hash] = [id_s]
                    translation_dict[id_s] = protein_translation
                    translation.write(f'>{id_s}\n{protein_translation}\n')
                # Remember CDS with that protein hash for future
                else:
                    protein_hashes[prot_hash].append(id_s)
            else:
                translation_dict[id_s] = protein_translation
                translation.write(f'>{id_s}\n{protein_translation}\n')
    if untras_seq and untras_path:
        with open(untras_path, 'w+') as untras_file:
            for id_s, exceptions in untras_seq.items():
                untras_file.write(">{}\n{}\n".format(id_s,'\n'.join(exceptions)))

    return translation_dict, protein_hashes, untras_seq

def fetch_fasta_dict(file_path, count_seq):
    """
    Fetches FASTAs from a file and adds the to a dict

    Parameters
    ----------
    file_path : str
        Path to the file with FASTAs
    count_seq : bool
        If count the number of processed sequences inside the file

    Returns
    -------
    fasta_dict : dict
        Returns dict with key as fasta header and value as fasta sequence.
    """
    
    fasta_dict = {}
    i = 1
    # Read FASTA files
    for rec in read_fasta_file_iterator(file_path):
        if count_seq:
            print(f"\rProcessed {i} CDS", end='', flush=True)
            i += 1
        fasta_dict[rec.id] = rec.seq

    return fasta_dict

def deduplicate_fasta_dict(fasta_dict):
    """
    Deduplicates a dictionary of FASTA sequences. The deduplication is based on the SHA256 hash of the sequences.

    Parameters
    ----------
    fasta_dict : dict
        A dictionary where the keys are sequence identifiers and the values are sequences.

    Returns
    -------
    fasta_dict : dict
        A dictionary where the keys are sequence identifiers and the values are sequences. Sequences that were duplicated in the input dictionary are removed.
    """

    deduplicated_list = []
    for key, sequence in fasta_dict.items():
        # Create a hash of the sequence
        sequence_hash = hashlib.sha256(sequence.encode('utf-8')).hexdigest()
        # If the hash is not already a key in the deduplicated_dict, add the sequence
        if sequence_hash not in deduplicated_list:
            deduplicated_list.append(sequence_hash)
        else:
            del fasta_dict[key]

    return fasta_dict

def sequence_generator(input_file):
	"""Create a SeqRecord iterator.

	Parameters
	----------
	input_file : str
		Path to a Fasta file.

	Returns
	-------
	records : Bio.SeqIO.FastaIO.FastaIterator
		SeqRecord iterator.
	"""
	# Useful to create the generator
	# Need to exhaust the generator to avoid high memory usage
	records = SeqIO.parse(input_file, 'fasta')

	return records

def fasta_str_record(record_template, record_data):
	"""Create the string representation of a FASTA record.

	Parameters
	----------
	record_template : str
		String template to construct the FASTA record.
	record_data : list
		List with the elements to add to the string.

	Returns
	-------
	record : str
		String representation of the FASTA record.
	"""
	record = record_template.format(*record_data)

	return record


def fasta_lines(template, records_data):
	"""Create a list with FASTA records.

	Parameters
	----------
	template : str
		String template to construct the FASTA record.
	records_data : list
		A list with one sublist per FASTA record.
		Each sublist contains the elements to insert
		inside the template placeholders.

	Returns
	-------
	seqs_lines : list
		A list with strings representing FASTA records.
	"""
	seqs_lines = [fasta_str_record(template, arg) for arg in records_data]

	return seqs_lines