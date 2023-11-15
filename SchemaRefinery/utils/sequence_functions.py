from Bio.Seq import Seq

try:
    from RefineSchema.constants import DNA_BASES
except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema.constants import DNA_BASES

def check_str_alphabet(input_string, alphabet):
    """
    Verifies if input string has valid characters.

    Parameters
    ----------
    input_string : str
        Input sequence.
    alphabet : list
        Contains list or set of letters for comparison.

    Returns
    -------
    return : bool
        if is or not multiple.
    """

    alphabet_chars = set(alphabet)
    string_chars = set(input_string)

    diff = string_chars - alphabet_chars

    return len(diff) == 0

def check_str_multiple(input_string, number):
    """
    Verifies is input string length is multiple of input number.
    
    Parameters
    ----------
    input_string : str
        Input sequence.
    number : int
        number to be multiple of.

    Returns
    -------
    return : bool
        if is or not multiple.    
    """

    return (len(input_string) % number) == 0

def reverse_str(string):
    """ Reverse character order in input string.

    Parameters
    ----------
    string : str
        String to be reversed.

    Returns
    -------
    revstr : str
        Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.

    Returns
    -------
    reverse_complement : str
        The reverse complement of the input DNA
        sequence.
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
    """ Translate a DNA sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        DNA sequence.
    table_id : int
        Translation table identifier.

    Returns
    -------
    protseq : str
        Protein sequence created by translating
        the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    # sequences must be a complete and valid CDS
    protseq = Seq.translate(myseq_obj, table=table_id, cds=cds)

    return protseq

def translate_dna_aux(dna_sequence, method, table_id, cds=True):
    """
    Translates DNA sequence to protein sequence.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence.
    table_id : string
        translation method.
    method : str
        Translation orientatio.
    cds : bool
        if has cds.

    Returns
    -------
    argh : exception
        returns exception error.
    protseq : str
        Translated sequence.
    myseq : str
        Original DNA sequence.
    """

    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id, cds=cds)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq, DNA_BASES)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq, DNA_BASES)
            protseq = translate_sequence(myseq, table_id, cds=cds)
        except Exception as argh:
            return argh

    return [protseq, myseq]

def translate_dna(dna_sequence, table_id, min_len, cds=True):
    """
    Translates DNA sequence to protein sequence.

    Parameters
    ----------
    dna_sequence : str
        DNA sequence.
    table_id : string
        translation method.
    min_len : int
        minimum length for the sequence.
    cds : bool
        if has cds

    Returns
    -------
    returns : str
        If failed verification steps, return string explaining what has failed.
    translated_seq : str
        Translated sequence.
    coding_strand : str
        DNA strand in correct orientation.
    exception_str : str
        String containing all the exceptions found.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the sequence has ambiguous bases
    valid_dna = check_str_alphabet(original_seq, DNA_BASES)
    if valid_dna is not True:
        return 'ambiguous or invalid characters'

    # check if sequence size is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return 'sequence length is not a multiple of 3'

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        translated_seq = translate_dna_aux(original_seq, translating_methods[i], table_id, cds=cds)
        if not isinstance(translated_seq, list):
            exception_collector.append('{0}({1})'.format(strands[i],
                                                         translated_seq.args[0]))

        i += 1
        if i == len(strands) or isinstance(translated_seq, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(translated_seq, list):
        return [translated_seq, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str
