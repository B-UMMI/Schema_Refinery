import subprocess


def make_blast_db(input_fasta, output_path, db_type):
    """
    Create a BLAST database.

    Parameters
    ----------
    input_fasta : str
        Path to a Fasta file.
    output_path : str
        Path to the output BLAST database.
    db_type : str
        Type of the database, nucleotide (nuc) or
        protein (prot).

    Returns
    -------
        returns : No return
            Creates blast database
    """
    blastdb_cmd = ['makeblastdb', '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stdout, stderr = makedb_cmd.communicate()
    print(stdout, stderr)

    makedb_cmd.wait()


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None):
    """
    Execute BLAST.

    Parameters
    ----------
    blast_path : str
        Path to the BLAST executable.
    blast_db : str
        Path to the BLAST database.
    fasta_file : str
        Path to the Fasta file that contains the sequences
        to align against the database.
    blast_output : str
        Path to the output file.
    max_hsps : int
        Maximum number of High-Scoring Pairs.
    threads : int
        Number of threads passed to BLAST.
    ids_file : path
        Path to a file with the identifiers of the sequences
        to align against. Used to specify the database sequences
        we want to align against.
    blast_task : str
        BLAST task. Allows to set default parameters for a specific
        type of search.
    max_targets : int
        Maximum number of targets sequences to align against.

    Returns
    -------
    stderr : list
        List with the warnings/errors reported by BLAST.
    """
    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])
    if blast_task is not None:
        blast_args.extend(['-task', blast_task])
    if max_targets is not None:
        blast_args.extend(['-max_target_seqs', str(max_targets)])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    return stderr

def run_blast_with_args_only(blast_args):
    """
    Runs BLAST based on input arguments.

    Parameters
    ----------
    blast_args : list
        Contains list with arguments used in subprocess.

    Returns
    -------
    returns : No return
        Generates BLAST files at output directory
    """
    blast_proc = subprocess.Popen(blast_args,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()
    if len(stderr) > 0:
        print(stderr)
