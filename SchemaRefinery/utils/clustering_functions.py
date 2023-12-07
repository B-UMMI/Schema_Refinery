import networkx as nx
from collections import Counter

try:
    from utils.kmers_functions import determine_minimizers
    from utils.list_functions import flatten_list
except:
    from SchemaRefinery.utils.kmers_functions import determine_minimizers
    from SchemaRefinery.utils.list_functions import flatten_list

def select_representatives(kmers, reps_groups, clustering_sim, clustering_cov,
                           prot_len_dict, protid):
    """Determine the clusters a sequence can be added to.

    Determines the set of clusters that a sequence can be
    added to based on the decimal proportion of shared
    distinct k-mers.

    Parameters
    ----------
    kmers : list or set
        Set of k-mers determined by decomposing a single
        sequence.
    reps_groups : dict
        Dictionary with k-mers as keys and sequence
        identifiers of sequences that contain each
        k-mer as values.
    clustering_sim : float
        Sequences are added to clusters if they
        share a minimum decimal proportion of
        distinct k-mers with a cluster representative.
    prot_len_dict : dict
        Containins the length of each CDS

    Returns
    -------
    selected_reps : list
        List with a tuple per cluster/representative
        that the sequence can be added to. Each tuple
        has the identifier of the cluster representative
        and the decimal proportion of shared distinct
        kmers.
    """
    # {start_pos: [[rep_cds1], [rep_cds2]]}
    current_reps = {k[1] : reps_groups[k[0]] for k in kmers if k[0] in reps_groups}
    
    # count number of kmer hits per representative
    counts = Counter(flatten_list(current_reps.values()))
    #selects reps_loci based on number of kmer hits/total number of kmers
    selected_reps = [(k, v/len(kmers))
                     for k, v in counts.items()
                     if v/len(kmers) >= clustering_sim]
            

    # sort by identifier and then by similarity to always get same order
    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)

    selected_reps_coverage = [rep[0] for rep in selected_reps]
    
    rep_coverage_all = {}
    print(protid)
    #Calculates the coverage of query kmers over rep prot sequence
    for rep in selected_reps_coverage:
        #get the pos of kmer if that kmer hit against the rep kmers
        rep_coverage = sorted([k for k, v in current_reps.items()
                                    if rep in v], key=lambda x: x)
        #calculate coverage
        rep_coverage_all[rep] = kmer_coverage(rep_coverage)/prot_len_dict[rep]
        print("\t",rep,kmer_coverage(rep_coverage)/prot_len_dict[rep])
    
    selected_reps = [(*rep,rep_coverage_all[rep[0]]) 
                     for rep in selected_reps 
                     if rep_coverage_all[rep[0]] >= clustering_cov]
        
            
    return selected_reps

def kmer_coverage(position):
    
    length_query = 0
    size = 0
    len_positions = len(position)-1
    for i, pos in enumerate(position):
        if i == 0:
            start = pos
            end = pos+4
        elif pos <= end:
            end = pos+4
        else:
            size += end-start+1
            start = pos
            end = pos+4
            
        if i == len_positions:
            end = pos+4
            size += end-start+1
    return size
    
def minimizer_clustering(sorted_sequences, word_size, window_size, position,
                         offset, clusters, reps_sequences, reps_groups,
                         seq_num_cluster, clustering_sim, clustering_cov, grow):
    """Cluster sequences based on shared distinct minimizers.

    Parameters
    ----------
    sorted_sequences : dict
        Dictionary with sequence identifiers as keys and
        sequences as values. Sorted by decreasing sequence
        length.
    word_size : int
        Value k for the kmer size.
    window_size : int
        Window size used to determine minimizers.
    position : bool
        If minimizer sequence position should be stored.
    offset : int
        Value to indicate offset of consecutive kmers.
    clusters : dict
        Dictionary with the identifiers of sequences
        that are cluster representatives as keys and
        a list with tuples as values. Each tuple has
        the identifier of a sequence that was added to
        the cluster, the decimal proportion of shared
        distinct minimizers and the length of the clustered
        sequence. This dictionary should be empty at the
        start of the clustering step during the CreateSchema
        process. For the AlleleCall process, the dictionary
        should contain the identifiers of the loci
        representatives.
    reps_sequences : dict
        Dictionary with the identifiers of sequences
        that are cluster representatives as keys and
        their sequences as values. This dictionary should
        be empty at the start of the clustering step
        during the CreateSchema process. For the AlleleCall
        process, the dictionary should contain the
        identifiers of the loci representatives as keys
        and their sequences as values.
    reps_groups : dict
        Dictionary with kmers as keys and a list with
        identifiers of sequences that contain that kmer
        as values. This dictionary should be empty at the
        start of the clustering step during the CreateSchema
        process. For the AlleleCall process, the dictionary
        should contain all kmers contained in the set of
        the schema's representatives sequences.
    seq_num_cluster : int
        Maximum number of clusters that a sequence can be
        added to.
    clustering_sim : float
        Similarity threshold to cluster a sequence into
        a cluster.

    Returns
    -------
    A list with the following elements:
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the decimal proportion of shared
            distinct minimizers and the length of the clustered
            sequence.
        reps_sequences : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            their sequences as values.
        reps_groups : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
    """
    # several = {}
    
    prot_len_dict = {protid: len(protein) for protid, protein 
                     in sorted_sequences.items()}
    
    for protid, protein in sorted_sequences.items():
        minimizers = determine_minimizers(protein, window_size,
                                             word_size, offset=offset,
                                             position=position)
        ##remove this because I am using start pos
        distinct_minimizers = set(minimizers)
        
        selected_reps = select_representatives(distinct_minimizers,
                                               reps_groups,
                                               clustering_sim, clustering_cov,
                                               prot_len_dict, protid)
        
        top = (len(selected_reps)
               if len(selected_reps) < seq_num_cluster
               else seq_num_cluster)

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for i in range(0, top):
                clusters[selected_reps[i][0]].append((protid,
                                                      selected_reps[i][1],
                                                      len(protein),
                                                      len(minimizers),
                                                      len(distinct_minimizers)))
####
            # if len(selected_reps) > 1:
            #     for i in range(0, top):
            #         several.setdefault(protid, []).append(selected_reps[i][0])
####
        else:
            if grow is True:
                for k in distinct_minimizers:
                    reps_groups.setdefault(k[0], []).append(protid)

                clusters[protid] = [(protid, 1.0, len(protein),
                                    len(minimizers), len(distinct_minimizers))]
                reps_sequences[protid] = protein

    # print(several)

    return [clusters, reps_sequences, reps_groups]

def cluster_based_on_ids(processed_representatives_dict):
    """
    Employs networkx to cluster loci based on their connection through ids. e.g x matches with y
    y matches with z, all three are put inside the same cluster x,y and z.

    Parameters
    ----------
    processed_representatives_dict : dict
        dict that contains all the info necessary to render the graphs.

    Returns
    -------
    connected : list
        list containing lists of the ids of the clustered loci.
    """

    pairs_list = set()

    for key in processed_representatives_dict.keys():

        pairs_list.add(tuple([locus.split("_")[0] for locus in key.split(";")]))

    G = nx.Graph()
    G.add_edges_from(pairs_list)

    connected = nx.connected_components(G)

    return connected