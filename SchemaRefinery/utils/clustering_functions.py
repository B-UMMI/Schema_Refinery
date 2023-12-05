import networkx as nx
from collections import Counter

try:
    from utils.kmers_functions import determine_minimizers
    from utils.list_functions import flatten_list
except:
    from SchemaRefinery.utils.kmers_functions import determine_minimizers
    from SchemaRefinery.utils.list_functions import flatten_list

def select_representatives(kmers, reps_groups, clustering_sim):
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

    Returns
    -------
    selected_reps : list
        List with a tuple per cluster/representative
        that the sequence can be added to. Each tuple
        has the identifier of the cluster representative
        and the decimal proportion of shared distinct
        kmers.
    """
    current_reps = [reps_groups[k] for k in kmers if k in reps_groups]
    current_reps = flatten_list(current_reps)
  
    # count number of kmer hits per representative
    counts = Counter(current_reps)
    selected_reps = [(k, v/len(kmers))
                     for k, v in counts.items()
                     if v/len(kmers) >= clustering_sim]

    # sort by identifier and then by similarity to always get same order
    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)
            
    return selected_reps

def minimizer_clustering(sorted_sequences, word_size, window_size, position,
                         offset, clusters, reps_sequences, reps_groups,
                         seq_num_cluster, clustering_sim, grow):
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
    protein_length = {}
    for protid, protein in sorted_sequences.items():
        protein_length[protid] = len(protein)
        
        minimizers = determine_minimizers(protein, window_size,
                                             word_size, offset=offset,
                                             position=position)
        coverage_kmers = []
        for i, kmer in enumerate(minimizers):
            if i == 1:
                start = kmer[1]
                shift = kmer[1]
            elif kmer[1] >= shift and kmer[1] <= shift+5:
                shift = kmer[1]
            else:
                coverage_kmers.append([start, shift])
                
                
        distinct_minimizers = set(minimizers)

        selected_reps = select_representatives(distinct_minimizers,
                                               reps_groups,
                                               clustering_sim)
        
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
                    reps_groups.setdefault(k, []).append(protid)

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