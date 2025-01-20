import networkx as nx
from collections import Counter
from typing import Dict, List, Set, Tuple, Union

try:
    from utils import (kmers_functions as kf,
                       iterable_functions as itf)
except:
    from SchemaRefinery. utils import (kmers_functions as kf,
                                       iterable_functions as itf)

def select_representatives(kmers: Union[List[Tuple[str, int]], Set[Tuple[str, int]]], 
                           reps_groups: Dict[str, List[str]], 
                           clustering_sim: float, 
                           clustering_cov: float,
                           prot_len_dict: Dict[str, int], 
                           protid: str, 
                           window_size: int, 
                           size_threshold: float) -> List[Tuple[str, float, float]]:
    """
    Determine the clusters a sequence can be added to.

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
    clustering_cov : float
        Coverage threshold to cluster a sequence into
        a cluster.
    prot_len_dict : dict
        Contains the length of each CDS.
    protid : str
        Identifier of the protein.
    window_size : int
        Window size used to determine minimizers.
    size_threshold : float
        Size threshold for clustering.

    Returns
    -------
    selected_reps : list
        List with a tuple per cluster/representative
        that the sequence can be added to. Each tuple
        has the identifier of the cluster representative,
        the decimal proportion of shared distinct kmers,
        and the coverage.
    """
    # {start_pos: [[rep_cds1], [rep_cds2]]}
    current_reps: Dict[int, List[str]] = {k[1]: reps_groups[k[0]] for k in kmers if k[0] in reps_groups}
    
    all_rep_pos: Dict[str, List[int]] = {}
    for k, v in current_reps.items():
        for value in v:
            all_rep_pos.setdefault(value, []).append(k)

    # count number of kmer hits per representative
    counts: Counter[str] = Counter(itf.flatten_list(current_reps.values()))
    
    # selects reps_loci based on number of kmer hits/total number of kmers
    selected_reps: List[Tuple[str, float]] = [(k, v / len(kmers))
                                              for k, v in counts.items()
                                              if v / len(kmers) >= clustering_sim]
            
    # sort by identifier and then by similarity to always get same order
    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)
    
    selected_reps_coverage: List[str] = [rep[0] for rep in selected_reps]
    
    rep_coverage_all: Dict[str, float] = {}
    # Calculates the coverage of query kmers over rep prot sequence
    for rep in selected_reps_coverage:
        # get the pos of kmer if that kmer hit against the rep kmers
        rep_coverage: List[int] = sorted(all_rep_pos[rep], key=lambda x: x)
        # calculate coverage
        rep_coverage_all[rep] = kf.kmer_coverage(rep_coverage, window_size) / prot_len_dict[protid]
    
    selected_reps = [(*representative, rep_coverage_all[representative[0]]) 
                     for representative in selected_reps 
                     if rep_coverage_all[representative[0]] >= clustering_cov]
    
    if size_threshold:
        max_limit: float = prot_len_dict[protid] + prot_len_dict[protid] * size_threshold
        min_limit: float = prot_len_dict[protid] - prot_len_dict[protid] * size_threshold
    
        selected_reps = [values for values in selected_reps if max_limit >= prot_len_dict[values[0]] >= min_limit]
            
    return selected_reps


def minimizer_clustering(sorted_sequences: Dict[str, str], 
                         word_size: int, 
                         window_size: int, 
                         position: bool,
                         offset: int, 
                         clusters: Dict[str, List[Tuple[str, float, int, int, int, float]]], 
                         reps_sequences: Dict[str, str], 
                         reps_groups: Dict[str, List[str]],
                         seq_num_cluster: int, 
                         clustering_sim: float, 
                         clustering_cov: float, 
                         grow: bool,
                         size_threshold: float) -> Tuple[Dict[str, List[Tuple[str, float, int, int, int, float]]], 
                                                        Dict[str, str], 
                                                        Dict[str, List[str]], 
                                                        Dict[str, int]]:
    """
    Cluster sequences based on shared distinct minimizers.

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
    clustering_cov : float
        Coverage threshold to cluster a sequence into
        a cluster.
    grow : bool
        If clusters are fixed or create new clusters.
    size_threshold : float
        Size threshold for clustering.
        
    Returns
    -------
    Tuple[Dict[str, List[Tuple[str, float, int, int, int, float]]], Dict[str, str], Dict[str, List[str]], Dict[str, int]]
        A tuple containing the clusters, representative sequences, representative groups, and protein length dictionary.
    """
    # get len of all the proteins
    prot_len_dict: Dict[str, int] = {protid: len(protein) for protid, protein in sorted_sequences.items()}

    for protid, protein in sorted_sequences.items():
        minimizers: List[Tuple[str, int]] = kf.determine_minimizers(protein, window_size,
                                                                    word_size, offset=offset,
                                                                    position=position, guarantee_tip=True)
        
        distinct_minimizers: Set[str] = {minimizer[0] for minimizer in minimizers}
        
        selected_reps: List[Tuple[str, float, float]] = select_representatives(minimizers,
                                                                               reps_groups,
                                                                               clustering_sim, clustering_cov,
                                                                               prot_len_dict, protid, 
                                                                               window_size,
                                                                               size_threshold)
        
        top: int = (len(selected_reps)
                    if len(selected_reps) < seq_num_cluster
                    else seq_num_cluster)

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for i in range(0, top):
                clusters[selected_reps[i][0]].append((protid,
                                                      selected_reps[i][1],
                                                      len(protein),
                                                      len(minimizers),
                                                      len(distinct_minimizers),
                                                      selected_reps[i][2]))
        else:
            if grow is True:
                for k in minimizers:
                    reps_groups.setdefault(k[0], set()).add(protid)

                clusters[protid] = [(protid, 1.0, 
                                     len(protein),
                                     len(minimizers), 
                                     len(distinct_minimizers),
                                     kf.kmer_coverage(
                                         sorted([i[1] for i in minimizers])
                                        ,window_size)/len(protein))]
                
                reps_sequences[protid] = protein

    return clusters, reps_sequences, reps_groups, prot_len_dict


def cluster_by_ids(list_of_ids: List[List[str]]) -> List[List[str]]:
    """
    Based on list of list containing pairs or more ids, merges them into larger
    lists if there are ids common between these sublists. This function maintains
    the order of appearance in list_of_ids.
    
    Parameters
    ----------
    list_of_ids : list
        List containing lists with ids.
    
    Returns
    -------
    connected : list
        List that contains sublists of ids merged from the initial list.
    """
    G = nx.Graph()
    G.add_edges_from(list_of_ids)

    connected_components: List[Set[str]] = list(nx.connected_components(G))

    # Create a mapping from node to its first appearance index in the input list
    node_order: Dict[str, int] = {}
    for index, node in enumerate(itf.flatten_list(list_of_ids)):
        if node not in node_order:
            node_order[node] = index

    # Sort connected components based on the first appearance of any node in the component
    connected_components.sort(key=lambda component: min(node_order[node] for node in component))

    # Convert sets to lists to maintain order
    connected: List[List[str]] = [sorted(list(component), key=lambda node: node_order[node]) for component in connected_components]

    return connected

