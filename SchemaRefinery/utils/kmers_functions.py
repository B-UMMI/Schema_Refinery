from typing import List, Optional, Union, Tuple

def string_kmerizer(input_string: str, k_value: int, offset: int = 1, position: bool = False) -> List[Union[str, Tuple[str, int]]]:
    """
    Decompose a string into k-mers.

    Parameters
    ----------
    input_string : str
        String to divide into k-mers.
    k_value : int
        Value for the size of k-mers.
    offset : int, optional
        Value to indicate offset of consecutive k-mers.
    position : bool, optional
        If the start position of the k-mers in the string
        should be stored.

    Returns
    -------
    kmers : list
        List that contains the k-mers determined for the
        input string. The list will contain strings if
        it is not specified that positions should be
        stored and tuples of k-mer and start position
        if the position is stored.
    """
    kmers: List[Union[str, Tuple[str, int]]]
    if position is False:
        kmers = [input_string[i:i+k_value]
                 for i in range(0, len(input_string)-k_value+1, offset)]
    elif position is True:
        kmers = [(input_string[i:i+k_value], i)
                 for i in range(0, len(input_string)-k_value+1, offset)]

    return kmers


def determine_minimizers(input_string: str, adjacent_kmers: int, k_value: int, offset: int = 1,
                         position: bool = False, guarantee_tip: bool = False) -> list:
    """
    Determine minimizers for a input string.

    Determines minimizers for a string based on lexicographical
    order. Skips windows that cannot have a minimizer based on
    the minimizer computed in the previous iteration.

    Parameters
    ----------
    input_string : str
        String representing the sequence.
    adjacent_kmers : int
        Window size value. Number of adjacent k-mers per group.
    k_value : int
        Value of k for the k-mer size.
    offset : int, optional
        Value to indicate offset of consecutive k-mers.
    position : bool, optional
        If the start position of the k-mers in the sequence
        should be stored.
    guarantee_tip : bool, optional
        If guarantee that the minimizer cover the start and end of the sequence.

    Returns
    -------
    minimizers : list
        A list with the set of minimizers determined
        for the input string.
    """

    # break string into k-mers
    kmers: List[Union[str, Tuple[str, int]]] = string_kmerizer(input_string, k_value, offset, position)

    i: int = 0
    previous: Optional[Union[str, Tuple[str, int]]] = None
    sell: bool = False
    minimizers: List[Union[str, Tuple[str, int]]] = []
    # determine total number of windows
    last_window: int = (len(kmers) - adjacent_kmers)
      
    while i <= last_window:
        # get kmers in current window
        window: List[Union[str, Tuple[str, int]]] = kmers[i:i+adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer: List[Union[str, Tuple[str, int]]] = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx: int = window.index(minimizer[0])
        # sliding window that does not included last minimizer
        if previous is None:
            # simply store smallest minimizer
            minimizers.extend(minimizer)
        # sliding window includes last minimizer because we
        # skipped some sliding windows
        else:
            # check if minimizer is the same as the one picked
            # in the last window
            # Do not store minimizer if it is the same
            if minimizer[0] != previous:
                # get kmers smaller than last minimizer
                skipped: List[Union[str, Tuple[str, int]]] = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal: Union[str, Tuple[str, int]] = previous
                for m in skipped:
                    if m < minimal:
                        minimizer.append(m)
                        minimal = m
                minimizers.extend(minimizer)

        # slide by 1 if minimizer has index 0 in window
        if minimizer_idx == 0:
            i += 1
            previous = None
        # skip sliding windows based on minimizer position
        else:
            i += minimizer_idx
            # if adding minimizer index surpasses last window value we
            # might miss one last minimizer because it will fail the condition
            # find a better way to control this condition!
            if i > last_window and sell is False:
                i = last_window
                sell = True
            previous = minimizer[0]

    #Guarantee that first and last kmer are part of the minimizers list
    if guarantee_tip:
        if kmers[0] not in minimizers:
            minimizers.append(kmers[0])
        if kmers[-1] not in minimizers:
            minimizers.append(kmers[-1])
        
    return minimizers


def kmer_coverage(position: list, window_size: int) -> int:
    """
    Determines the sum size of kmers for coverage calculation
    
    Parameters
    ----------
    position : list
        List with starting positions sorted from initial to last.
    window_size : int
        Size of each minimizers window.
        
    Returns
    -------
    size : int
        Total sum of the kmers size based on starting position.
    """
    
    size: int = 0
    len_positions: int = len(position)-1
    for i, pos in enumerate(position):
        if i == 0:
            start = pos
            end = pos+window_size-1
        elif pos <= end:
            end = pos+window_size-1
        else:
            size += end-start+1
            start = pos
            end = pos+window_size-1
            
        if i == len_positions:
            end = pos+window_size-1
            size += end-start+1
    return size
