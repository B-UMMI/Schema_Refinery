def string_kmerizer(input_string, k_value, offset=1, position=False):
    """Decompose a string into k-mers.

    Parameters
    ----------
    input_string : str
        String to divide into k-mers.
    k_value : int
        Value for the size of k-mers.
    offset : int
        Value to indicate offset of consecutive k-mers.
    position : bool
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
    if position is False:
        kmers = [input_string[i:i+k_value]
                 for i in range(0, len(input_string)-k_value+1, offset)]
    elif position is True:
        kmers = [(input_string[i:i+k_value], i)
                 for i in range(0, len(input_string)-k_value+1, offset)]

    return kmers

def determine_minimizers(input_string, adjacent_kmers, k_value, offset=1,
                         position=False):
    """Determine minimizers for a input string.

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
    offset : int
        Value to indicate offset of consecutive k-mers.
    position : bool
        If the start position of the k-mers in the sequence
        should be stored.

    Returns
    -------
    minimizers : list
        A list with the set of minimizers determined
        for the input string.
    """
    # break string into k-mers
    kmers = string_kmerizer(input_string, k_value, offset, position)

    i = 0
    previous = None
    sell = False
    minimizers = []
    # determine total number of windows
    last_window = (len(kmers)-adjacent_kmers)
    while i <= last_window:
        # get kmers in current window
        window = kmers[i:i+adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx = window.index(minimizer[0])
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
                skipped = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal = previous
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

    return minimizers