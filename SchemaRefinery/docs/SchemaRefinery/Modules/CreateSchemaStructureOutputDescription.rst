Output files and folders description for the CreateSchemaStructure module
==========================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

    **schema**
        The folder where the final created schema will be stored.

        **x.fasta**
            The fasta file containing the loci.

        **y.fasta**
            The fasta file containing the loci.

        **z.fasta**
            The fasta file containing the loci.

        **...**
            Other fasta files containing the loci.

        **short**
            The folder containing the short loci.
            
            **x_short.fasta**
                The short fasta file containing the loci.
            
            **y_short.fasta**
                The short fasta file containing the loci.
            
            **z_short.fasta**
                The short fasta file containing the loci.
            
            **...**
                Other short fasta files containing the loci.

    **temp_fasta**
        Folder where the temporary fasta files (already altered but not structured) will be stored.

        **x.fasta**
            The fasta file containing the loci.

        **y.fasta**
            The fasta file containing the loci.

        **z.fasta**
            The fasta file containing the loci.

        **...**
            Other fasta files containing the loci.
    
    **schema_invalid_alleles.txt**
        List of all the identifiers of the alleles that were excluded and the reason for the exclusion of each allele.
        If none is found then the file is empty.

    ** schema_invalid_loci.txt**
        List of genes that had no valid alleles, one gene identifier per line.

    ** schema_summary_stats.tsv**
        Summary statistics for each gene (number of alleles in the external schema, number of valid alleles included in the adapted schema and number of representative alleles chosen by chewBBACA).
