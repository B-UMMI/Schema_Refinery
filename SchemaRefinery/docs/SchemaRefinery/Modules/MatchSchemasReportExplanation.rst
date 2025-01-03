Output files and folders description for the MatchSchemas module
================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

    Subject_Translation
        Folder containing the translation files.
        
        subject_translations_x.fasta
            FASTA file containing the translation for the locus x.
        
        subject_translations_y.fasta
            FASTA file containing the translation for the locus y.
        
        subject_translations_z.fasta
            FASTA file containing the translation for the locus z.
        
        ...
            All of the other translation files.

    self_score_folder
        Folder containing the self-score BLAST results.
        
        blast_results_x.tsv
            TSV file containing the BLASTp results for self-score for the locus x.
        
        blast_results_y.tsv
            TSV file containing the BLASTp results for self-score for the locus y.
        
        blast_results_z.tsv
            TSV file containing the BLASTp results for the locus z.
        
        ...
            All of the other TSV BLASTp for self-score results files.

    **best_blast_matches.tsv**
        TSV file containing the best BLAST matches for the query and subject schemas.
