Output files and folders description for the MatchSchemas module
================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

    blast_processing
        Folder containing BLASTp database, BLASTp output files, master file, self-score folder, and translation files.

        subject_reps_vs_reps_blastdb
            Folder containing the BLASTp database for the rep vs rep analysis.
            
            subject_reps_vs_reps_blastdb.pdb
                Position-specific Data Base file. Contains position-specific scoring matrices (PSSMs) used in PSI-BLAST searches.
            
            subject_reps_vs_reps_blastdb.phr
                Protein Header Record file. Contains the header information for each sequence in the protein database.
            
            subject_reps_vs_reps_blastdb.pin
                Protein Index file. Contains the index of the sequences in the protein database.
            
            subject_reps_vs_reps_blastdb.pog
                Protein Organism Group file. Contains information about the taxonomic grouping of the sequences in the protein database.
            
            subject_reps_vs_reps_blastdb.pos
                Protein Organism Sequence file. Contains the actual sequence data for the protein database.
            
            subject_reps_vs_reps_blastdb.pot
                Protein Organism Taxonomy file. Contains taxonomic information for the sequences in the protein database.
            
            subject_reps_vs_reps_blastdb.psq
                Protein Sequence Query file. Contains the sequence data in a format optimized for BLAST searches.
            
            subject_reps_vs_reps_blastdb.ptf
                Protein Taxonomy File. Contains taxonomy information for the sequences in the protein database.
            
            subject_reps_vs_reps_blastdb.pto
                Protein Taxonomy Organism file. Contains organism-specific taxonomy information for the sequences in the protein database.

        subject_reps_vs_alleles_blastdb
            Folder containing the BLASTp database for the rep vs rep analysis.
            
            subject_reps_vs_alleles_blastdb.pdb
                Position-specific Data Base file. Contains position-specific scoring matrices (PSSMs) used in PSI-BLAST searches.
            
            subject_reps_vs_alleles_blastdb.phr
                Protein Header Record file. Contains the header information for each sequence in the protein database.
            
            subject_reps_vs_alleles_blastdb.pin
                Protein Index file. Contains the index of the sequences in the protein database.
            
            subject_reps_vs_alleles_blastdb.pog
                Protein Organism Group file. Contains information about the taxonomic grouping of the sequences in the protein database.
            
            subject_reps_vs_alleles_blastdb.pos
                Protein Organism Sequence file. Contains the actual sequence data for the protein database.
            
            subject_reps_vs_alleles_blastdb.pot
                Protein Organism Taxonomy file. Contains taxonomic information for the sequences in the protein database.
            
            subject_reps_vs_alleles_blastdb.psq
                Protein Sequence Query file. Contains the sequence data in a format optimized for BLAST searches.
            
            subject_reps_vs_alleles_blastdb.ptf
                Protein Taxonomy File. Contains taxonomy information for the sequences in the protein database.
            
            subject_reps_vs_alleles_blastdb.pto
                Protein Taxonomy Organism file. Contains organism-specific taxonomy information for the sequences in the protein database.

        blastp_results
            Folder containing the BLASTp output files.
            
            blast_results_x.tsv
                TSV file containing the BLASTp results for the locus x.
            
            blast_results_y.tsv
                TSV file containing the BLASTp results for the locus y.
            
            blast_results_z.tsv
                TSV file containing the BLASTp results for the locus z.
            
            ...
                All of the other TSV BLASTp results files.

        Query_Translation
            Folder containing the query translation files.
            
            query_translations_x.fasta
                FASTA file containing the translation for the locus x.
            
            query_translations_y.fasta
                FASTA file containing the translation for the locus y.
            
            query_translations_z.fasta
                FASTA file containing the translation for the locus z.
            
            ...
                All of the other translation files.
        
        Query_Translation_Rep
            Folder containing the representative query translation files.
            
            query_translations_x.fasta
                FASTA file containing the translation for the locus x.
            
            query_translations_y.fasta
                FASTA file containing the translation for the locus y.
            
            query_translations_z.fasta
                FASTA file containing the translation for the locus z.
            
            ...
                All of the other translation files.

        Subject_Translation
            Folder containing the subject translation files.
            
            subject_translations_x.fasta
                FASTA file containing the translation for the locus x.
            
            subject_translations_y.fasta
                FASTA file containing the translation for the locus y.
            
            subject_translations_z.fasta
                FASTA file containing the translation for the locus z.
            
            ...
                All of the other translation files.

        Subject_Translation_Rep
            Folder containing the representative subject translation files.
            
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

    **hashes_dna_matches.tsv**
        TSV file with the matches from the DNA hashes analysis.
    
    **hashes_prot_matches.tsv**
        TSV file with the matches from the proteins hashes analysis.

    **reps_vs_reps_matches.tsv**
        TSV file with the matches from the reps vs reps Blastp analysis.
    
    **reps_vs_alleles_matches.tsv**
        TSV file with the matches from the reps vs alleles Blastp analysis.

    **unmatched.tsv**
        TSV file with the loci that were not matched analysis.

    **Match_Schemas_Results.tsv**
        TSV file containing the best BLAST matches for the query and subject schemas sorted by locus name. Also contains the final non matched locus from both Query and Subject.
