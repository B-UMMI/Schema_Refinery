Output files and folders description for the MatchSchemas module
================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

    blast_processing
        Folder containing BLASTp database, BLASTp output files, master file, self-score folder, and translation files.

        subject_master_file.fasta
            FASTA file containing all of the protein sequences that have not found a match yet used in the analysis (used to create BLAST DB).

        subject_master_rep_file.fasta
            FASTA file containing all of the representative protein sequences that have not found a match yet used in the analysis (used to create BLAST DB).

        Blast_db_prot
            Folder containing the BLASTp database.
            
            Blast_db_protein.pdb
                Position-specific Data Base file. Contains position-specific scoring matrices (PSSMs) used in PSI-BLAST searches.
            
            Blast_db_protein.phr
                Protein Header Record file. Contains the header information for each sequence in the protein database.
            
            Blast_db_protein.pin
                Protein Index file. Contains the index of the sequences in the protein database.
            
            Blast_db_protein.pog
                Protein Organism Group file. Contains information about the taxonomic grouping of the sequences in the protein database.
            
            Blast_db_protein.pos
                Protein Organism Sequence file. Contains the actual sequence data for the protein database.
            
            Blast_db_protein.pot
                Protein Organism Taxonomy file. Contains taxonomic information for the sequences in the protein database.
            
            Blast_db_protein.psq
                Protein Sequence Query file. Contains the sequence data in a format optimized for BLAST searches.
            
            Blast_db_protein.ptf
                Protein Taxonomy File. Contains taxonomy information for the sequences in the protein database.
            
            Blast_db_protein.pto
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

        Query_Not_Translated
            Folder containing the query translation files that were not translated due to any error. If this folder is empty (all translations were sucessful) the folder is deleted. 
            
            query_translations_x.fasta
                FASTA file containing the translation for the locus x.
            
            query_translations_y.fasta
                FASTA file containing the translation for the locus y.
            
            query_translations_z.fasta
                FASTA file containing the translation for the locus z.
            
            ...
                All of the other translation files.
        
        Query_Not_Translated_Rep
            Folder containing the representative query translation files. that were not translated due to any error. If this folder is empty (all translations were sucessful) the folder is deleted. 
            
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

    **existing_matches.txt
        TXT file containing the blast matches from all processes.

    **Match_Schemas_Results.tsv**
        TSV file containing the best BLAST matches for the query and subject schemas sorted by locus name. Also contains the final non matched locus from both Query and Subject.
