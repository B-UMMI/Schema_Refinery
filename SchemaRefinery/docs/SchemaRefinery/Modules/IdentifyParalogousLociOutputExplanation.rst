Output files and folders description IdentifyParalogousLoci module
===================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

    Blast
        Folder containing BLASTp database, BLASTp output files, master file, self-score folder, and translation files.

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

        Blast_output
            Folder containing the BLASTp output files.
            
            blast_results_x.tsv
                TSV file containing the BLASTp results for the locus x.
            
            blast_results_y.tsv
                TSV file containing the BLASTp results for the locus y.
            
            blast_results_z.tsv
                TSV file containing the BLASTp results for the locus z.
            
            ...
                All the other TSV BLASTp results files.

        master_file.fasta
            FASTA file containing all the protein sequences used in the analysis (used to create BLAST DB).

        self_score_folder
            Folder containing the self-score BLAST results.
            
            blast_results_x.tsv
                TSV file containing the BLASTp results for self-score for the locus x.
            
            blast_results_y.tsv
                TSV file containing the BLASTp results for self-score for the locus y.
            
            blast_results_z.tsv
                TSV file containing the BLASTp results for self-score for the locus z.
            
            ...
                All the other TSV BLASTp for self-score results files.

        Translation
            Folder containing the translation files.
            
            x_translation.fasta
                FASTA file containing the translation for the locus x.
            
            y_translation.fasta
                FASTA file containing the translation for the locus y.
            
            z_translation.fasta
                FASTA file containing the translation for the locus z.
            
            ...
                All the other translation files.

    **paralogous_annotations.tsv**
        TSV file with the clusters that should be joined as well as the annotations for each locus.
        Only with the -ann option.
    
    **paralogous_loci_final_recommendations.tsv**
        TSV with the loci clusters (one locus per row) with the action to be taken (always 'Join'). The clusters are separated by a row with a '#'.
        The clusters here differ from the ones in paralogous_loci_report_cluster_by_id.tsv because it consists of only the loci that passed certain checks.

    **paralogous_loci_report.tsv**
        TSV file containing the report of the paralogous loci. Here we can see which loci pass the checks to be considered paralogs that can be joined.

    **paralogous_loci_report_cluster_by_id.tsv**
        TSV file containing the report of the paralogous loci clustered by ID.
        Only with the --nocleanup option.