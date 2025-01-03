Output files and folders description IdentifyParalogousLoci module
===================================================================
::

    **OutputFolderName**
        The folder where the output files are stored.

    **Blast**
        Folder containing BLASTp database, BLASTp output files, master file, self-score folder, and translation files.

        **Blast_db_prot**
            Folder containing the BLASTp database.
            
            **Blast_db_protein.pdb**
                Position-specific Data Base file. Contains position-specific scoring matrices (PSSMs) used in PSI-BLAST searches.
            
            **Blast_db_protein.phr**
                Protein Header Record file. Contains the header information for each sequence in the protein database.
            
            **Blast_db_protein.pin**
                Protein Index file. Contains the index of the sequences in the protein database.
            
            **Blast_db_protein.pog**
                Protein Organism Group file. Contains information about the taxonomic grouping of the sequences in the protein database.
            
            **Blast_db_protein.pos**
                Protein Organism Sequence file. Contains the actual sequence data for the protein database.
