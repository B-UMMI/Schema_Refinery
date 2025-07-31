Output files and folders description for the SchemaAnnotation module
====================================================================

::

    **OutputFolderName**
        The folder where the output files are stored.

        **genebank_annotations.tsv**
            Final tsv file with all the annotations from the Genbank and IDs to add.
            What is left after the cleanup.

        genbank_annotations
            Folder containing Genbank annotations.

            genebank_annotations.tsv 
                file with all the annotations from the Genbank and IDs to add.

            best_annotations_all_genbank_files
                Folder containing the best Genbank annotations.

                best_genbank_annotations.tsv
                    Best Genbank annotations.

            best_annotations_per_genbank_file
                Folder containing the best Genbank annotations per file.

                genbank_file_x_annotations.tsv
                    Genbank file x annotations.

                genbank_file_y_annotations.tsv
                    Genbank file y annotations.

                ...
                    Other Genbank file annotations.

            blast_processing
                Folder containing BLASTp database, BLASTp output files, and translation files.

                selected_genbank_proteins.fasta
                    Selected Genbank proteins.

                blast_db
                    Folder containing the BLASTp database.

                    blast_db_protein.pdb
                        Position-specific Data Base file.

                    blast_db_protein.phr
                        Protein Header Record file.

                    blast_db_protein.pin
                        Protein Index file.

                    blast_db_protein.pog
                        Protein Organism Group file.

                    blast_db_protein.pos
                        Protein Organism Sequence file.

                    blast_db_protein.pot
                        Protein Organism Taxonomy file.

                    blast_db_protein.psq
                        Protein Sequence Query file.

                    blast_db_protein.ptf
                        Protein Taxonomy File.

                    blast_db_protein.pto
                        Protein Taxonomy Organism file.

                blastp_results
                    Folder containing BLASTp results.

                    blast_results_x.tsv
                        BLAST results for x.

                    blast_results_y.tsv
                        BLAST results for y.

                    ...
                        Other BLAST results.

                self_score_folder
                    Folder containing self-score results.

                    blast_results_x.tsv
                        BLAST results for x.

                    blast_results_y.tsv
                        BLAST results for y.

                    ...
                        Other BLAST results.

            reps_translations
                Folder containing translations.

                x_translation.fasta
                    Translation for x.

                y_translation.fasta
                    Translation for y.

                ...
                    Other translations.


        **matched_annotations.tsv**
            Final file containing all annotations from the match schema option.
            What is left after cleanup.

        **uniprot_annotations.tsv**
            Final file containing all the UniProt annotations.
            What is left after cleanup.

        uniprot_annotations
            Folder containing UniProt annotations.

            best_proteomes_annotations_swiss_prot.tsv
                Best annotations for Swiss-Prot proteomes.

            best_proteomes_annotations_trEMBL.tsv
                Best annotations for TrEMBL proteomes.

            proteome_matcher_output
                Folder containing proteome matcher output.

                best_annotations_per_proteome_file
                    Folder containing the best annotations per proteome file.

                    Swiss-Prot
                        Folder containing Swiss-Prot annotations.

                        proteome_file_x_Swiss-Prot_annotations.tsv
                            Swiss-Prot annotations for proteome file x.

                        proteome_file_y_Swiss-Prot_annotations.tsv
                            Swiss-Prot annotations for proteome file y.

                        ...
                            Other Swiss-Prot annotations.

                    TrEMBL
                        Folder containing TrEMBL annotations.

                        proteome_file_x_TrEMBL_annotations.tsv
                            TrEMBL annotations for proteome file x.

                        proteome_file_y_TrEMBL_annotations.tsv
                            TrEMBL annotations for proteome file y.

                        ...
                            Other TrEMBL annotations.

                reps_translations
                    Folder containing translations.

                    x_translation.fasta
                        Translation for x.

                    y_translation.fasta
                        Translation for y.

                    ...
                        Other translations.

                self_score_folder
                    Folder containing self-score results.

                    blast_results_x.tsv
                        BLAST results for x.

                    blast_results_y.tsv
                        BLAST results for y.

                    ...
                        Other BLAST results.

                swiss_prots_processing
                    Folder containing Swiss-Prot processing results.

                    blast_processing
                        Folder with the files of the Blast database and the BLASTp results.
                        
                        swiss_prots.fasta
                            Swiss-Prot protein sequences.

                    swiss_prots_annotations.tsv
                        Swiss-Prot annotations.

                trembl_prots_processing
                    Folder containing TrEMBL processing results.

                    blast_processing
                        Folder with the files of the Blast database and the BLASTp results.

                        trembl_prots.fasta
                            TrEMBL protein sequences.

                    trembl_prots_annotations.tsv
                        TrEMBL annotations.
            
            Proteomes
                Folder containing the zipped fasta files of each proteome.
                
                Proteome_x.fasta.gz
                    Fasta file with the sequence of the proteome x.
            
            split_proteomes
                folder with the files with all the sequences for each database.

                prots_descriptions
                    File with the information of each sequence.
                
                swiss_prot.fasta
                    File with all the proteome sequences that are from the Swiss-Prot database.
                
                trembl_prot.fasta
                    File with all the proteome sequences that are from the TrEMBL database.
        
        **consolidated_annotations**
            Final file containing all the annotations of the inputs in a single file.