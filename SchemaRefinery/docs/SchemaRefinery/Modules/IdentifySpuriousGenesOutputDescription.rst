Output files and folders description for IdentifySpuriousGenes module
=====================================================================

**For --run-mode schema:**

::

    **OutputFolderName**
        The folder where the output files are stored.

    1_schema_processing
        Folder containing schema processing results.

        master_nucleotide.fasta
            Master FASTA file.

        schema
            Folder containing schema files.

            loci_x.fasta
                FASTA file for locus x.

            new_loci_y.fasta
                FASTA file for new locus y.

            ...
                Other loci files.

            short
                Folder containing short loci files.

                loci_x.fasta
                    Short FASTA file for locus x.

                new_loci_y.fasta
                    Short FASTA file for new locus y.

                ...
                    Other short loci files.

        schema_translation_folder
            Folder containing schema translations.

            loci_x.fasta
                Translation for locus x.

            new_loci_y.fasta
                Translation for new locus y.

            ...
                Other translations.

    2_BLAST_processing
        Folder containing BLAST processing results.

        1_BLASTn_processing
            Folder containing BLASTn processing results.

            blast_db_nucl
                Folder containing BLASTn database.

                Blast_db_nucleotide.ndb
                    BLASTn nucleotide database file.

                Blast_db_nucleotide.nhr
                    BLASTn nucleotide header file.

                Blast_db_nucleotide.nin
                    BLASTn nucleotide index file.

                Blast_db_nucleotide.nog
                    BLASTn nucleotide organism group file.

                Blast_db_nucleotide.nsd
                    BLASTn nucleotide sequence data file.

                Blast_db_nucleotide.nsi
                    BLASTn nucleotide sequence index file.

                Blast_db_nucleotide.nsq
                    BLASTn nucleotide sequence query file.

                Blast_db_nucleotide.ntf
                    BLASTn nucleotide taxonomy file.

                Blast_db_nucleotide.nto
                    BLASTn nucleotide taxonomy organism file.

            BLASTn_results
                Folder containing BLASTn results.

                blast_results_x.tsv
                    BLASTn results for x.

                blast_results_y.tsv
                    BLASTn results for y.

                ...
                    Other BLASTn results.

        2_BLASTp_processing
            Folder containing BLASTp processing results.

            blastn_results_matches_translations
                Folder containing BLASTn results matches translations.

                cluster_rep_translation
                    Folder containing BLASTn results of representative matches translations.

                    cluster_matches_translation_x.tsv
                        Cluster representatives matches translation for x.

                    cluster_matches_translation_y.tsv
                        Cluster representatives matches translation for y.

                    ...
                        Other cluster representatives matches translations.

                cluster_matches_translation_x.tsv
                    Cluster matches translation for x.

                cluster_matches_translation_y.tsv
                    Cluster matches translation for y.

                ...
                    Other cluster matches translations.

            BLASTp_results
                Folder containing BLASTp results.

                blast_results_x.tsv
                    BLASTp results for x.

                blast_results_y.tsv
                    BLASTp results for y.

                ...
                    Other BLASTp results.

            BLASTp_results_self_score_results
                Folder containing BLASTp self-score results.

                blast_results_x.tsv
                    BLASTp self-score results for x.

                blast_results_y.tsv
                    BLASTp self-score results for y.

                ...
                    Other BLASTp self-score results.

    3_processing_results
        Folder containing processing results.

        blast_results
            Folder containing BLAST results.

            blast_all_matches.tsv
                TSV file containing all BLAST matches.

            blast_by_cluster
                Folder containing BLAST results by cluster.

                blast_joined_cluster_x.tsv
                    BLAST results for cluster x.

                blast_retained_y.tsv
                    BLAST results for cluster y.

                ...
                    Other cluster results.

            blast_results_by_class
                Folder containing BLAST results by class.

                class_1a.tsv
                    BLAST results for class 1.

                class_2a.tsv
                    BLAST results for class 2.

                ...
                    Other class results.

        **count_results_by_cluster.tsv**
            TSV file containing count results by cluster.

        **drop_loci_reason.tsv**
            TSV file containing reasons for dropping loci.

        **recommendations.tsv**
            TSV file containing recommendations.

        **recommendations_annotations.tsv**
            TSV file containing the annotations of the recommendations loci.

        **related_matches.tsv**
            TSV file containing related matches.

**For --run-mode unclassified_cds:**

::

    **OutputFolderName**
        The folder where the output files are stored.

    1_CDS_processing
        Folder containing CDS processing results.

        CDS_not_found.fasta
            FASTA file containing CDS not found.

        CDS_not_found_translation.fasta
            FASTA file containing translations of CDS not found.

    2_BLAST_processing
        Folder containing BLAST processing results.

        1_BLASTn_processing
            Folder containing BLASTn processing results.

            blast_db_nucl
                Folder containing BLASTn database.

                Blast_db_nucleotide.ndb
                    BLASTn nucleotide database file.

                Blast_db_nucleotide.nhr
                    BLASTn nucleotide header file.

                Blast_db_nucleotide.nin
                    BLASTn nucleotide index file.

                Blast_db_nucleotide.nog
                    BLASTn nucleotide organism group file.

                Blast_db_nucleotide.nsd
                    BLASTn nucleotide sequence data file.

                Blast_db_nucleotide.nsi
                    BLASTn nucleotide sequence index file.

                Blast_db_nucleotide.nsq
                    BLASTn nucleotide sequence query file.

                Blast_db_nucleotide.ntf
                    BLASTn nucleotide taxonomy file.

                Blast_db_nucleotide.nto
                    BLASTn nucleotide taxonomy organism file.

            BLASTn_results
                Folder containing BLASTn results.

                blast_results_x.tsv
                    BLASTn results for x.

                blast_results_y.tsv
                    BLASTn results for y.

                ...
                    Other BLASTn results.

        2_BLASTp_processing
            Folder containing BLASTp processing results.

            blastn_results_matches_translations
                Folder containing BLASTn results matches translations.

                cluster_rep_translation
                    Folder containing BLASTn results of representative matches translations.

                    cluster_matches_translation_x.tsv
                        Cluster representatives matches translation for x.

                    cluster_matches_translation_y.tsv
                        Cluster representatives matches translation for y.

                    ...
                        Other cluster representatives matches translations.

                cluster_matches_translation_x.tsv
                    Cluster matches translation for x.

                cluster_matches_translation_y.tsv
                    Cluster matches translation for y.

                ...
                    Other cluster matches translations.

            BLASTp_results
                Folder containing BLASTp results.

                blast_results_x.tsv
                    BLASTp results for x.

                blast_results_y.tsv
                    BLASTp results for y.

                ...
                    Other BLASTp results.

            BLASTp_results_self_score_results
                Folder containing BLASTp self-score results.

                blast_results_x.tsv
                    BLASTp self-score results for x.

                blast_results_y.tsv
                    BLASTp self-score results for y.

                ...
                    Other BLASTp self-score results.

    3_processing_results
        Folder containing processing results.

        blast_results
            Folder containing BLAST results.

            blast_all_matches.tsv
                TSV file containing all BLAST matches.

            blast_by_cluster
                Folder containing BLAST results by cluster.

                blast_joined_cluster_x.tsv
                    BLAST results for cluster x.

                blast_retained_y.tsv
                    BLAST results for cluster y.

                ...
                    Other cluster results.

            blast_results_by_class
                Folder containing BLAST results by class.

                class_1a.tsv
                    BLAST results for class 1.

                class_2a.tsv
                    BLAST results for class 2.

                ...
                    Other class results.

        cds_id_changes.tsv
            TSV file containing changes in CDS IDs.

        dropped_cds.tsv
            TSV file containing dropped CDS.

        Graph_folder
            Folder containing graphs.

            All_of_CDS_graphs.html
                HTML file containing all CDS graphs.

            graphs_class_1a.html
                HTML file containing class 1a graphs.

            ...
                Other graph files.

        **count_results_by_cluster.tsv**
            TSV file containing count results by cluster.

        **drop_loci_reason.tsv**
            TSV file containing reasons for dropping loci.

        **recommendations.tsv**
            TSV file containing recommendations.

        **recommendations_annoations.tsv**
            TSV file containing the annotations of the recommendations loci.

        **related_matches.tsv**
            TSV file containing related matches.

        **temp_fastas**
            Folder containing temporary FASTA files.

            **cluster_x.fasta**
                Temporary FASTA file for cluster x.

            **cluster_y.fasta**
                Temporary FASTA file for cluster y.

            **...**
                Other temporary FASTA files.