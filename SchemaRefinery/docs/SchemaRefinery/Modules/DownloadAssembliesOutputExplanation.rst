Output files and folders description for the DownloadAssemblies module
=======================================================================

**OutputFolderName**
    The folder where the output files are stored.

**assemblies_ncbi.zip**
    Zip file containing all the assemblies and extra information that the user wants downloaded from NCBI.

**ena661k_assemblies**
    Folder containing the assemblies downloaded from ENA661K.
    
    **x.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly x.
    
    **y.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly y.
    
    **z.contigs.fa.gz**
        Gzipped FASTA file containing the contigs for the assembly z.
    
    **...**
        Other gzipped FASTA files for the assemblies.

**metadata_all**
    Folder containing all the metadata downloaded from NCBI and ENA661K.
    
    **biosamples_ids.tsv**
        TSV file containing the BioSample IDs for the assemblies.
    
    **id_matches.tsv**
        TSV file containing the matches between the BioSample IDs and the assembly IDs and SRA IDs.
    
    **all_ids_fetched.tsv**
        TSV file containing all the IDs fetched during the download process.
