#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import sys
import csv
import time
import argparse
import urllib.request
from zipfile import ZipFile
import hashlib


maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


# ftp path
ebi_ftp = 'http://ftp.ebi.ac.uk'

# checklist.chk
url_hash_file = 'http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/checklist.chk'

# Dictionary of hash values
hashes_dict = {}

def HandleProgress(block_num, block_size, total_size):
        read_data= 0
        # calculating the progress
        # storing a temporary value  to store downloaded bytesso that we can add it later to the overall downloaded data
        temp = block_num * block_size
        read_data = temp + read_data
        #calculating the remaining size
        remaining_size = total_size - read_data
        if(remaining_size<=0):
            downloaded_percentage = 100
            remaining_size = 0
        else:
            downloaded_percentage = int(((total_size-remaining_size) / total_size)*(100))

        if downloaded_percentage == 100:
            print( f'Downloaded: {downloaded_percentage}%  ' , end="\n")
        else:
            print(" ", end="\r")   
            print( f'Downloaded: {downloaded_percentage}%  ' , end="\r")


def checkDownload(download:str, original_hash:str):
    """
        Function to check the integrity of a downloaded file.
    """

    with open(download, 'rb') as infile:
        data = infile.read()
        md5 = hashlib.md5(data).hexdigest()


    if md5 != original_hash:
        print('Hashes do not match')
        return False
    
    return True


def read_table(file_path, delimiter='\t'):
    """ Reads a tabular file.

        Parameters
        ----------
        file_path : str
            Path to the tabular file.
        delimiter : str
            Field delimiter.

        Returns
        -------
        lines : list
            List that contains one sublist per
            line in the input tabular file.
    """

    with open(file_path, 'r') as infile:
        lines = list(csv.reader(infile, delimiter=delimiter))

    return lines


def download_ftp_file(file_url, out_file, original_hash=''):
    """ Downloads a file from a FTP server.

        Parameter
        ---------
        file_url : str
            FTP path to the file to download.
        out_file : str
            Local path to the file that will be downloaded.

        Returns
        -------
        res : tup
            A tuple with the path to the local file
            and the response object.
    """

    tries = 0
    downloaded = False
    while downloaded is False and tries <= 5:
        try:
            res = urllib.request.urlretrieve(file_url, out_file, HandleProgress)
        except:
            time.sleep(1)
        tries += 1
        if os.path.isfile(out_file) is True:
            if original_hash != '':
                if checkDownload(out_file, original_hash):
                    downloaded = True
            else:
                downloaded = True

    return downloaded


def download_assemblies(sample_ids, sample_paths, output_directory):
    """ Downloads a set of assemblies from the FTP server of
        the "ENA2018-bacteria-661k" study.

        Parameters
        ----------
        sample_ids : list
            List with the identifiers of the samples/assemblies
            to download.
        sample_paths : dict
            Dictionary with sample/assemblies identifiers as
            keys and FTP paths as values.
        output_directory : str
            Path to the output directory.

        Returns
        -------
        failed : int
            Number of failed downloads.
    """

    # list files in output directory
    local_files = os.listdir(output_directory)

    # create URLs to download
    remote_urls = []
    for sid in sample_ids:
        sample_basename = sample_paths[sid].split('/')[-1]
        # do not download files that have already been downloaded
        if sample_basename not in local_files and sample_basename.split('.gz')[0] not in local_files:
            sample_file = os.path.join(output_directory, sample_basename)
            sample_url = ebi_ftp + sample_paths[sid]
            remote_urls.append([sample_url, sample_file, hashes_dict[sid]])

    if len(remote_urls) < len(sample_ids):
        print('{0} assemblies had already been downloaded.'
              ''.format(len(sample_ids)-len(remote_urls)))

    print('\nDownloading {0} assemblies...'.format(len(remote_urls)))
    failed = 0
    downloaded = 0
    for url in remote_urls:
        res = download_ftp_file(*url)
        if res is True:
            downloaded += 1
            if (downloaded == len(remote_urls)):
                print('Downloaded {0}/{1}'.format(downloaded,
                                                    len(remote_urls)), end="\n")
            else:
                print('Downloaded {0}/{1}'.format(downloaded,
                                                    len(remote_urls)), end="\033[F")
        else:
            failed += 1

    return failed

def download_assemblies_stride(sample_ids, sample_paths, output_directory, 
                                stride):
    """ Downloads a set of assemblies from the FTP server of
        the "ENA2018-bacteria-661k" study, in a specified interval.

        Parameters
        ----------
        sample_ids : list
            List with the identifiers of the samples/assemblies
            to download.
        sample_paths : dict
            Dictionary with sample/assemblies identifiers as
            keys and FTP paths as values.
        output_directory : str
            Path to the output directory.
        stride : str
            Interval of indexes to download.

        Returns
        -------
        failed : int
            Number of failed downloads.
    """


    interval_list = stride.split(':')
    low = int(interval_list[0]) -1
    high = int(interval_list[1])

    # make sure the interval doesnt go outside of the sample_ids list bounds
    if high > len(sample_ids):
        high = len(sample_ids)
        #update the stride string for the zip archive name
        stride = str(low + 1) + ':' + str(high)
        

    # create URLs to download
    remote_urls = []
    for i in range(low, high):
        sample_basename = sample_paths[sample_ids[i]].split('/')[-1]
        
        sample_file = os.path.join(output_directory, sample_basename)
        sample_url = ebi_ftp + sample_paths[sample_ids[i]]
        remote_urls.append([sample_url, sample_file, hashes_dict[sample_ids[i]]])


    print('\nDownloading {0} assemblies...'.format(len(remote_urls)))
    failed = 0
    downloaded = 0
    for url in remote_urls:
        res = download_ftp_file(*url)
        if res is True:
            downloaded += 1
            if (downloaded == len(remote_urls)):
                print('Downloaded {0}/{1}'.format(downloaded,
                                                    len(remote_urls)), end="\n")
            else:
                print('Downloaded {0}/{1}'.format(downloaded,
                                                    len(remote_urls)), end="\033[F")
        else:
            failed += 1

    # create a ZipFile object
    with ZipFile(output_directory + '/' + stride + '.zip', 'w') as zipObj:
        for file in os.listdir(output_directory):
            if '.gz' in file:
                zipObj.write(output_directory + '/' + file, file)
                #remove the .gz file
                os.remove(output_directory + '/' + file)


    return failed


def main(metadata_table, paths_table, species_name, output_directory,
         ftp_download, abundance, genome_size, size_threshold,
         max_contig_number, mlst_species, known_st, any_quality, stride):

    # read file with metadata
    metadata_lines = read_table(metadata_table)

    metadata_header = metadata_lines[0]

    # select lines based on species name
    species_index = metadata_header.index('species')
    species_lines = [line
                     for line in metadata_lines[1:]
                     if line[species_index] == species_name]

    print('\nFound {0} samples for species={1}.'
          ''.format(len(species_lines), species_name))

    if len(species_lines) == 0:
        print('Did not find matches for {0}.'.format(species_name))
        print('Please provide a valid species name.')
        sys.exit(0)

    # filter based on genome size
    if genome_size is not None:
        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)
        size_index = metadata_header.index('total_length')
        species_lines = [line
                         for line in species_lines
                         if int(line[size_index]) >= bot_limit
                         and int(line[size_index]) <= top_limit]

        print('{0} with genome size >= {1} and '
              '<= {2}.'.format(len(species_lines), bot_limit, top_limit))

    # filter based on species abundance
    if abundance is not None:
        abundance_index = metadata_header.index('adjust_abundance')
        species_lines = [line
                         for line in species_lines
                         if float(line[abundance_index]) >= abundance]

        print('{0} with abundance >= {1}.'.format(len(species_lines),
                                                  abundance))

    # filter based on number of contigs
    if max_contig_number is not None:
        contigs_index = metadata_header.index('total_contigs')
        species_lines = [line
                         for line in species_lines
                         if int(line[contigs_index]) <= max_contig_number]

        print('{0} with <= {1} contigs.'.format(len(species_lines),
                                                max_contig_number))

    # filter based on MLST species
    if mlst_species is not None:
        st_species_index = metadata_header.index('mlst-species')
        species_lines = [line
                         for line in species_lines
                         if line[st_species_index] == mlst_species]

        print('{0} with MLST species == {1}.'.format(len(species_lines),
                                                     mlst_species))

    # filter based on known ST
    if known_st is True:
        st_index = metadata_header.index('mlst')
        species_lines = [line
                         for line in species_lines
                         if line[st_index] != '-']

        print('{0} with known ST.'.format(len(species_lines)))

    # filter based on quality level
    if any_quality is False:
        quality_index = metadata_header.index('high_quality')
        species_lines = [line
                         for line in species_lines
                         if line[quality_index] == 'TRUE']

        print('{0} with high quality.'.format(len(species_lines)))

    # get sample identifiers
    sample_ids = [line[0] for line in species_lines]
    print('Selected {0} samples/assemblies that meet filtering '
          'criteria.'.format(len(sample_ids)))

    if len(sample_ids) == 0:
        sys.exit('Did not find samples/assemblies that passed '
                 'filtering criteria.')

    # create output directory
    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    selected_file = os.path.join(output_directory, 'selected_samples.tsv')
    with open(selected_file, 'w') as outfile:
        selected_lines = ['\t'.join(line)
                          for line in [metadata_header]+species_lines]
        selected_text = '\n'.join(selected_lines)
        outfile.write(selected_text)

    # download hashes file
    print('Downloading checklist.chk...')
    download_ftp_file(url_hash_file, 'checklist.chk')

    
    # Putting checksums in dictionary
    with open('checklist.chk', 'r') as table:

        lines = table.readlines()

        # hashes_dict = {rows[0]:rows[1] for rows in csv.reader(table, delimiter='\t')}
        for line in lines:
            if 'contigs.fa.gz' in line:
                vars = line.split('  ')
                hashes_dict[vars[1].split('/')[3].split('.')[0]] = vars[0]
    

    if ftp_download is True:

        # read table with FTP paths
        ftp_lines = read_table(paths_table)
        sample_paths = {l[0]: l[1].split('/ebi/ftp')[1] for l in ftp_lines}

        # download assemblies
        if stride:
            download_assemblies_stride(sample_ids, sample_paths, output_directory, stride)
        else:
            download_assemblies(sample_ids, sample_paths, output_directory)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', '--metadata-table', type=str,
                        required=True, dest='metadata_table',
                        help='Summarised QC and general characterisation for '
                             'each assembly in the "File4_QC_characterisati'
                             'on_661K" file from the study "Exploring '
                             'bacterial diversity via a curated and '
                             'searchable snapshot of archived DNA '
                             'sequences" (https://doi.org/10.1101/2021.03'
                             '.02.433662).')

    parser.add_argument('-p', '--paths-table', type=str,
                        required=True, dest='paths_table',
                        help='File with sample identifier to FTP path '
                             'mapping (available at http://ftp.ebi.ac.uk/'
                             'pub/databases/ENA2018-bacteria-661k/).')

    parser.add_argument('-s', '--species-name', type=str,
                        required=True, dest='species_name',
                        help='Name of the species. Must match one of the '
                             'species names in the "species" column in the '
                             'metadata table.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('--ftp-download', action='store_true',
                        required=False, dest='ftp_download',
                        help='If the assemblies from the selected samples'
                             'should be downloaded.')

    parser.add_argument('-a', '--abundance', type=float,
                        required=False,
                        dest='abundance',
                        help='Minimum species abundance. Samples with species'
                        ' abundance below this value are not selected.')

    parser.add_argument('-gs', '--genome-size', type=int,
                        required=False,
                        dest='genome_size',
                        help='Expected genome size.')

    parser.add_argument('-st', '--size-threshold', type=float,
                        required=False,
                        dest='size_threshold',
                        help='Genome size can vary in size +/- this value.')

    parser.add_argument('-mc', '--max_contig_number', type=int,
                        required=False,
                        dest='max_contig_number',
                        help='Maximum number of contigs. Assemblies with '
                             'a number of contigs greater than this value '
                             'are not selected.')

    parser.add_argument('--mlst-species', type=str,
                        required=False,
                        dest='mlst_species',
                        help='The species predicted by the MLST tool.')

    parser.add_argument('--known-st', action='store_true',
                        required=False,
                        dest='known_st',
                        help='If the samples must have a known ST.'
                             'Invalid or unkown STs will be "-".')

    parser.add_argument('--any-quality', action='store_true',
                        required=False,
                        dest='any_quality',
                        help='Download all assemblies, even the ones '
                             'that are not high quality.')

    parser.add_argument('-stride', '--stride', type=str,
                        required=False,
                        dest='stride',
                        help='Interval specifying which sample ids to download.'
                        'Example: "1:2000" - This will download the first 2000 samples.'
                        'Note: If you want to download from the first id, you have to put '
                        '"1", not "0" in the lower value.')


    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    print(args)
    main(**vars(args))
