#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Inês Almeida
"""

import os
import urllib
import time
import hashlib
import argparse
import requests
import json
import csv
import subprocess

base_path = 'https://'


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


def HandleProgress(block_num, block_size, total_size):
    read_data = 0
    # calculating the progress
    # storing a temporary value  to store downloaded bytesso that we can add it later to the overall downloaded data
    temp = block_num * block_size
    read_data = temp + read_data
    # calculating the remaining size
    remaining_size = total_size - read_data
    if(remaining_size <= 0):
        downloaded_percentage = 100
        remaining_size = 0
    else:
        downloaded_percentage = int(
            ((total_size-remaining_size) / total_size)*(100))

    if downloaded_percentage == 100:
        print(f'Downloaded: {downloaded_percentage}%  ', end="\n")
    else:
        print(" ", end="\r")
        print(f'Downloaded: {downloaded_percentage}%  ', end="\r")


def checkDownload(download: str, original_hash: str):
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


def download_ftp(file_url, out_file, original_hash):
    downloaded = False
    while downloaded is False:
        try:
            res = urllib.request.urlretrieve(
                file_url, out_file, HandleProgress)
        except:
            time.sleep(1)
        if os.path.isfile(out_file) is True:
            if original_hash != '':
                if checkDownload(out_file, original_hash):
                    downloaded = True
            else:
                downloaded = True

    return downloaded


def download_from_ena(sample, outfile1, outfile2):
    json_fields = json.loads(requests.get(
        f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sample}&format=JSON&result=read_run&fields=fastq_ftp,fastq_md5').text)[0]
    download_links = json_fields['fastq_ftp'].split(';')
    md5_hashes = json_fields['fastq_md5'].split(';')

    print('Downloading 1')
    while not download_ftp(f'{base_path}{download_links[0]}', outfile1, md5_hashes[0]):
        continue

    print('Downloading 2')
    while not download_ftp(f'{base_path}{download_links[1]}', outfile2, md5_hashes[1]):
        continue

    print("Uncompressing files..")
    res = subprocess.run(["gzip", "-d", outfile1])
    if (res.returncode != 0):
        print("The exit code was: %d" % res.returncode)
        print("Something went wrong...")

    res = subprocess.run(["gzip", "-d", outfile2])
    if (res.returncode != 0):
        print("The exit code was: %d" % res.returncode)
        print("Something went wrong...")


def download_from_sra(sample, outfile1, outfile2):

    fastq_command = ["fastq-dump","--gzip", sample, "--split-files", "-v"]
    res = subprocess.run(fastq_command)
    if (res.returncode != 0):
        print("The exit code was: %d" % res.returncode)
        print("Something went wrong...")

    print(f"Moving file 1")
    mv_command = ["mv", f'{sample}_1.fastq.gz', outfile1[:-3]]
    print()
    res = subprocess.run(mv_command)
    if (res.returncode != 0):
        print("The exit code was: %d" % res.returncode)
        print("Something went wrong...")

    print(f"Moving file 2")
    mv_command = ["mv", f'{sample}_2.fastq.gz', outfile2[:-3]]
    res = subprocess.run(mv_command)
    if (res.returncode != 0):
        print("The exit code was: %d" % res.returncode)
        print("Something went wrong...")


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-a', '--accession_ids', type=str,
                        required=True, dest='accession_ids_file',
                        help='File with accession identifiers.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    args = parser.parse_args()

    return args


def main(accession_ids_file, output_directory):

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Reading through the acession id's file
    accession_ids = read_table(accession_ids_file)

    for id in accession_ids:
        sample = id[0]
        # Accession id directory
        accession_id_dir = f'{output_directory}/{sample}'

        if not os.path.exists(accession_id_dir):
            os.mkdir(accession_id_dir)

        outfile1 = f'{accession_id_dir}/{sample}_1.fastq.gz'
        outfile2 = f'{accession_id_dir}/{sample}_2.fastq.gz'

        # check_link = f'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=ena&format=JSON&includeAccessions={sample}&result=read_run&sortDirection=asc'

        # r = requests.get(check_link)
        # if (r.status_code == 200):
        #     print(f"Downloading {sample} from ENA...")
        #     download_from_ena(sample, outfile1, outfile2)
        # elif(r.status_code == 204):
        #     print(f"Downloading {sample} with SRA Toolkit...")
        #     download_from_sra(sample, outfile1, outfile2)
        # else:
        #     print(f'A different status code returned for sample: {sample}')

        download_from_sra(sample, outfile1, outfile2)


if __name__ == "__main__":

    args = parse_arguments()
    print(args)
    main(**vars(args))
