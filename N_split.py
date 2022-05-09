#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script serves to split assemblies by their consecutive missing data (N).

Authors : Mykyta Forofontov ,Raquel Romão
Code documentation
------------------
"""

"""
Inputs:

    -i folder with assemblies

    -o output folder where modified assemblies are stored

    -min_N threshold for number N, default = 2, includes 2 consecutive Ns.

    -t number of threads to use, default = 1

Outputs:
    Assemblies separated by min_N consecutive Ns, devided in to new contigs.
"""

import os
import argparse
import re
import concurrent.futures
from itertools import repeat
from Bio import SeqIO



def split(fasta, min_N):

    """
    Splits by min_N consecutive Ns in to new contigs

    Arguments:
        fasta: str
            assembly path

        min_N : int
            cuts by >= min_N consecutive Ns

    Returns :
        locus_mean : list
            containing ordered contigs, fasta path and identifier to add to new contigs
    """
    if sum([rec.seq.count('N') for rec in SeqIO.parse(fasta, 'fasta')]) != 0:

        records = SeqIO.parse(fasta, 'fasta')

        contigs = {rec.id : re.split("N{"+ str(min_N) + ",100000}", str(rec.seq.upper())) for rec in records}

        return [contigs,fasta]

    else:
        return [None,fasta]

def main(inputs, outputs, min_N,threads):
    """
    main body of the script
    """

    if not os.path.exists(outputs):
        os.mkdir(outputs)

    print("Reading folder")

    fastas = []

    for path in os.listdir(inputs):

        if os.path.isfile(os.path.join(inputs,path)):

            fastas.append(os.path.join(inputs,path))

    print("Removing N and splitting into contigs")

    """
    Multithreading:
        arguments:
            function : split
            input 1 : list
                all assemblies path
            input 2 : int
                repeats of min_N
    
        output: fasta file format
            creates fasta file format at output folder
    """

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        
        for res in executor.map(split, fastas,repeat(min_N)):


            newcontigs = 0

            if res[0] != None:
            
                with open(os.path.join(outputs,res[1].split("/")[-1]),'w+') as file:

                    for (id,contigs) in res[0].items():
                        index = 0
                        for fasta in contigs:

                            if len(fasta) > 0 :
                                    
                                file.write(f">{id}_{index}\n")
                                file.write(fasta + "\n")

                                if index != 0:
                                    newcontigs += 1
                                
                                index += 1

                print(res[1].split("/")[-1])
                print("new contigs = {}".format(newcontigs))
            
            else:

                print(res[1].split("/")[-1])
                print("No Ns detected.")

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                            dest='inputs',
                            help='input assemblies')
    
    parser.add_argument('-o', type=str, required=True,
                            dest='outputs',
                            help='output')

    parser.add_argument('-min_N', type=int, required=False,
                            dest='min_N',
                            help='minimum consecutive N to exclude',
                            default = 2)

    parser.add_argument('-t', type=int, required=False,
                            dest='threads',
                            help='number of cpu to use',
                            default = 1)
        

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    
    main(**vars(args))
