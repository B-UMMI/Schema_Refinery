import os
import argparse
from numpy import False_
import pandas as pd
import re
from Bio import SeqIO



def main(inputs, outputs, min_N):

    fastas = []

    for path in os.listdir(inputs):

        if os.path.isfile(os.path.join(inputs,path)):

            fastas.append(os.path.join(inputs,path))

    for fasta in fastas:
        records = SeqIO.parse(fasta, 'fasta')
        seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

        contigs = [re.split("N{"+ str(min_N) + ",100000}", rec) for rec in seqs_dict.values()]

        contigs = [j for i in contigs for j in i]

        with open(os.path.join(outputs,fasta.split("/")[-1]),'w+') as file:
            for number, contig in enumerate(contigs):
                file.write(f">contig_{number}\n")
                file.write(contig + "\n")
            




def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                            dest='inputs',
                            help='input assemblies')
    
    parser.add_argument('-o', type=str, required=True,
                            dest='outputs',
                            help='output')

    parser.add_argument('-min_N', type=str, required=False,
                            dest='min_N',
                            help='minimum consecutive N to exclude',
                            default = 6)
        

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    
    main(**vars(args))
