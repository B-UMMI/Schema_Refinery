import os

def main(schema_path):
    # Get all of the schema loci FASTA files path
    schema_loci = {loci_path.replace(".fasta", ""): loci_path 
                   for loci_path in os.listdir(schema_path) 
                   if loci_path.endswith('.fasta')}
    
    # Get all of the schema loci short FASTA files path
    schema_short_path = os.apth.join(schema_path, 'short')
    schema_loci_short = {loci_path.replace(".fasta", ""): loci_path 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}