#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 20 17:05:00 2022

@author: Catarina Mendes, Mykyta Forofontov

"""

"""
Purpose
-------
This script serves to analyse loci from different groups,
identifying alleles that are shared or exlcusive.
Code documentation
------------------
"""

"""
Inputs:
    -p --profile tsv file containing cgMLST profiles

    -g --group tsv file with groups, each column is a group and each
        line contains the id of assembly and 0(not in group) or 1(in group)

    -o , path to output file

    -a --absent how missing data is identified in profile file

Outputs:
    two txt files, one containing ids of loci that are  not shared between desired groups
    and another txt file containing shared loci with their respective shared allele ids, as
    in input profile matrix.
"""


import argparse
from asyncore import write 
import csv 
import sys
import os
from pickle import TRUE
import itertools

class Logger(object):    
    '''Class to create a log file of all that is printed onto the console''' 

    def __init__(self, out_directory):        
        self.logfile = os.path.join(out_directory, "run.log")        
        if os.path.isfile(self.logfile):            
            print(("Logfile already exists! It will be overwritten..." + "\n"))        
        self.terminal = sys.stdout        
        self.log = open(self.logfile, "w")    

    def write(self, message):        
        self.terminal.write(message)        
        self.log.write(message)        
        self.log.flush()    
    
    def flush(self):        
        pass
    
def loadGroups(filename):    
    '''Load group file into a dictionary'''

    group_dic = {}    
    
    with open(filename, 'r') as groups:        
        reader = list(zip(*csv.reader(groups, delimiter='\t')))        

        if len(reader) == 1:

            for i in range(1,len(reader[0])):

                new_group = list(0 for a in range(0,len(reader[0])-2))

                reader.append(new_group)

                reader[i].insert(0,reader[0][i])

                reader[i].insert(i,1)

        for group in range(1, len(reader)):            
            p = dict(list(zip(reader[0], reader[group])))

            if "" in p:                
                name_trait = p[""]                
                del p[""]            
            
            elif "Name" in p:                
                name_trait = p["Name"]                
                del p["Name"]            

            else:                
                sys.exit("Make sure the top-left cell in the traits file is either empty or 'Name'.")            

            for key, value in list(p.items()):                
                if int(value) == 0:                    
                    del p[key]            
            
            group_dic[name_trait] = list(p.keys())    
            
    print("Comparing groups:")

    for key, value in list(group_dic.items()):        
        print(key)        
        print((" ...containing %s isolates" % str(len(value))))

    print('')    
    
    groups = list(group_dic.keys())    
    
    return group_dic, groups
    
def loadProfile(filename, dataStructure):    
    '''Load profile file into the program, separating the samples into the groups'''    
    
    with open(filename, 'r') as profile:        
        reader = csv.reader(profile, delimiter='\t')        
        genes = reader.__next__()[1:]        
        profileDic = {}        
        
        for line in reader:            
            skipIsolate = True            
            strain = line[0]            
            profile = line[1:]  # has the same order as in genes            
            
            for key, value in list(dataStructure.items()):                
                
                if strain in value:                   
                    group = key                    
                    skipIsolate = False                    
                    break                
                
                else:                    
                    pass            
            
            if not skipIsolate:                
                
                for i in range(0, len(genes)):                    
                    
                    if genes[i] not in list(profileDic.keys()):                        
                        temp = {}                        
                        temp[group] = [profile[i]]                        
                        profileDic[genes[i]] = temp                    
                    
                    else:                        
                        temp = profileDic[genes[i]]                       
                        
                        if group in list(temp.keys()):                            
                            temp[group].append(profile[i])                        
                        
                        else:                            
                            temp[group] = [profile[i]]                        
                            profileDic[genes[i]] = temp            
            
            else:                
                pass    
    
    print(("profile size: " + str(len(profileDic))))        
    
    return profileDic

def proCompare(dataStructure, key1, key2,absent):    
    '''Comparison of the the two groups in the profile''' 

    shared_allele = []    
    exclusive_alleles = []
    allele_id_shared = []
    allele_id_exclusive = []
    # performing counts

    for gene, value in list(dataStructure.items()):

        if absent in dataStructure[gene][key1] or absent in dataStructure[gene][key2]:			
            pass

        else:
            # group one        
            group_one = set(value[key1])        
            
            # group two        
            group_two = set(value[key2])        
            
            # common alleles in the two groups        
            common_alleles = len(group_one.intersection(group_two))        
            
            if common_alleles != 0:            
                shared_allele.append(gene)        
                allele_id_shared.append(group_one.intersection(group_two))

            else:            
                exclusive_alleles.append(gene)
                allele_id_exclusive.append(group_one.difference(group_two))        
                
            # Alleles present only in group one        
            group_one_diff = len(group_one.difference(group_two))        
            
            # Alleles present only in group two        
            group_two_diff = len(group_two.difference(group_one))        
            
            value['Stats'] = [group_one, group_two, common_alleles, group_one_diff, group_two_diff] # shows alleles id for each group, how many common, and how many different
        
    print(('\n Comparing: ' + str(key1) + ' + ' + str(key2)))    
    print(("\tLoci with shared alleles: " + str(len(shared_allele))))    
    print(("\tLoci with exclusive alleles: " + str(len(exclusive_alleles))))    
    
    return dataStructure, shared_allele, exclusive_alleles,allele_id_shared,allele_id_exclusive
    
def main(profile,group,absent,outdir):       
    
    if not os.path.isdir(outdir):        
        os.makedirs(outdir)    
    
    sys.stdout = Logger("./")    
    
    # Loads group file   
    groups_dic, groups_list = loadGroups(group)    
      
    # Loads profile file    
    profiles = loadProfile(profile, groups_dic)    
    
    # performs profile comparison for two groups at a time     
    for a, b in itertools.combinations(groups_list, 2):

        _, shared_allele, exclusive_alleles,allele_id_shared,allele_id_exclusive = proCompare(profiles, a, b,absent)        
    
        with open(os.path.join(outdir,                               
                                '{type}_{a}_{b}.txt'.format(type='shared', a=a, b=b)),'w+') as writer:

            for index,i in enumerate(shared_allele):
                writer.writelines(i+'\n')
                writer.writelines("{}\n".format(str(allele_id_shared[index]).replace("'","")))
        
        with open(os.path.join(outdir,                               
                                '{type}_{a}_{b}.txt'.format(type='exclusive', a=a, b=b)),'w+') as writer:

            for index,i in enumerate(exclusive_alleles):
                writer.writelines(i+'\n')
                writer.writelines("{}\n".format(str(allele_id_exclusive[index]).replace("'",""))) 

        with open(os.path.join(outdir,                               
                                'stats_procompare.tsv'),'a') as writer:
    
            writer.writelines('Comparing: ' + str(a) + ' + ' + str(b) + '\n')  
            writer.writelines("Loci with shared alleles: " + str(len(shared_allele)) + '\n')   
            writer.writelines("Loci with exclusive alleles: " + str(len(exclusive_alleles)) + '\n')
                    
        
    print("\nFinished")    
    
    
    sys.exit(0)
    
def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p','--profile', type=str, required=True,
                        dest='profile',
                        help='Input profile tab file')

    parser.add_argument('-g','--group', type=str, required=True,
                        dest='group',
                        help='group tsv file')

    parser.add_argument('-o','--outdir', type=str, required=True,
                        dest='outdir',
                        help='outdir path')

    parser.add_argument('-a','--absent', type=str, required=True,
                        dest='absent',
                        help='character representing missing data')
    

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))