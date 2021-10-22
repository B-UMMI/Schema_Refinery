import argparse
import os
import csv


def main(input_table_1, input_table_2, input_table_3, output_directory, bsr_threshold):

    processed_table_1 = []
    dict = {}
    final_table = []

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    

    # open input_table_1 "matches.tsv" outputed from "match_schemas-.py" script
    with open(input_table_1, 'r') as table1:
        lines_table_1 = list(csv.reader(table1, delimiter='\t'))

    # open input_table_2 with loci annotations
    with open(input_table_2, 'r') as table2:
        lines_table_2 = list(csv.reader(table2, delimiter='\t'))
    
     # open input_table_2 with loci annotations
    with open(input_table_3, 'r') as table3:
        lines_table_3 = list(csv.reader(table3, delimiter='\t'))

    #processing first table
    for row in lines_table_1:
        if ( float(row[2]) >= bsr_threshold ):
            new_list = []
            new_list.append(row[0])

            strs = row[1].split('_')
            new_list.append(strs[0])
            new_list.append(row[2])

            processed_table_1.append(new_list)
 

    #processing second table
    for row in lines_table_2:
        row[0] = row[0].split('.fasta')[0]
        dict[row[0]] = row


    for row in processed_table_1:

        entry = dict.get(row[1])

        #this check is done so the script completes without errors, there are still unmatched genes in the tables
        if entry:
            entry[0] = row[0]
            final_table.append(entry)
        else:
            final_table.append([row[0]])


    for row in lines_table_3:
        final_table.append([row[0]])


    
    with open(output_directory + '/new_annotations.tsv', 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        # writing the data rows 
        csvwriter.writerow(lines_table_2[0]) #put the collumn names
        csvwriter.writerows(final_table)



    with open(output_directory + '/processed_matches.tsv', 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        csvwriter.writerows(processed_table_1)

      



def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i1', '--input_table_1', type=str,
                        required=True, dest='input_table_1',
                        help='')
    
    parser.add_argument('-i2', '--input_table_2', type=str,
                        required=True, dest='input_table_2',
                        help='')

    parser.add_argument('-i3', '--input_table_3', type=str,
                        required=True, dest='input_table_3',
                        help='')

    parser.add_argument('-bsrt', '--bsr_threshold', type=float,
                        required=False, dest='bsr_threshold', default=0.7,
                        help='')

    parser.add_argument('-o', '--output_directory', type=str,
                    required=True, dest='output_directory',
                    help='Path to the directory where downloaded '
                            'files will be stored.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))