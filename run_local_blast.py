#!/usr/bin/python

# This script runs the local version of (PSI-)BLAST.

import argparse
import subprocess
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import itertools


def blast(db, query, query_folder="./queries/", v_blast="blast"):
    """
    This function executes blast or psi-blast for the given query and db.
    :param db: database filename
    :param query: query filename
    :param query_folder: query folder name
    :param v_blast: blast if BLAST, psiblast if PSI-BLAST
    :return: result from blast run
    """
    if v_blast == "blast":
        ##########################
        ### START CODING HERE ####
        ##########################
        # Define the variable 'cmd' as a string with the command for BLASTing 'query' against
        # the specified database 'db'.
        cmd = f"blastp -query {query_folder}{query}.fasta -db {db} -outfmt 6"
        print(cmd)

        ##########################
        ###  END CODING HERE  ####
        ##########################
    else:
        ##########################
        ### START CODING HERE ####
        ##########################
        # Define the variable 'cmd' as a string with the command for PSI-BLASTing 'query' against
        # the specified database 'db'.
        # Note that it is is easier to parse the output if it is in tabular format.
        # For that use can use the option -outfmt 6. Check this link https://www.ncbi.nlm.nih.gov/books/NBK279682/
        cmd = f"psiblast -query {query_folder}{query}.fasta -db {db} -outfmt 6 -num_iterations 100"
        print(cmd)

        ##########################
        ###  END CODING HERE  ####
        ##########################

    # Running shell command in python script. See https://docs.python.org/2/library/subprocess.html#popen-constructor
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         close_fds=True)
    blast_result = p.stdout.read().decode("utf8")
    return blast_result


def parse_blast_result(blast_result, blast_dict):
    """
    This function parses the output of (PSI-)BLAST and stores the result in blast_dict (defined in main()).
    :param blast_result: output  obtained after running (PSI-)BLAST
    :param blast_dict: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values]
    :return: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values]
    """
    for line in blast_result.split("\n"):
        if line and line[0] != "#" and line[0] != "[" and line[0] != 'W' and line[0] != 'S':
            try:
                splitted_line = line.split('\t')
                query = splitted_line[0].split("|")[1]
                subject = splitted_line[1]
                ##########################
                ### START CODING HERE ####
                ##########################
                # Create the variable evalue and assign the right value obtained from the splitted_line variable
                evalue = splitted_line[-2]

                ##########################
                ###  END CODING HERE  ####
                ##########################
                blast_dict[(query, subject)] = float(evalue)

            except IndexError:
                print("Could not parse (psi-)blast response line", line)
    return blast_dict
            

def write_output(uniprot_ids, output_filename, blast_dict):
    """
    This function writes the scores of all-against-all protein pairs to the output file.
    :param uniprot_ids: list of all uniprot IDs
    :param output_filename: file to write output to
    :param blast_dict: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values]
    """
    with open(output_filename, "w") as f:
        # Generate all possible combinations of IDs
        combinations_ids = itertools.product(uniprot_ids, uniprot_ids)

        for pair in combinations_ids:
            # Check if the pair has unique elements
            if len(pair) == len(set(pair)):
                pair_str = "\t".join(pair)
                if pair in blast_dict:
                    f.write(pair_str + "\t" + str(blast_dict[pair]) + "\n")
                else:
                    f.write(pair_str + "\t" + "NA\n")


def plot_evalue_distribution(blast_dict, output_figure):
    """
    This function plots the distribution of the log(e-value). pseudocount is added to avoid log(0).
    The pseudocount in this case is the smallest non-zero e-value divided by 1000.
    :param blast_dict: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values]
    :param output_figure: location to save the evalue distribution figure
    """
    sorted_e_val = sorted(blast_dict.values())
    nonzero_indices = numpy.nonzero(sorted_e_val)[0]
    pseudo_count = sorted_e_val[nonzero_indices[0]] / 1000.0
    plt.hist(list(map(lambda x: math.log10(x + pseudo_count), blast_dict.values())))
    plt.xlabel("log(e-value)")
    plt.ylabel("Frequency")
    plt.savefig(output_figure)


    ##########################
    ### START CODING HERE ####
    ##########################
    # (Question 5)
    # Calculate the number of e-values lower than threshold.

    num_evalues_lower_than_threshold = len(list(filter(lambda x: x < 0.002, blast_dict.values())))
    print(f"The number of e-values lower than threshold is: {num_evalues_lower_than_threshold}")

    ##########################
    ###  END CODING HERE  ####
    ##########################


def main():

    parser = argparse.ArgumentParser(description="Automatically running BLAST and PSI-BLAST")
    parser.add_argument("-ids", help="the list of UniProt IDs", required=True)
    parser.add_argument("-q", help="the query folder", required=True)
    parser.add_argument("-db", help="the fasta file of the database", required=True)
    parser.add_argument("-outfile", help="output file", required=True)
    parser.add_argument("-outpng", help="output figure", required=True)
    parser.add_argument("-vblast", choices=['blast', 'psiblast'], default='blast',
                        help="Type of BLAST to be run: blast = BLASTP (default); psiblast = PSI-BLAST")

    args = parser.parse_args()

    # Assign the parsed arguments to the corresponding variables.
    uniprot_id_list = args.ids
    query_folder = args.q
    db = args.db
    version_blast = args.vblast
    output_filename = args.outfile
    output_figure = args.outpng
    
    # The blast_dict dictionary will be used to store protein pair and the corresponding e-value.
    # Keys for blast_dict are the combination of query and subject/hit, e.g.:
    # key             = (query, subject)
    blast_dict = dict()
    # uniprot_ids is a list to store all UniProt IDs contained in uniprot_id_list.
    uniprot_ids = list()

    uniprot_ids_file = open(uniprot_id_list)
    for line in uniprot_ids_file:
        query = line.strip()
        blast_result = blast(db, query, query_folder=query_folder, v_blast=version_blast)
        
        print(blast_result)

        ##########################
        ### START CODING HERE ####
        ##########################
        # Store all the uniprot IDs in the uniprot_ids variable.
        # Hint: the 'query' variable has to be used
        # Furthermore, run the parse_blast_result function with the blast_result
        # and blast_dict variable

        uniprot_ids.append(query)
        #print(line)

        parse_blast_result(blast_result, blast_dict)

        ##########################
        ###  END CODING HERE  ####
        ##########################

    print(uniprot_ids)
    print(len(uniprot_ids))
    uniprot_ids_file.close()
    write_output(uniprot_ids, output_filename, blast_dict)
    plot_evalue_distribution(blast_dict, output_figure)
    print(uniprot_ids)
    print(len(uniprot_ids))
    
if __name__ == "__main__":
    main()
