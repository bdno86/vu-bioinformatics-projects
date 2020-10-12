#!/usr/bin/python

import os
import urllib.request
import argparse
import itertools
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from shapely.geometry import LineString



def retrieve_go_terms(protein_ids):
    """
    # Get the GO data, from EBI's webservice or from the local file if it was already retrieved and saved.
    :param protein_ids: list of protein IDs.
    :param go_terms_file:
    :return:
    """
    # download GO terms for all uniprot ID
    url = 'http://www.uniprot.org/uniprot/?query=accession:%s&format=tab&columns=id,organism-id,go-id,pathway' % \
          '+OR+'.join(protein_ids)
    link = urllib.request.urlopen(url)
    page = link.read()
    link.close()
    # print(page)
    go_file_name = 'GOterms.tsv'
    out = open('GOterms.tsv', 'wb')
    out.write(page)

    return go_file_name


def get_go_terms_dict(protein_ids, go_data):
    """
    Returns a dictionary, in which each item has a UniProt ID as its key, and
    the corresponding set of GO terms as its value. The code:
    s = go_dict["Q9BS40"]
    would store the set of GO terms associated to protein Q9BS40 in 's'.
    :param go_data: data retrieved from GO database.
    :return: dictionary with UniProt ID as its key, and the corresponding set of GO terms as its value.
    """
    # Define the output dictionary.
    go_dict = {}

    # Loop over the lines in 'go_data' and add every GO annotation to the appropriate set.

    for line in open(go_data, 'r'):
        if len(line) > 0  and not line.startswith('Entry'):
            terms = line.split("\t")

            protein_id = terms[0]
            go_annotations = terms[2].split('; ')
            go_dict[protein_id] = set(go_annotations)

    return go_dict


def read_protein_ids_file(prot_id_file):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    protein_ids = list()
    ##########################
    ### START CODING HERE ####
    ##########################
    # 1: Open the file (prot_id_file)
    # 2: Store all the protein identifiers in the list 'protein_ids'
    # Don't forget to remove all the newlines!
    # for line in open(filename, 'r'):
    #print(type())
    #with open(prot_id_file, 'r') as f:
    for line in prot_id_file.readlines():
        protein_ids.append(line.strip())


    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids


def compute_similarity_score(prot1_go_terms, prot2_go_terms):
    """
    Computes the score for two proteins on the basis of their GO terms
    :param prot1_go_terms:
    :param prot2_go_terms:
    :return:
    """
    intersection_len = len(prot1_go_terms.intersection(prot2_go_terms))
    union_len = len(prot1_go_terms.union(prot2_go_terms))
    score = float(intersection_len) / union_len
    ##########################
    ### START CODING HERE ####
    ##########################
    # Round the score variable to three decimals.
    score = round(score,3)

    ########################
    ### END CODING HERE ####
    ########################
    return score


def check_similarity_for_protein_pair(score, threshold1, threshold2):
    """
    Given a score and two thresholds (threshold1 < threshold2), the following
    function returns the string "different" if the score is less than threshold1;
    it returns the string "ambiguous" if the score is greater than threshold1, but
    less than threshold2; it returns the string "similar" if the score is greater
    than threshold2.
    :param score:
    :param threshold1:
    :param threshold2:
    :return: "different", "similar" or "ambiguous"
    """
    output = str()
    ##########################
    ### START CODING HERE ####
    ##########################
    # This function has three variables; score, threshold1 and threshold2.
    # Your task is to assign the values 'different', 'similar' and 'ambiguous'.
    # If threshold1 is bigger than the score, assign 'different' to output.
    # if the score is between threshold1 and threshold2, assign 'ambiguous' to output.
    # otherwise, assign 'similar' to the output variable.

    score = float(score)

    # between inclusive or exclusive?
    if threshold1 > score:
        output = 'different'
    elif threshold1 <= score <= threshold2:
        output = 'ambiguous'
    else:
        output = 'similar'

    # prob: if th1 > score AND elif th1 > score > th2, if th1==score neither will be TRUE, so similar would be output

    ########################
    ### END CODING HERE ####
    ########################
    return output


def generate_all_possible_protein_pairs(protein_ids):
    """
    Returns a list containing all unique protein pairs.
    :param proteins_list:
    :return:
    """
    pairs = list()
    combinations_ids = itertools.product(protein_ids, protein_ids)
    for pair in combinations_ids:
        # Check if the pair has unique elements
        if len(pair) == len(set(pair)):
            pairs.append(pair)

    return pairs


def assign_homology(go_dict, pairs, threshold1, threshold2):
    """
    :param go_dict:
    :param pairs:
    :param threshold1:
    :param threshold2:
    :return:
    """
    go_homology = {}
    for (p1, p2) in pairs:
        # print(go_dict[p1], go_dict[p2])
        try:
            score = str(compute_similarity_score(go_dict[p1], go_dict[p2]))
            similarity = check_similarity_for_protein_pair(score, threshold1, threshold2)
            #print(f"Similarity: {similarity} for score {score} and th1 {threshold1} and th2 {threshold2}") # EXTRA
            go_homology[(p1, p2)] = [similarity, score]
        except KeyError as e:
            print("No GO terms available for protein", e.args[0])

    return go_homology


def write_results(filename, go_homology):
    """
    Writes the pairs, similarity and, score to the output file.
    :param filename:
    :param go_dict:
    :param pairs_list:
    :param threshold1:
    :param threshold2:
    :return:
    """
    with open(filename, "w") as f:
        for pair, value in go_homology.items():
            try:
                f.write("\t".join(pair) + "\t" + "\t".join(value) + "\n")
            except KeyError as e:
                print("No GO terms available for protein", e.args[0])


def plot_evolution_scores(homology_results):
    """
    Plots the accumulated percentage of sequences (y-axis) with a Jaccard
    score equal or higher than certain value (x-axis).
    :param homology_results:
    :return:
    """
    score = []
    for instance in homology_results.values():
        score.append(float(instance[1]))

    score = np.array(score)
    sorted_scores = np.sort(score)
    # print(np.sort(score))

    accumulated_count = 0
    list_accumulated = []
    unique_sorted = np.unique(sorted_scores)
    for value in unique_sorted:
        accumulated_count += np.sum(sorted_scores == value)
        list_accumulated.append(accumulated_count)

    list_accumulated = np.array(list_accumulated)
    # print(list_accumulated)
    # print(len(sorted_scores))

    norm_accumulated = list_accumulated/float(len(sorted_scores))*100.

    plt.plot(unique_sorted, norm_accumulated)
    plt.xlabel("Jaccard score")
    plt.ylabel("Accumulated % of scores")

    ##########################
    ### START CODING HERE ####
    ##########################
    # Plot two horizontal lines to help you decide suitable threshods
    # as described in the assignment


    #plt.axhline(y=5, color='r', linestyle='-')
    #plt.axhline(y=95, color='b', linestyle='-')

    #first_line = LineString(np.column_stack((unique_sorted, norm_accumulated)))
    #second_line = LineString(np.column_stack((unique_sorted, np.repeat(5, len(score)))))
    #intersection = first_line.intersection(second_line)
    #plt.plot(*intersection.xy, 'o')


    #non_zero = sorted_scores[sorted_scores != 0]
    print(np.percentile(sorted_scores, 90))
    print(np.percentile(sorted_scores, 95))

    #th1 = (1-np.percentile(sorted_scores, 95))*100
    #th2 = (1-np.percentile(sorted_scores, 90))*100
    th1 = np.percentile(sorted_scores, 90)
    th2 = np.percentile(sorted_scores, 95)

    plt.axhline(y=90, color='r', linestyle='-')
    plt.axhline(y=95, color='b', linestyle='-')
    plt.axvline(x=0.054, color='r')
    plt.axvline(x=0.111)
    plt.plot(th1, 90, 'ro')
    plt.plot(th2, 95, 'bo')
    plt.annotate(f"threshold 1: {th1}",xy=(th1, 90), xytext=(0.20 , 65), arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=-90"))
    plt.annotate(f"threshold 2: {th2}",xy=(th2, 95), xytext=(0.30 , 75), arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=-90"))

    ########################
    ### END CODING HERE ####
    ########################

    plt.axis([0, 1 , 0, 100])
    plt.savefig('find_threshold.png')

    # print(np.unique(sorted_scores))
    # print(len(np.unique(sorted_scores)))

    # print(np.sum(sorted_scores == 0.))
    # print(np.sum(sorted_scores == 0.133))

def main():
    parser = argparse.ArgumentParser(description='')

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument("-uniprot_ids", type=argparse.FileType("r"),
                        help="File with the list of protein IDs", required=True)

    required_arguments.add_argument("-output_file", help="", required=True)
    required_arguments.add_argument("-threshold1", help="", required=True)
    required_arguments.add_argument("-threshold2", help="", required=True)

    args = parser.parse_args()

    input_file = args.uniprot_ids
    output_file = args.output_file
    threshold_1 = float(args.threshold1)
    threshold_2 = float(args.threshold2)

    # Parse the input file and retrieve the GO terms associated with each protein.
    protein_ids = read_protein_ids_file(input_file)

    go_data = retrieve_go_terms(protein_ids)
    go_dict = get_go_terms_dict(protein_ids, go_data)

    # Compute and prints the scores and the similarity string of each protein pair.
    pairs = generate_all_possible_protein_pairs(protein_ids)
    go_homology = assign_homology(go_dict, pairs, threshold_1, threshold_2)
    write_results(output_file, go_homology)

    # Plot the evolution of pairs with equal or better Jaccard score
    plot_evolution_scores(go_homology)




if __name__ == "__main__":
    main()