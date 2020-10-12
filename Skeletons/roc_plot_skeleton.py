#!/usr/bin/python

import argparse
import math
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from operator import itemgetter
import pylab


def calculate_rate(tpr, fpr, num_similar, num_different,
                   total_similar, total_different):
    """
    Calculates the True Positive Rate (TPR) and False Positive Rate (FPR)
    based on the number of similar and different values.
    :param tpr: (list) TPR rates
    :param fpr: (list) FPR rates
    :param num_similar: (integer) Number of similar values
    :param num_different: (integer) Number of different values
    :param total_similar: (integer) Total number of similar values
    :param total_different: (integer) Total number of different values
    :return: Nothing
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # Calculate the TPR and FPR and append the value to the tpr and fpr
    # variables.

    ##########################
    ###  END CODING HERE  ####
    ##########################


def plot_figure(tpr, fpr, evals, file_name):
    """
    Plot the ROC-plot based on the
    :param tpr: (list) all the TPR values
    :param fpr: (list) all the FPR values
    :param evals: (list) all the sorted e-values
    :return:
    """
    if len(tpr) < 2 or len(fpr) < 2:
        return
    pylab.figure()
    lw = 2
    pylab.scatter(tpr, fpr, cmap=plt.cm.coolwarm, s=20,
                  lw=lw, c=evals, edgecolors='none')
    pylab.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--', label='Random')
    pylab.xlim([0.0, 1.0])
    pylab.ylim([0.0, 1.0])
    pylab.xlabel('False Positive Rate')
    pylab.ylabel('True Positive Rate')
    pylab.title('Receiver operating characteristic')
    pylab.legend(loc='lower right')
    color_bar = pylab.colorbar()
    color_bar.ax.set_ylabel('E-value score')
    pylab.savefig(fname = file_name)


def count_total_results(blast_results, go_results):
    """
    Return the total number (PSI-)BLAST and GO results.
    :param blast_results: (list) all the (PSI-)BLAST results
    :param go_results: (list)  all the GO results
    :return: similar; total similar results. different; total different results
    """
    similar = 0
    different = 0
    for key, value in blast_results.items():
        if value != 'NA':
            if key not in go_results:
                protein1, protein2 = key.split('_')
                new_key = protein2 + '_' + protein2
                if new_key in go_results:
                    key = new_key
                else:
                    print('No GO results for this combination: ' + key)
                    continue
            result = go_results[key]
            if result == 'similar':
                similar += 1
            elif result == 'different':
                different += 1
    return [similar, different]


def get_roc_values(blast_results, go_results):
    """
    Function for calculating the TPR and FPR for the x-y axis ROC-plot
    :param blast_results: (list) all the (PSI-)BLAST results
    :param go_results: (list) all the GO results
    :return:
    FPR & TPR list of sorted False positive rates
    evals, list of all the sorted evalues
    """
    tpr = [0.0]
    fpr = [0.0]
    evals = [min(blast_results.values())]
    num_different = 0.0
    num_similar = 0.0
    total_similar, total_different = count_total_results(blast_results, go_results)
    print('number of similar GO results: ' + str(total_similar))
    print('number of different GO results: ' + str(total_different))
    sorted_blast_scores = sorted(blast_results.items(), key=itemgetter(1))

    for i, (protein_comb, score) in enumerate(sorted_blast_scores):
        if protein_comb in go_results:
            result = go_results[protein_comb]

            if result == 'different':
                num_different += 1.0
                if num_different == 1.0:
                    print(protein_comb, score, result)
            elif result == 'similar':
                num_similar += 1.0

            calculate_rate(tpr, fpr, num_similar, num_different,
                           total_similar, total_different)

            evals.append(score)

    return tpr, fpr, evals


def parse_benchmark_results(go_file):
    """
    Put all the GO results into a dictionary.
    :param go_file: (File object) go results
    :return: benchmark_dict: (dict) key: protein combinations, value: go result
    """
    benchmark_dict = dict()

    for line in go_file:
        splitted_line = line.rstrip().split("\t")[0:3]
        protein1, protein2, result = splitted_line
        if protein1 != protein2:
            key = protein1 + "_" + protein2
            if key not in benchmark_dict:
                benchmark_dict[key] = result

    go_file.close()

    return benchmark_dict


def parse_blast_results(blast_file):
    """
    Put all the (PSI-)BLAST results into a dictionary.
    :param blast_file: (File object) All the (PSI-)BLAST results
    :return: blast_evalues: (dict) dict with all the blast e-values as value.
    """
    blast_evalues = dict()

    for line in blast_file:
        if line.strip():  # If line is not empty

            splitted_line = line.rstrip().split("\t")[0:3]
            protein1, protein2, result = splitted_line
            if protein1 != protein2:
                key = protein1 + "_" + protein2
                duplicate = protein2 + "_" + protein1
                if key not in blast_evalues and duplicate not in blast_evalues:
                    if result != "NA":
                        try:
                            blast_evalues[key] = math.log(float(result))
                        except ValueError:
                            pass

    blast_file.close()

    return blast_evalues


def main():
    """
    Main function
    parses arguments:
    -blast_results, tab-separated file with e-values
    -go_results, tab-separated file with GO results
    -outpng, path to ROC plot figure
    """
    parser = argparse.ArgumentParser(description='')

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument("-blast_results", type=argparse.FileType("r"),
                                    help="File with all the BLAST results", required=True)
    required_arguments.add_argument("-go_results", help="GO results file", required=True,
                                    type=argparse.FileType("r"))
    required_arguments.add_argument("-outpng", help="ROC plot output", required=False)

    args = parser.parse_args()
    blast_results_file = args.blast_results
    go_results_file = args.go_results
    roc_plot = args.outpng

    blast_evalues = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(go_results_file)
    tpr, fpr, evals = get_roc_values(blast_evalues, benchmark_results)

    plot_figure(fpr, tpr, evals, roc_plot)


# initialise main function
if __name__ == "__main__":
    main()