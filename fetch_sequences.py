#!/usr/bin/python

# This script downloads all the sequences in the list provided by the user, and puts them into a single file to be used as a database.
# It also stores the fetched fasta sequences into individual files which will be used as queries for (PSI-)BLAST search.

import urllib.request
import argparse
import os.path


def fetch_one_fasta(uniprot_id):
    """
    Fetch the fasta formatted sequence of the uniProtID supplied in the argument.
    :param uniProtID: id of the protein to fetch in uniprot database
    :return: sequence data in fasta format
    """
    ### DONE ###

    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"

    print(f"Fetching fasta repr of uniprot id:{uniprot_id}")
    fh = urllib.request.urlopen(url)
    result = fh.read()
    print(result)
    print()
    fh.close()
    return result


def check_query_folder(query_folder):
    if not os.path.exists(query_folder):
        print("q: " + query_folder)
        print("Query folder does not exist. "
              "Please make sure that the specified folder exists before you run this script.")
        return False
    return True


def fetch_all_sequences(query_folder, uniprot_filename, database_filename):
    """
    Fetch the fasta formatted sequence for each uniProtID in the list.
    :param query_folder: folder containing fasta files.
    :param uniprot_filename: path to the file containing the uniprot ids
    :param database_filename: path to the file saving the database
    """

    if not query_folder.endswith("/"):
        query_folder += "/"

    print("Processing the list of ids...")

    uniprot_file = open(uniprot_filename, "r")
    database_file = open(database_filename, "w")

    for line in uniprot_file:
        ### DONE ###
        # Fetch the fasta formatted sequence for each uniProtID.

        identifier = line.strip()
        fasta_sequence = fetch_one_fasta(identifier).decode('utf-8').strip()

        # Store the fasta sequences as individual fasta file in your query directory.

        with open(f"{query_folder}/{identifier}.fasta", 'w') as f:
            f.write(fasta_sequence)

        # Store all the fasta sequences in one single fasta file as well. These individual files will be used

        database_file.write(fasta_sequence + "\n")

        # as (PSI-)BLAST queries later on.

    print("Processing finished.")
    uniprot_file.close()
    database_file.close()


def main():
    parser = argparse.ArgumentParser(description='This script downloads all the sequences in the list provided'
                                                 'by the user, and puts them into a single file to be used'
                                                 'as a database. It also stores the fetched fasta sequences into'
                                                 'individual files which will be used as queries for (PSI-)BLAST search')

    parser.add_argument("-ids", "--uniprot_filename", help="the list of UniProt IDs", required=True)
    parser.add_argument("-db", "--database_filename", help="Output file containing all fetched sequences",
                        required=True)
    parser.add_argument("-q", "--qfolder", help="Output folder for your queries", required=True)
    args = parser.parse_args()

    query_folder = args.qfolder
    uniprot_filename = args.uniprot_filename
    database_filename = args.database_filename

    if check_query_folder(query_folder):
        fetch_all_sequences(query_folder, uniprot_filename, database_filename)
    else:
        print("Query folder does not exist. "
              "Please make sure that the specified folder exists before you run this script.")

if __name__ == "__main__":
    main()
