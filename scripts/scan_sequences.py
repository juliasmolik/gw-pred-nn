from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import numpy as np
import pandas as pd
import os, argparse, sys


def read_sequences(file, filetype = "fasta"):
    """
    Function used to read sequences from a file.

    :param file: file with sequences
    :param filetype: file type; default: fasta
    :return: list of sequences from the file
    """
    
    sequences = []
    
    for seq_record in SeqIO.parse(file, filetype):
        sequences.append(seq_record.seq.upper())
    
    return sequences


def one_hot_encode(seq):
    """
    Function used to encode nucleotide sequences using one-hot encoding.

    :param seq: sequence to be encoded
    :return: one-hot encoded sequence
    """
    mapping = dict(zip("ACGT", range(4)))    
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]


def scan_sequences(sequences_file, filters_file):    
    """
    Function used to scan a sequence with a filter. It calculates a convolution
    of filter weights with a one-hot encoded matrix representing a sequence.
    It returns a dictionary, the key of which is the name of a filter and the value
    is a list of the best sub-sequences, where each of them comes from a different 
    sequence (for each sequence the function chooses the best sub-sequence). 
    
    :param sequences_file: file with sequences from training dataset
    :param filters_file: file with filter weights
    :return: dictionary, the key of which is the name of a filter and the value
    is a list of the best sub-sequences, where each of them comes from a different 
    sequence
    """
    
    # dictionary in which the key is the names of subsequent filters,
    # and the value is a list of tuples, where at zero position there
    # is the best sub-sequence for each sequence, and at first position
    # there is the convolution value of this subsequence and the filter
    best_sequences = {}
    
    sequences = read_sequences(sequences_file)
    
    network_name = os.path.basename(filters_file).replace("_filters.txt", "")
    df_above = pd.read_csv("./statistics/{}_filters_above_tresh.csv".format(network_name))
    filters_above_tresh = df_above["Filter"].tolist()
    
    nucleotides = {0: "A", 1: "C", 2: "G", 3: "T"}

    
    for sequence in sequences:
        
        # one-hot encoded sequence
        matrix = one_hot_encode(sequence).transpose()
        
        with open(filters_file, "r") as file:
            lines = file.readlines()
            f = 0
            for i in range(1,len(lines)-1,5):
                # load the filter
                filter = np.loadtxt(lines[i:i+4])
                best = float('-inf')
                if "filter_{}".format(f) in filters_above_tresh:
                    # scanning the matrix of a given sequence with a filter of size (4,19)
                    for j in range(np.shape(matrix)[1]):
                        if np.shape(matrix[:,j:j+np.shape(filter)[1]])[1] == np.shape(filter)[1]:
                            # sequence submatrix of size (4,19)
                            submatrix = matrix[:,j:j+np.shape(filter)[1]]
                            # checking which subsequence multiplied by the filter 
                            # values will give the highest value
                            convolution = filter*submatrix
                            total = np.sum(convolution)
                            if total > best:
                                best = total
                                best_seq = submatrix
                                where = np.argwhere(best_seq)
                                seq = "*" * np.shape(best_seq)[1]
                                # sequence reconstruction based on one hot encoding matrix
                                for index in where:
                                    seq = seq[:index[1]] + nucleotides[index[0]] + seq[index[1]+1:]
                    if "filter_"+str(f) in best_sequences.keys():
                        tmp = best_sequences["filter_"+str(f)]
                        tmp.append((seq, float(best)))
                    else:
                        best_sequences["filter_"+str(f)] = [(seq, float(best))]
                f += 1
    return best_sequences


def choose_best_subseq(dictionary, n):
    """
    Function used to find n best sub-sequences in the dictionary,
    in which the key is the filter name and the value is a list of
    tuples, where at zero position there is the best sub-sequence
    for each sequence that has been scanned with this filter, and
    at the first position there is the convolution value of this 
    sub-sequence and the filter. The function from all found
    sub-sequences, chooses n best sub-sequences that will be
    used to create a PFM matrix. The function returns a dictionary,
    the key of which is the filter name, and the value is a list
    of n best sub-sequences. 
    
    :param dictionary: dictionary, the key of which is the filter
    name, and the value is a list of the best sub-sequences along
    with their convolution value with this filter
    :param n: number of the best sub-sequences; default: 100
    :return: dictionary, the key of which is the filter name, and 
    the value is a list of n best sub-sequences
    """

    dictionary_new = {}
    for key in dictionary:
        sorted_list = sorted(dictionary[key], key=lambda tup: tup[1], reverse=True)
        tmp = []
        try:
            for i in range(n):
                tmp.append(sorted_list[i][0])
        except:
            print("The chosen n value is too big. Setting the value: {}".format(len(sorted_list)))
            for i in range(len(sorted_list)):
                tmp.append(sorted_list[i][0])
        dictionary_new[key]=tmp

    return dictionary_new


def create_pfm_matrices(dictionary, dataset, network):
    """
    Function used to create a PFM matrix based on a dictionary
    the key of which is the filter name, and the value is a list
    of n best sub-sequences. It saves the created PFM matrices
    for each filter separately to a file. 
    
    :param dictionary: dictionary the key of which is the filter 
    name, and the value is a list of n best sub-sequences
    """
    
    path = "./pfm_filters/{}/".format(dataset)
    if not os.path.exists(path):
        os.makedirs(path)
    
    with open(path+"{}_pfm.txt".format(network), "w") as fileout:
        fileout.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        # for each filter
        for filter in dictionary:
            fileout.write("MOTIF {}\nletter-probability matrix: alength= 4 w= 19\n".format(str(filter)))
            instances = []
            # creating Seq instance for each sub-sequence
            for seq in dictionary[filter]:
                instances.append(Seq(seq))
            # creating a PFM matrix for the set of sub-sequences for each filter
            m = motifs.create(instances)
            pfm_matrix = m.counts
            # replace Motif instance with numpy matrix to change
            # dimentions of the pfm matrix to be 19x4 and not 4x19
            pfm_list = []
            for nucleotide in ["A", "C", "G", "T"]:
                pfm_list.append(list(pfm_matrix[nucleotide]))
            pfm_list = np.matrix(pfm_list).transpose()
            np.savetxt(fileout, pfm_list)
            fileout.write("\n")

    
    

def main():
    
    parser = argparse.ArgumentParser(description='Scan training dataset sequences using chosen filters and create PFM matrices based on n best sub-sequences for each filter')
    parser.add_argument('--sequences', action='store', dest='seq', help='The path where the training dataset sequences are located', required=True)
    parser.add_argument('--filters', action='store', dest='filters', help='The path where the files with filters weights are located', required=True)
    parser.add_argument('--dataset', action='store', dest='dataset', help='The dataset name from which the training sequences came', required=True)
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    parser.add_argument('--n', action='store', dest='number', type= int, help='Number of sub-sequences used to create PFM matrices; default=100', default=100)
    args = parser.parse_args()

    if args.seq == None: 
        print("No sequences file provided (see help; -h)")
        sys.exit(-1)

    if args.filters == None: 
        print("No filters files provided (see help; -h)")
        sys.exit(-1)

    if args.dataset == None: 
        print("No dataset name provided (see help; -h)")
        sys.exit(-1)

    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)
        
    try:
        int(args.number)
    except ValueError:
        print("The number must be a value that can be converted to int")
        sys.exit(-1)
    
    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # scan sequences
    print("Scanning sequences...")
    scanned_seqs = scan_sequences(args.seq, args.filters)
    print("Done.")
    
    # choose best sub-sequences
    print("Choosing {} best sequences...".format(str(args.number)))
    best_subseq = choose_best_subseq(scanned_seqs, args.number)
    print("Done.")
    
    
    # create PFM matrices
    print("Creating PFM matrices...")
    create_pfm_matrices(best_subseq, args.dataset, args.network)
    print("Done.")
    
    print("ALL DONE.")
    

if __name__ == "__main__":
    main()