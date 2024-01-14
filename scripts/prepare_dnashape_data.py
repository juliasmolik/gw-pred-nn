import pandas as pd
import numpy as np
import os, argparse, sys
from itertools import product


def all_repeat(alphabet, n):
    """
    Function used to create all possible sequences of length n.
    
    :param alphabet: the alphabet from which sequences will be created
    :param n: sequence length
    :return: list of all possible sequences of length n
    """
    
    chars = list(alphabet)
    results = []
    for c in product(chars, repeat = n):
        string = "".join(list(c))
        results.append(string)
    results = sorted(list(set(results)))
    return results

def kmers_weights(network, k):
    """
    Function used to create a csv file with information about analyzed one 
    (out of all 1024) sequences and the sum of the weight values from the ("raw") 
    filter. The sum is calculated based on the analyzed k-mer from the filter.
    
    :param network: analyzed network
    :param k: kmer length
    """
    
    nucleotides = {"A":0, "C":1, "G":2, "T":3}
    
    # a directory containing files with the weights of each possible k-mer
    path_weights = "./dnashape/kmer_weights/"
    if not os.path.exists(path_weights):
        os.makedirs(path_weights)

    # a directory containing the sequences of every possible k-mer
    path_seq = "./dnashape/sequences/"
    if not os.path.exists(path_seq):
        os.makedirs(path_seq)

    # filters that have an average weight above the threshold ("significant" filters)
    df_above = pd.read_csv("./statistics/{}_filters_above_tresh.csv".format(network))
    filters_above_tresh = df_above["Filter"].tolist()
    
    # loading the raw filter weight matrix
    filter_file = "./filters/{}_filters.txt".format(network)
    
    # all possible (1024) k-mer sequences
    all_seqs = all_repeat("ACGT", k)

    # saving all possible sequences to separate files 
    for index, seq in enumerate(all_seqs):
        with open(path_seq+"sequences_{}.fa".format(index), "w") as file_seq:
            file_seq.write(">seq_{}\n{}\n".format(index, seq))
    
    data = []

    with open(filter_file, "r") as file:
        lines = file.readlines()
        f = 0
        for i in range(1,len(lines)-1,5):
            # loading raw filter weights
            weights = np.loadtxt(lines[i:i+4]).transpose()
            filter = "filter_{}".format(f)
            # if the filter belongs to the filters above the threshold (there are about 40 of them)
            if filter in filters_above_tresh:
                # looking for kmers
                for i in range(0, 19-k+1):
                    kmer_weights = weights[i:(i+k)]
                    kmer = "filter_{}_{}_{}".format(f, i, i+k)
                    # every possible sequence of length k
                    for seq in all_seqs: 
                        sum_weight = 0
                        # for each nucleotide in the sequence, adding its weight
                        # index = row in filter matrix
                        # [nucleotides[nucleotide]] = column in the filter matrix
                        for index, nucleotide in enumerate(seq):
                            sum_weight += kmer_weights[index][nucleotides[nucleotide]]
                        data.append([seq, kmer, sum_weight])          
            f += 1
    # saving data to a csv file
    df = pd.DataFrame(data, columns=["Seq", "K-mer", "Suma wag z filtra"])
    df.to_csv(path_weights+"{}.csv".format(network), index=False)
    
    
            

def main():

    parser = argparse.ArgumentParser(description='Create all possible sequences of length k. Calculate the convolution of each of them with the weight values in the analyzed filter k-mer')
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    parser.add_argument('--k', action='store', type=int, dest='k', help='Kmer value (length of the divided filter matrix)', default=5)
    args = parser.parse_args()

    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)

    try:
        int(args.k)
    except ValueError:
        print("The k value must be a value that can be converted to int")
        sys.exit(-1)

    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # perform all the analysis
    print("Creating all possible sequences. Calculating the convolution value...")
    kmers_weights(args.network, args.k)
    print("Done.")
    
    print("ALL DONE.")

if __name__ == "__main__":
    main()