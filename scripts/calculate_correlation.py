import pandas as pd
import scipy.stats, os, argparse, sys

def calculate_correlation(network):
    """
    Function used to create a correlation matrix in which
    columns contain DNAshape parameters (7 (5+2) parameters), and rows
    these are the names of the filter kmers.
    
    :param network: analyzed network
    """

    path_corr = "./dnashape/correlation/"
    if not os.path.exists(path_corr):
        os.makedirs(path_corr)

    path_pval = "./dnashape/pvalues/"
    if not os.path.exists(path_pval):
        os.makedirs(path_pval)
    
    data_corr = []
    data_pvals = []
    
    # loading weights for filter k-mers
    kmer_weights = pd.read_csv("./dnashape/kmer_weights/{}.csv".format(network))
    kmers = sorted(list(set(kmer_weights["K-mer"])))
    
    
    for kmer in kmers:
        # only those lines that concern the analyzed kmer
        sub_df = kmer_weights[kmer_weights["K-mer"] == kmer]
        # all possible k-mers
        sequences = sorted(list(set(sub_df["Seq"])))
        parameters = {"MGW": [], "HelT1": [], "HelT2": [], "ProT": [], "Roll1":[], "Roll2":[], "EP": []}
        weights = []
        for index, seq in enumerate(sequences):
            # a row in which a given kmer and sequence are analyzed (weight of a given sequence based on a given kmer)
            sub_weights = kmer_weights[(kmer_weights["K-mer"] == kmer) & (kmer_weights["Seq"] == seq)]
            # loading DNAshape parameters for a given sequence
            params = pd.read_csv("./dnashape/dna_parameters/sequences_{}.csv".format(index))
            # adding k-mer weight to the list
            weights.append(sub_weights.values[0][2])
            for idx, row in params.iterrows():
                if row["Parameter"] == "HelT":
                    if row["Vector"] == "V2":
                        parameters["HelT1"].append(row["Value"])
                    else:
                        parameters["HelT2"].append(row["Value"])
                elif row["Parameter"] == "Roll":
                    if row["Vector"] == "V2":
                        parameters["Roll1"].append(row["Value"])
                    else:
                        parameters["Roll2"].append(row["Value"])
                else:
                    parameters[row["Parameter"]].append(row["Value"])
            
        # a list to store all 7 (5+2) correlations
        data_tmp_corr = [] 
        # a list to store all 7 (5+2) p-values
        data_tmp_pval = []
        for parameter in parameters:
            spearman = scipy.stats.spearmanr(parameters[parameter], weights)
            corr = spearman[0]
            pval = spearman[1]
            data_tmp_corr.append(corr)
            data_tmp_pval.append(pval)
        
        data_corr.append(data_tmp_corr)
        data_pvals.append(data_tmp_pval)
        
    # dataframe in which the row names are the names of k-mers, and the column names 
    # are the names of the parameters, calculated correlations are entered in the cells
    df_corr = pd.DataFrame(data_corr, columns = list(parameters.keys()), index=kmers)
    df_corr.to_csv(path_corr+"{}.csv".format(network))
    
    # dataframe in which the row names are the names of k-mers, and the column names 
    # are the names of the parameters, calculated p-values are entered in the cells
    df_pvals = pd.DataFrame(data_pvals, columns = list(parameters.keys()), index=kmers)
    df_pvals.to_csv(path_pval+"{}.csv".format(network))


def main():
    
    parser = argparse.ArgumentParser(description='Calculate the correlation of the convolution value depending on the analyzed sequence with the DNAshape parameters of a given sequence')
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    args = parser.parse_args()
    
    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)

    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # perform all the analysis
    print("Calculating correlation...")
    calculate_correlation(args.network)
    print("Done.")
    
    print("ALL DONE.")
    

if __name__ == "__main__":
    main()