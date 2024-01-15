import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import os, argparse, sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def calculate_qvals(network, alpha):
    """
    Function used to prepare a matrix of q-values based on the 
    p-value matrix. The q-value matrix is analogous to p-values, 
    but in the cells there are q-values.
    
    :param network: analyzed network
    :param alpha: q-value cut-off threshold
    """

    path = "./dnashape/qvalues/"
    if not os.path.exists(path):
        os.makedirs(path)
    
    data = {}
    
    # p-value table
    df_pvals = pd.read_csv("./dnashape/pvalues/{}.csv".format(network))
    data["Filter"] = list(df_pvals["Unnamed: 0"])
    # dla kazdego parametru obliczanie poprawionych q-wartosci na podstawie 
    # for each parameter, calculating corrected q-values based on p-values for a given parameter for a given filter
    for parameter in list(df_pvals.columns)[1:]: # no column with kmers
        pvals = df_pvals[parameter]
        reject, qvals = fdrcorrection(pvals, alpha=alpha)
        data[parameter] = qvals
    
    df_qvals = pd.DataFrame(data)
    df_qvals = df_qvals.set_index("Filter")
    df_qvals.index.name = None
    df_qvals.to_csv(path+"{}.csv".format(network))
    

def significant_qvals_filters(network, alpha):
    """
    Function used to search for filters that have at least 1 kmer significantly 
    (q-value <= alpha) correlating with each parameter. The function saves to 
    a csv file dataframe in which the index is the filters names and the cells 
    contain values 1 (some k-mer correlates significantly with a given parameter) 
    or 0 (otherwise). Columns mean DNAshape parameters.
    
    :param network: analyzed network
    :param alpha: q-value cut-off threshold
    """
    
    path = "./dnashape/significant/"
    if not os.path.exists(path):
        os.makedirs(path)

    # q-values table
    df_qvals = pd.read_csv("./dnashape/qvalues/{}.csv".format(network))
    # all kmers in file
    kmers = df_qvals["Unnamed: 0"].tolist()
    # all filter numbers to which the kmers belong
    filters = sorted(list(set([x.split("_")[1] for x in kmers])))
    # dataframe which contains information about the filter and its significance (or lack thereof) for each parameter
    df_significant = pd.DataFrame()
    # for each filter
    for filter in filters:
        data = {}
        data["Filter"] = ["filter_{}".format(filter)]
        # q-value dataframe in which the rows refer to kmers from one filter
        sub_df_qvals = df_qvals[df_qvals["Unnamed: 0"].str.contains("filter_{}_".format(filter))]
        for parameter in list(sub_df_qvals.columns)[1:]: # no column with kmers
            qvals = list(sub_df_qvals[parameter])
            # if any q-value is less than or equal to the alpha threshold
            if any(item <= alpha for item in qvals):
                data[parameter] = [1]
            else:
                data[parameter] = [0]
        df_significant_tmp = pd.DataFrame(data)
        df_significant = pd.concat([df_significant,df_significant_tmp])
        
    df_significant = df_significant.set_index("Filter")
    df_significant.index.name = None
    df_significant.to_csv(path+"{}.csv".format(network))    


def main():
    
    parser = argparse.ArgumentParser(description='Calculate q-values and look for filter kmers significantly correlated with each parameter')
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    parser.add_argument('--thresh', action='store', dest='thresh', help='Threshold for q-value cut-off', default="1e-15")
    args = parser.parse_args()
    
    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)
        
    try:
        float(args.thresh)
    except ValueError:
        print("The threshold must be a numeric value that can be converted to float")
        sys.exit(-1)
    
    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # calculate q-values
    print("Calculating q-values...")
    calculate_qvals(args.network, alpha=float(args.thresh))
    print("Done.")
    
    # search for filters kmers
    print("Searching for significantly correlated filters kmers...")
    significant_qvals_filters(args.network, alpha=float(args.thresh))
    print("Done.") 
    
    print("ALL DONE.")



if __name__ == "__main__":
    main()