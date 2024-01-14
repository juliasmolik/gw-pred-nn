import numpy as np
import os, csv, argparse, sys

def create_pseudo_counts_matrix(pfm_matrix, epsilon): 
    """
    Function used to create a pseudo counts matrix based on a PFM matrix
    
    :param pfm: PFM matrix
    :param epsilon: pseudo counts value
    :return: pseuto counts matrix
    """
    
    pseudo_counts_matrix = np.empty([pfm_matrix.shape[0], pfm_matrix.shape[1]])
    
    for (x,y), value in np.ndenumerate(pfm_matrix):
        pseudo_counts_matrix[x][y] = (value+(epsilon/4))/(100+epsilon)
    
    return pseudo_counts_matrix


def create_ppm_matrix(pseudo_counts_matrix):
    """
    Function used to create a PPM matrix based on a pseudo counts matrix.

    :param pseudo_counts_matrix: pseudo counts matrix
    :return: PPM matrix
    """
    
    ppm_matrix = np.empty([pseudo_counts_matrix.shape[0], pseudo_counts_matrix.shape[1]])
    for (x,y), value in np.ndenumerate(pseudo_counts_matrix):
        total = np.sum(pseudo_counts_matrix[x])
        ppm_matrix[x][y] = value/total
    return ppm_matrix
    

def ic_content(ppm_matrix):
    """
    Function used to calculate IC content in a PPM matrix.
    
    :param ppm_matrix: PPM matrix
    :return: IC content value for a given matrix
    """
    
    ic_matrix = np.empty([ppm_matrix.shape[0], ppm_matrix.shape[1]])
    
    
    for r in range(len(ppm_matrix)):
        uncertainty = 0
        for value in ppm_matrix[r]:
            uncertainty += value * np.log2(value)
        uncertainty = (-1)*uncertainty
        ic_final = 2 - uncertainty
        row = []
        for value in ppm_matrix[r]:
            row.append(value*ic_final)
        ic_matrix[r] = row
    
    total = np.sum(ic_matrix)
    
    return total
    

def calculate_ic(pfm_file, dataset, network, epsilon, k, ic_tresh):
    """
    Function used to calculate IC content for a given PFM matrix.
    From a PFM matrix it creates a pseudo counts matrix, and then
    a PPM matrix. It divides the PPM matrix intro k-mers matrices
    and for each of them calculates IC value. If IC >= ic_tresh,
    it saves that kmer (PPM matrix) into a file that will be used for
    Tomtom analysis. If IC < ic_tresh, it saves the IC value for that
    kmer into a different file (below_tresh). 
    
    :param pfm_file: txt file with a PFM matrix representing a filter
    :param dataset: class considered 
    :param network: network considered
    :param k: kmer length
    :param ic_tresh: IC value threshold
    """
    
    path_input = "./tomtom/{}/{}/tomtom_input/".format(dataset, network)
    path_kmers = "./tomtom/{}/{}/kmers_below_tresh/".format(dataset, network)  
        
    if not os.path.exists(path_input):
        os.makedirs(path_input)

    if not os.path.exists(path_kmers):
        os.makedirs(path_kmers)
    
    with open(pfm_file) as file:
        lines_tmp = file.readlines()
        lines = lines_tmp[9:]
        for l in lines:
            if "letter-" in l:
                lines.remove(l)
            if l == "\n":
                lines.remove(l)

        for i in range(0,len(lines)-1,20):
            # filter weights
            weights = np.loadtxt(lines[i+1:i+20])
            # filter name
            filter = lines[i].split(" ")[1].strip()
            with open(path_input+filter+".txt", "w") as filterkmers:
                filterkmers.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
                with open(path_kmers+filter+".csv", "w") as belowtresh:
                    writer = csv.writer(belowtresh)
                    writer.writerow(["Filter", "IC"])
                    # looking for kmers
                    for i in range(0, 19-k+1):
                        # kmer weights
                        kmer_freq = weights[i:(i+k)]
                        kmer_pseudo_counts = create_pseudo_counts_matrix(kmer_freq, epsilon)
                        kmer_ppm = create_ppm_matrix(kmer_pseudo_counts)
                        info_content = ic_content(kmer_ppm)
                        if info_content >= ic_tresh:
                            filterkmers.write("MOTIF {}_{}_{}\nletter-probability matrix: alength= 4 w= {}\n".format(filter, i, i+k, k))
                            np.savetxt(filterkmers, kmer_ppm)
                            filterkmers.write("\n")
                        else:
                            writer.writerow(["{}_{}_{}".format(filter, i, i+k), info_content])
                            

def main():

    parser = argparse.ArgumentParser(description='Divide filters into kmers and for each of them create a PFM matrix. Calculate the IC value and save the results to a file.')
    parser.add_argument('--pfm', action='store', dest='pfm', help='The path where the PFM matrices are located', required=True)
    parser.add_argument('--dataset', action='store', dest='dataset', help='The dataset name which the training sequences came from', required=True)
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    parser.add_argument('--epsilon', action='store', type=float, dest='epsilon', help='The pseudo counts value', default=1)
    parser.add_argument('--k', action='store', type=int, dest='k', help='Kmer value (length of the divided filter matrix)', default=7)
    parser.add_argument('--thresh', action='store', type=float, dest='threshold', help='IC value cutoff threshold', default=6.5)
    args = parser.parse_args()

    if args.pfm == None: 
        print("No PFM file provided (see help; -h)")
        sys.exit(-1)

    if args.dataset == None: 
        print("No dataset name provided (see help; -h)")
        sys.exit(-1)

    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)

    try:
        float(args.epsilon)
    except ValueError:
        print("The epsilon value must be a value that can be converted to float")
        sys.exit(-1)
        
    try:
        int(args.k)
    except ValueError:
        print("The k value must be a value that can be converted to int")
        sys.exit(-1)

    try:
        float(args.threshold)
    except ValueError:
        print("The threshold value must be a value that can be converted to float")
        sys.exit(-1)
    
    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # perform all the calculations
    print("Dividing filters into kmers, creating a PPM matrix, and calculating the IC value...")
    calculate_ic(args.pfm, args.dataset, args.network, args.epsilon, args.k, args.threshold)
    print("Done.")
    
    print("ALL DONE.")
    

if __name__ == "__main__":
    main()