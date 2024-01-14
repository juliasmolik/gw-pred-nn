import os, subprocess, csv, shutil, argparse, sys
import pandas as pd
from pathlib import Path

def tomtom(kmers_files, dataset, network):
    """
    Function used to perform Tomtom analysis. It loads the
    file with kmers from PPM matrices and saves the Tomtom
    results into folders corresponding to the filters.
    
    :param kmers_files: file with kmers from PPM matrices
    :param dataset: class analyzed 
    :param network: network analyzed
    """

    path = "./tomtom/{}/{}/tomtom_output/".format(dataset, network) 

    if not os.path.exists(path):
        os.makedirs(path)
    else:
        shutil.rmtree(path)
        os.makedirs(path)
    
    for entry in os.scandir(kmers_files):
        if entry.path.endswith(".txt") and entry.is_file():
            filter = os.path.basename(entry.path).replace(".txt", "")
            subprocess.run("tomtom -incomplete-scores -dist kullback -o {} {} human_hocomoco.meme".format(path+filter, entry.path), shell=True)

def count_found_motifs(dataset, network, k=7):
    """
    Function used to create csv files containing the analysis of Tomtom results.
    It counts how many TF motifs were found for each kmer and creates a file with 
    kmer in the rows and in the columns names of all TFs. The intersection is 1 when 
    a given TF was a hit for this one kmer, and 0 otherwise.   
    
    :param dataset: class analyzed 
    :param network: network analyzed
    :param k: kmer length
    """

    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    path2 = "./tomtom/{}/{}/tomtom_analysis/count_tf_hits/".format(dataset, network)
    
    if not os.path.exists(path2):
        os.makedirs(path2)

    # list of all found TFs for all kmers in each filter
    tf_names = []
    
    # dictionary in which the key is the name of the kmer and the value 
    # is a list of found TF motifs for that kmer
    dictionary = {}
    # iterating directories corresponding to filters
    for subdir, dirs, files in os.walk("./tomtom/{}/{}/tomtom_output/".format(dataset, network)):
        for file in files:
            if file.endswith(".tsv"):
                df = pd.read_csv(os.path.join(subdir, file), sep='\t')
                for index, row in df.iterrows():
                    if "H11" in str(row["Target_ID"]):
                        if row["Query_ID"] not in dictionary:
                            dictionary[row["Query_ID"]] = [row["Target_ID"]]
                        else:
                            dictionary[row["Query_ID"]].append(row["Target_ID"])
                            
                # creating TF found count files for each kmer
                with open(path2+Path(subdir).stem+".csv", "w") as count_hits:
                    writer = csv.writer(count_hits)
                    writer.writerow(["Filter kmer", "Number of TF motifs found"])
                    for i in range(0, 19-k+1):
                        kmer = "{}_{}_{}".format(Path(subdir).stem, i, i+k)
                        if kmer in list(dictionary.keys()):
                            writer.writerow([kmer, len(dictionary[kmer])])
                        else:
                            writer.writerow([kmer, 0])
                            
                for key in dictionary:
                    for value in dictionary[key]:
                        if value not in tf_names:
                            tf_names.append(value)
    # creating a file showing which TF motifs occurred 
    # in a given kmer (1 if true, 0 otherwise)    
    with open(path+ "tf_motifs_in_kmers.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [i for i in tf_names]
        row.insert(0, "K-mer")
        writer.writerow(row)
        for key in dictionary:
            row2 = [key]
            for r in row:
                if r != "K-mer":
                    if r in dictionary[key]:
                        row2.append(1)
                    else:
                        row2.append(0)
            writer.writerow(row2)
                        

def aggregate_filters(dataset, network):
    """
    Function used to aggregate Tomtom results into filters.
    It gets hits for kmer csv file and creates a new csv file,
    in which the data is aggregated into filters.
    
    :param dataset: class analyzed 
    :param network: network analyzed
    """
    
    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    # loading a file with the kmer in the rows and the names of all TFs in the columns.
    # at the intersection there is 1 if a given TF was a hit for this kmer, and 0 otherwise
    df = pd.read_csv("./tomtom/{}/{}/tomtom_analysis/tf_motifs_in_kmers.csv".format(dataset, network))
    
    # a dictionary in which the key is the name of the filter and the value is a list of TF motifs found for it
    dictionary = {}
    # TF names from the kmer file
    columns = list(df.keys())
    for i in range(300):
        # analyzed filter
        filter = "filter_{}".format(i)
        # for each kmer from the file
        for index, row in df.iterrows():
            # loading the kmer name from a file
            filter_tmp = row['K-mer'].split("_")
            # for which filter was this kmer
            filter_tmp_2 = "_".join(filter_tmp[:2])
            # if the analyzed filter was in the kmer file
            if filter == filter_tmp_2:
                # checking which TF motifs found for this filter's kmer
                indices = [i for i, x in enumerate(row) if x == 1]
                # adding the TF names found for the analyzed filter
                if filter in dictionary:
                    for ind in indices:
                        dictionary[filter].append(columns[ind])
                else:
                    dictionary[filter] = [columns[ind] for ind in indices]
    # avoiding repetition of TF names for the filter
    for key in dictionary:
        dictionary[key] = list(set(dictionary[key]))
    

    # saving to file with aggregation to filter
    with open(path+"tf_motifs_in_filters.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [col for col in columns]
        row[0] = "Filter"
        row.insert(len(row), "Sum")
        writer.writerow(row)
        for key in dictionary:
            # for each filter, creating 1s and 0s corresponding to the occurrences of TF motifs
            values = []
            # dla kazdego TF z pliku z kmerami
            for col in columns:
                if "K-mer" != col:
                    # if a given TF motif was found for the kmer of a given filter
                    if col in dictionary[key]:
                        values.append(1)
                    # otherwise
                    else:
                        values.append(0)
            suma = sum(values)
            values.append(suma)
            values.insert(0, key)
            writer.writerow(values)

def annotate_filters(dataset, network):
    """
    Function used to aggregate Tomtom results into families.
    It loads hits for a csv file from aggregation to filters and creates
    a new csv file in which the data is aggregated into families.
    
    :param dataset: class analyzed 
    :param network: network analyzed
    """

    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    # dictionary in which the key is the name TF and the 
    # value is the name of the family to which it belongs
    dict_annotate = {}
    df = pd.read_csv("./hocomoco_annotate.tsv", sep="\t")
    for index, row in df.iterrows():
        # for some TFs there is no information about the family
        if type(row["TF family"]) == str:
            dict_annotate[row["Model"]] = row["TF family"]
        else:
            dict_annotate[row["Model"]] = "No TF family info for this motif"
    
    vals = list(dict_annotate.values())
    for v in vals:
        if type(v) != str:
            print(v)
    
    # loading a file with aggregation into filters
    df = pd.read_csv(path+"tf_motifs_in_filters.csv")
    df = df.drop(columns=['Sum'])
    # dictionary in which the key is the filter name and the value 
    # is a list of TF families that were found for this filter
    dictionary = {}
    # TF names present in the filter aggregation file
    columns = list(df.keys())
    # list of family names of these TFs that were in the filter aggregation file
    columns_annotated = [dict_annotate[col] for col in columns if col != "Filter"]
    columns_annotated = list(set(columns_annotated))
    # for each filter
    for index, row in df.iterrows():
        # which TF motifs found for this filter
        indices = [i for i, x in enumerate(row) if x == 1]
        # adding information about TFs families for each filter to the dictionary
        if row["Filter"] in dictionary:
            for ind in indices:
                dictionary[row["Filter"]].append(dict_annotate[columns[ind]])
        else:
            dictionary[row["Filter"]] = [dict_annotate[columns[ind]] for ind in indices]
            
    # removing repetitions of family names for a given filter
    for key in dictionary:
        dictionary[key] = list(set(dictionary[key]))
    # saving to a file with aggregation to families
    with open(path+"tf_families_in_filters.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [col for col in columns_annotated]
        row.insert(0, "Filter")
        row.insert(len(row), "Sum")
        writer.writerow(row)
        for key in dictionary:
            values = []
            for col in columns_annotated:
                if col != "Filter":
                    if col in dictionary[key]:
                        values.append(1)
                    else:
                        values.append(0)
            suma = sum(values)
            values.append(suma)
            values.insert(0, key)
            writer.writerow(values)


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--ppm', action='store', dest='ppm', help='The path where the filter kmer PPM matrices are located', required=True)
    parser.add_argument('--dataset', action='store', dest='dataset', help='The dataset name which the training sequences came from', required=True)
    parser.add_argument('--network', action='store', dest='network', help='The network name where the filters came from', required=True)
    args = parser.parse_args()
    
    if args.ppm == None: 
        print("No PFM file provided (see help; -h)")
        sys.exit(-1)

    if args.dataset == None: 
        print("No dataset name provided (see help; -h)")
        sys.exit(-1)

    if args.network == None: 
        print("No network name provided (see help; -h)")
        sys.exit(-1)
    
    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # perform Tomtom analysis
    print("Performing Tomtom analysis...")
    tomtom(args.ppm, args.dataset, args.network)
    print("Done.")
    
    # hit counts and summary of found TFs for kmers
    print("Analysing Tomtom results...")
    count_found_motifs(args.dataset, args.network)
    aggregate_filters(args.dataset, args.network)
    annotate_filters(args.dataset, args.network)
    print("Done.")
    
    print("ALL DONE.")


if __name__ == "__main__":
    main()