import torch
from torch import nn
import math, os, csv, argparse, sys
import numpy as np
import pandas as pd
from pathlib import Path

class CustomNetwork(torch.nn.Module):

    def __init__(self, seq_len=2000, num_channels=[300, 200, 200], kernel_widths=[19, 11, 7], pooling_widths=[3, 4, 4],
                 num_units=[2000, 4], dropout=0.5):
        super(CustomNetwork, self).__init__()
        paddings = [int((w-1)/2) for w in kernel_widths]
        self.seq_len = seq_len
        self.dropout = dropout
        self.params = {
            'input sequence length': seq_len,
            'convolutional layers': len(num_channels),
            'fully connected': len(num_units),
            'number of channels': num_channels,
            'kernels widths': kernel_widths,
            'pooling widths': pooling_widths,
            'units in fc': num_units,
            'dropout': dropout

        }

        conv_modules = []
        num_channels = [1] + num_channels
        for num, (input_channels, output_channels, kernel, padding, pooling) in \
                enumerate(zip(num_channels[:-1], num_channels[1:], kernel_widths, paddings, pooling_widths)):
            k = 4 if num == 0 else 1
            conv_modules += [
                nn.Conv2d(input_channels, output_channels, kernel_size=(k, kernel), padding=(0, padding)),
                nn.BatchNorm2d(output_channels),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(1, pooling), ceil_mode=True)
            ]
            seq_len = math.ceil(seq_len / pooling)
        self.conv_layers = nn.Sequential(*conv_modules)

        fc_modules = []
        self.fc_input = 1 * seq_len * num_channels[-1]
        num_units = [self.fc_input] + num_units
        for input_units, output_units in zip(num_units[:-1], num_units[1:]):
            fc_modules += [
                nn.Linear(in_features=input_units, out_features=output_units),
                nn.ReLU(),
                nn.Dropout(p=self.dropout)
            ]
        self.fc_layers = nn.Sequential(*fc_modules)

    def forward(self, x):
        x = self.conv_layers(x)
        x = x.view(-1, self.fc_input) 
        x = self.fc_layers(x)
        return torch.sigmoid(x)


def get_filters(modelfile, fileout):
    """
    Function to extract filters from the first network layer
    and saving them to a file.
    
    :param modelfile: trained model file (network)
    :param fileout: file to which the filter weights are written
    """
    
    path = './filters/'
    if not os.path.exists(path):
        os.makedirs(path)
    
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    model = CustomNetwork()
    # loading the trained model
    model.load_state_dict(torch.load(modelfile, map_location=torch.device(device)))
    for name, param in model.named_parameters():
        # zero (first) convolutional layer
        if name == "conv_layers.0.weight":
            with open(path+fileout, "w") as f:
                for i in range(len(param)):
                    f.write(">filter_"+str(i)+"\n")
                    # each filter has a size of 4x19 and there are 300 of them (from 0 to 299)
                    filter = param[i].detach().numpy().reshape((4,19))
                    np.savetxt(f,filter)


def calculate_statistics_filters(modelfile):
    """
    Function that calculates filters statistics. Creates a csv file, which 
    contains information about the max, min, mean, std and median value 
    of each filter.
    
    :param modelfile: trained model file (network)
    """
    
    path = "./statistics/"
    if not os.path.exists(path):
        os.makedirs(path)
        
    with open("./filters/{}_filters.txt".format(modelfile), "r") as file:                 
        lines = file.readlines()
        f = 0
        # creating the filename
        filename = "{}_stats.csv".format(modelfile)
        with open(path+filename, 'w') as fileout:
            writer = csv.writer(fileout)
            writer.writerow(["Filter", "Max", "Min", "Mean", "Std", "Median"])
            for i in range(1,len(lines)-1,5):
                # loading a filter
                filter = np.loadtxt(lines[i:i+4])
                max = np.max(filter)
                min = np.min(filter)
                mean = np.mean(filter)
                std = np.std(filter)
                median = np.median(filter)
                fil = "filter_{}".format(f)
                writer.writerow([fil, max, min, mean, std, median])
                f += 1

def divide_filters_weights(filter_stats, thresh):
    """
    Function used to divide a set of filters into those that
    have average values of weights above the threshold and those that have
    these values below. They are saved to separate files:
    {modelfile}_filter_above_tresh.csv and {modelfile}_filter_below_tresh.csv respectively.

    :param filter_stats: file with statistics for filters in a given network
    :param tresh: filter average value threshold
    """

    path = "./statistics/"
    if not os.path.exists(path):
        os.makedirs(path)
    
    stats = pd.read_csv(filter_stats)
    
    below_tresh = []
    above_tresh = []
    
    for index, row in stats.iterrows():
        # if the average value in the filter is above the threshold
        if abs(row["Mean"]) > thresh:
            above_tresh.append([row["Filter"], row["Mean"]])
        # otherwise
        else:
            below_tresh.append([row["Filter"], row["Mean"]])
    
    df_below = pd.DataFrame(below_tresh, columns=["Filter", "Mean"])
    df_above = pd.DataFrame(above_tresh, columns=["Filter", "Mean"])
    
    filter = os.path.basename(filter_stats).replace("_stats.csv", "")
    
    df_below.to_csv(path+filter+"_filters_below_tresh.csv")
    df_above.to_csv(path+filter+"_filters_above_tresh.csv")


def main():
    
    
    parser = argparse.ArgumentParser(description='Extract filters from the network, calculate their stats and then based on that divide them into those with high and low average weight values')
    parser.add_argument('--model', action='store', dest='model', help='The path where the network file is located', required=True)
    parser.add_argument('--thresh', action='store', dest='thresh', help='The threshold above which the average filter weight is considered high', default="1e-5")
    args = parser.parse_args()

    if args.model == None: 
        print("No model file provided (see help; -h)")
        sys.exit(-1)
    else:
        if not args.model.endswith('.model'):
            print("Invalid model format")
            sys.exit(-1)
    
    try:
        float(args.thresh)
    except ValueError:
        print("The threshold must be a numeric value that can be converted to float")
        sys.exit(-1)
    
    # get filters from the network
    get_filters(args.model, "{}_filters.txt".format(Path(args.model).stem))
    
    # calculate filter statistics
    calculate_statistics_filters(Path(args.model).stem)
    
    # get filters the weights of which meet threshold
    divide_filters_weights("./statistics/{}_stats.csv".format(Path(args.model).stem), float(args.thresh))

if __name__ == "__main__":
    main()