import torch
from torch import nn
import math, os, random, csv, sys, argparse, glob 
import pandas as pd
from Bio import SeqIO
from sklearn.preprocessing import OneHotEncoder as Encoder
from torch.utils.data import Dataset
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from statistics import mean
from itertools import product

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


def change_weights(modelfile):
    """
    Function used to zero the weights in the first layer filters,
    which have an average weight value below the preset threshold.

    :param modelfile: trained network in which filters needs to be changed
    """
    
    path = "./changed_models/"
    if not os.path.exists(path):
        os.makedirs(path)
    
    network_name = Path(modelfile).stem
        
    # load filters that have an average weight below the threshold
    df = pd.read_csv("./statistics/{}_filters_below_tresh.csv".format(network_name))
    below_tresh = df["Filter"].tolist()
    
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    model = CustomNetwork()
    # load the trained model
    model.load_state_dict(torch.load(modelfile, map_location=torch.device(device)))
    for name, param in model.named_parameters():
        # first convolutional layer
        if name == "conv_layers.0.weight":
            weights = param.detach()

    # change filters weights 
    weights_changed = weights
    for i in range(len(weights_changed)):
        if "filter_{}".format(i) in below_tresh:
            weights_changed[i] = torch.zeros(weights_changed[i].size())
            
    # substitute of changed filters to the network     
    for name, param in model.named_parameters():
        if name == "conv_layers.0.weight":
            param = weights_changed
            
    # save the changed network
    torch.save(model.state_dict(), path+network_name+"_changed.model")



def read_sequences(file, chrom_list, filetype = "fasta"):
    """
    Function used to load sequences from a file

    :param file: file with sequences from the test set
    :param chrom_list: chromosomes from which the sequences for the test set are to come
    :param filetype: file type, default is fasta
    :return: test set sequence list
    """
    
    sequences = []
    
    for seq_record in SeqIO.parse(file, filetype):
        # only chromosomes from the test set
        if chrom_list:
            if seq_record.id in chrom_list:
                sequences.append(seq_record)
        else:
            sequences.append(seq_record)
    
    return sequences

class OHEncoder:

    def __init__(self, categories=np.array(['A', 'C', 'G', 'T'])):
        self.encoder = Encoder(sparse=False, categories=[categories])
        self.dictionary = categories
        self.encoder.fit(categories.reshape(-1, 1))

    def __call__(self, seq, info=False):
        seq = list(seq)
        if 'N' in seq:
            pos = [i for i, el in enumerate(seq) if el == 'N']
            if len(pos) <= 0.05*len(seq):
                if info:
                    print('{} unknown position(s) in given sequence - changed to random one(s)'.format(len(pos)))
                for p in pos:
                    seq[p] = random.choice(self.dictionary)
            else:
                return None
        if info:
            return True
        else:
            s = np.array(seq).reshape(-1, 1)
            return self.encoder.transform(s).T

    def decode(self, array):
        return ''.join([el[0] for el in self.encoder.inverse_transform(array.T)])


def create_dataset(sequences_path, chrom_list):
    """
    Function used to create a set of sequences to be used
    for testing original and modified models. From the sequence path, 
    it selects all FASTA files and creates a set of all test sequences.
    Sequences belonging to the test set are saved to a file.
    
    :param sequences_path: the path containing the sequences that are to be included in the test set
    :param chrom_list: chromosomes from which the sequences for the test set are to come
    :return: shuffled list of sequences from every collection data
    """

    np.random.seed(123)
    random.seed(123)

    def is_fasta(filename):
        """
        Function used to check whether a file is in FASTA format
        
        :param filename: file to be checked
        :return: boolean value whether the file is in FASTA format
        """
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    
    # all files in the path with sequences
    files = [f for f in os.listdir(sequences_path) if os.path.isfile(os.path.join(sequences_path, f)) and not f.startswith(".")]
    # only FASTA files
    fasta_files = [f for f in files if is_fasta(sequences_path+f) == True]
    
    dictionary = {}
    lengths = []

    for fasta_file in fasta_files:
        try:
            sequences = read_sequences(sequences_path+fasta_file, chrom_list)
        except:
            pass
        dictionary[Path(fasta_file).stem] = sequences
        lengths.append(len(sequences))

    output = []

    for key in dictionary:
        output.extend(dictionary[key])

    # shuffle all of the test sequences
    shuffled = sorted(output, key=lambda k: random.random())
    
    path = './prediction/'
    if not os.path.exists(path):
        os.makedirs(path)

    # save test set sequences to a file
    with open(path+"testset.txt", "w") as file:
        for seq in shuffled:
            SeqIO.write(seq, file, 'fasta')
        
    return shuffled


class SeqsDataset(Dataset):

    def __init__(self, sequences_path, chrom_list, seq_len=2000):

        np.random.seed(123)

        self.classes = ['promoter active', 'nonpromoter active', 'promoter inactive', 'nonpromoter inactive']
        self.num_classes = len(self.classes)
        self.seq_len = seq_len
        self.encoder = OHEncoder()
        self.data = create_dataset(sequences_path, chrom_list)
        for i in range(len(self.data)):
            seq, label = self.__getitem__(i, info=False)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index, info=False):

        seq = self.data[int(index)].seq.upper()
        id = self.data[int(index)].description.split(" ")
        label = self.classes.index(id[len(id)-2]+" "+id[len(id)-1])
        
        encoded_seq = self.encoder(seq, info=info)
        #encoded_seq = one_hot_encode(seq).transpose()
        X = torch.tensor(encoded_seq)
        X = X.reshape(1, *X.size())
        y = torch.tensor(label)
        return X, y


def test_network(modelfile, sequences_path, chrom_list):
    """
    Function used for testing the original and modified models.
    Creates prediction results and saves them to a file. Checks
    fluctuations between predictions. Performs prediction evaluation using
    AUC and other metrics. Compares predictions for the changed and original models
    
    :param modelfile: model file (network)
    :param sequences_path: path to the sequences used to train and test the network
    :param chrom_list: chromosomes from which the sequences for the test set are to come
    """
    
    networkname = Path(modelfile).stem

    batch_size = 64

    classes = ['promoter active', 'nonpromoter active', 'promoter inactive', 'nonpromoter inactive']

    num_classes = len(classes)
    
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    model = CustomNetwork()
    model.load_state_dict(torch.load(modelfile, map_location=torch.device(device)))

    dataset = SeqsDataset(sequences_path, chrom_list)
    testloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size)
    
    prediction = []
    true, scores = [], []
    pred_vals = []
    with torch.no_grad():
        model.eval()
        confusion_matrix = np.zeros((num_classes, num_classes))
        loss_neurons = [[] for _ in range(num_classes)]
        # saving prediction fluctuations within the model
        path = "./prediction/diff_prediction/"
        if not os.path.exists(path):
            os.makedirs(path)
        with open(path+"{}_pred_diff.csv".format(networkname), "w") as pred_diff_file:
            writer = csv.writer(pred_diff_file)
            writer.writerow(["Tensor", "Seq", "Class", "Prediction", "Pred value {}".format(classes[0]), "Pred value {}".format(classes[1]), "Pred value {}".format(classes[2]), "Pred value {}".format(classes[3])])
            for i, (seqs, labels) in enumerate(testloader):
                if use_cuda:
                    seqs = seqs.cuda()
                    labels = labels.cuda()
                seqs = seqs.float()
                labels = labels.long()
                outputs = model(seqs)
                pred_vals.append(outputs)
                for o, l in zip(outputs, labels):
                    loss_neurons[l].append(-math.log((math.exp(o[l])) / (sum([math.exp(el) for el in o]))))
                _, predicted = torch.max(outputs, 1)
                for ind, label, outp in zip(predicted, labels.cpu(), outputs):
                    confusion_matrix[ind][label] += 1
                prediction.append(predicted)
                for j in range(len(labels.tolist())):
                    if labels.tolist()[j] != predicted.tolist()[j]:
                        out_vals = outputs.tolist()[j]
                        writer.writerow(["tensor_{}".format(i), "seq_{}".format(i*batch_size+j), classes[labels.tolist()[j]], classes[predicted.tolist()[j]], out_vals[0], out_vals[1], out_vals[2], out_vals[3]])
                true += labels.tolist()
                scores += outputs.tolist()
                
   # saving metrics
    path = "./prediction/metrics/"
    if not os.path.exists(path):
        os.makedirs(path)
    with open(path+"{}_metrics.csv".format(networkname), "w") as file_metrics:
        losses, sens, spec = calculate_metrics(confusion_matrix, loss_neurons)
        auc = calculate_auc(true, scores)
        writer = csv.writer(file_metrics)
        writer.writerow(["AUC", "Sensitivity", "Specificity", "Losses"])
        writer.writerow([auc, sens, spec, losses])
    
    # saving prediction values
    path = "./prediction/outputs/"
    if not os.path.exists(path):
        os.makedirs(path)
    with open(path+"{}_outputs.txt".format(networkname), "w") as file_outputs:
        for out in pred_vals:
            file_outputs.write(str(out.tolist()))
            file_outputs.write("\n")

    # saving output predictions
    path = "./prediction/"
    if not os.path.exists(path):
        os.makedirs(path)
    save_to_file(prediction, path+"{}_pred.txt".format(networkname))


def calculate_metrics(confusion_matrix, losses):
    """
    Function used to calculate metrics
    
    :param confusion_matrix: confusion matrix
    :param losses: losses
    :return: losses, sensitivity, specificity values
    """
    num_classes = confusion_matrix.shape[0]
    sens, spec = [], []
    for cl in range(num_classes):
        tp = confusion_matrix[cl][cl]
        fn = sum([confusion_matrix[row][cl] for row in range(num_classes) if row != cl])
        tn = sum([confusion_matrix[row][col] for row, col in product(range(num_classes), repeat=2)
                  if row != cl and col != cl])
        fp = sum([confusion_matrix[cl][col] for col in range(num_classes) if col != cl])
        sens += [float(tp) / (tp + fn) if (tp + fn) > 0 else 0.0]
        spec += [float(tn) / (tn + fp) if (tn + fp) > 0 else 0.0]
    loss = [mean(el) if el else None for el in losses]
    return loss, sens, spec

def calculate_auc(true, scores):
    """
    Function used to calculate AUC
    
    :param true: true classes
    :param scores: predicted classes
    :return: AUC value
    """
    num_classes = len(scores[0])
    auc = [0 for _ in range(num_classes)]
    for neuron in range(num_classes):
        y_true = [1 if el == neuron else 0 for el in true]
        y_score = [el[neuron] for el in scores]
        auc[neuron] = roc_auc_score(y_true, y_score)
    return auc


def save_to_file(prediction, filename):
    """
    Function used to save the output predictions (assigned classes) to a file
    
    :param prediction: the output predictions
    :param filename: the output file name
    """

    classes = ['promoter active', 'nonpromoter active', 'promoter inactive', 'nonpromoter inactive']
    
    with open(filename, "w") as file:
        for pred in prediction:
            pred_list = pred.tolist()
            class_list = [classes[i] for i in pred_list]
            file.write(str(class_list))
            file.write("\n")

def main():
    
    parser = argparse.ArgumentParser(description='Perform a network test and save the output predictions, prediction fluctuations, and metrics to files')
    parser.add_argument('--model', action='store', dest='model', help='The path where the network file is located', required=True)
    parser.add_argument('--sequences', action='store', dest='seqs', help='The path to the sequences that will be used for the test set. The sequences need to be in FASTA format', required=True)
    parser.add_argument('--chrom_list', nargs='*', dest='chrom', default=["chr21", "chr22", "chrX", "chrY"], help='Chromosomes from which the sequences for the test set are to come; default=["chr21", "chr22", "chrX", "chrY"]')
    args = parser.parse_args()

    if args.model == None: 
        print("No model file provided (see help; -h)")
        sys.exit(-1)
    else:
        if not args.model.endswith('.model'):
            print("Invalid model format")
            sys.exit(-1)

    if args.seqs == None: 
        print("No sequences path provided (see help; -h)")
        sys.exit(-1)
    
    if not args.chrom == None:
        if not all(isinstance(item, str) for item in args.chrom):
            print("All elements (chromosome names) must be string instances")
            sys.exit(-1)
    
    print("Arguments used: "+', '.join(f'{k}={v}' for k, v in vars(args).items()))
    
    # setting the filter weights for the network to 0
    print("Setting the filter weights to 0...")
    change_weights(args.model)
    print("Done.")
    
    # testing the original and changed model on the same data
    print("Testing the original model...")
    test_network(args.model, args.seqs, args.chrom)
    print("Done.")
    print("Testing the changed model...")
    test_network("./changed_models/{}_changed.model".format(Path(args.model).stem), args.seqs, args.chrom)
    print("Done.")
    
if __name__ == "__main__":
    main()