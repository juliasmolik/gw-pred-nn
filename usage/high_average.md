### Genome wide prediction of regulatory elements using neural networks
###### Interpretation of filters from convolutional layers of neural networks based on DNA motifs
--------------------------------------------------------------------------------
## Analysis of filters with high average weight value

<h3>Creating PFM matrices</h3>

<a name="scan_sequences.py"/>

<a href=“scripts/scan_sequences.py”><h4>scan_sequences.py</h4></a>

Description

| Arguments | Type | Description |
| --- | --- | --- |
| --sequences | string | The path where the training dataset sequences are located |
| --filters | string | The path where the files with filters weights are located |
| --dataset | string | The dataset name from which the training sequences came |
| --network | string | The network name where the filters came from |
| --n | ing | Number of sub-sequences used to create PFM matrices [Default: 100] |
