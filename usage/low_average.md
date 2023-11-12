### Genome wide prediction of regulatory elements using neural networks
###### Interpretation of filters from convolutional layers of neural networks based on DNA motifs
--------------------------------------------------------------------------------
## Analysis of filters with low average weight value

<a name="test_changed_models.py"/>

<a href=“scripts/test_changed_models.py”><h4>test_changed_models.py</h4></a>

Extract filters from the network, calculate their stats and then based on that divide them into those with high and low average weight values.

| Arguments | Type | Description |
| --- | --- | --- |
| --model | string | The path where the network file is located |
| --sequences | string | The path to the sequences (FASTA) that will be used for the test set. |
| --chrom_list | string | Chromosomes (zero or more) from which the sequences for the test set are to come [Default: ["chr21", "chr22", "chrX", "chrY"]] |

