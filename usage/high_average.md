### Genome wide prediction of regulatory elements using neural networks
###### Interpretation of filters from convolutional layers of neural networks based on DNA motifs
--------------------------------------------------------------------------------
## Analysis of filters with high average weight value

<h3>Creating PFM matrices</h3>

<a name="scan_sequences.py"/>

<a href=“scripts/scan_sequences.py”><h4>scan_sequences.py</h4></a>

Scan training dataset sequences using chosen filters and create PFM matrices based on n best sub-sequences for each filter.

| Arguments | Type | Description |
| --- | --- | --- |
| --sequences | string | The path where the training dataset sequences are located |
| --filters | string | The path where the files with filters weights are located |
| --dataset | string | The dataset name from which the training sequences came |
| --network | string | The network name where the filters came from |
| --n | int | Number of sub-sequences used to create PFM matrices [Default: 100] |


<h3>Scanning the transcription factor motif database</h3>

<a name="pfm_to_kmers.py"/>

<a href=“scripts/pfm_to_kmers.py”><h4>pfm_to_kmers.py</h4></a>

Divide filters into kmers and for each of them create a PFM matrix. Calculate the IC value and save the results to a file.

| Arguments | Type | Description |
| --- | --- | --- |
| --pfm | string | The path where the PFM matrices are located |
| --dataset | string | The dataset name from which the training sequences came |
| --network | string | The network name where the filters came from |
| --epsilon | float |The pseudo counts value [Default: 1] |
| --k | int | Kmer value (length of the divided filter matrix) [Default: 7] |
| --thresh | float | IC value cutoff threshold [Default: 6.5] |


<a name="tomtom_analysis.py"/>

<a href=“scripts/tomtom_analysis.py”><h4>tomtom_analysis.py</h4></a>

Perform Tomtom analysis. Create statistics files summarizing Tomtom results. 

| Arguments | Type | Description |
| --- | --- | --- |
| --ppm | string | The path where the filter kmer PPM matrices are located |
| --dataset | string | The dataset name from which the training sequences came |
| --network | string | The network name where the filters came from |