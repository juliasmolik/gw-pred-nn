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


<h3>Analysis of structural features of DNA sequences</h3>

<a name="prepare_dnashape_data.py"/>

<a href=“scripts/prepare_dnashape_data.py”><h4>prepare_dnashape_data.py</h4></a>

Create all possible sequences of length k. Calculate the convolution of each of them with the weight values in the analyzed filter k-mer.
| Arguments | Type | Description |
| --- | --- | --- |
| --network | string | The network name where the filters came from |
| --k | int | Kmer value (length of the divided filter matrix as well as sequences) [Default: 5] |


<a name="dnashape_analysis.Rmd"/>

<a href=“dnashape_analysis.Rmd”><h4>dnashape_analysis.Rmd</h4></a>

Calculate DNAshape parameters for each of the possible sequences of previously selected length
| Variables | Type | Description |
| --- | --- | --- |
| path_to_sequences | string | The path where the sequences of previously selected length are located |
| output_dir | string | The path to write DNAshape parameters for each sequence |


<a name="calculate_correlation.py"/>

<a href=“scripts/calculate_correlation.py”><h4>calculate_correlation.py</h4></a>

Calculate the correlation of the convolution value depending on the analyzed sequence with the DNAshape parameters of a given sequence.
| Arguments | Type | Description |
| --- | --- | --- |
| --network | string | The network name where the filters came from |

<a name="significant_qvalues.py"/>

<a href=“scripts/significant_qvalues.py”><h4>significant_qvalues.py</h4></a>

Calculate q-values and look for filter kmers significantly correlated with each parameter.
| Arguments | Type | Description |
| --- | --- | --- |
| --network | string | The network name where the filters came from |
| --thresh | value that can be converted to float | Threshold for q-value cut-off [Default: 1e-15] |
