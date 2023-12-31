### Genome wide prediction of regulatory elements using neural networks
###### Interpretation of filters from convolutional layers of neural networks based on DNA motifs
--------------------------------------------------------------------------------
## Extracting and grouping filters

<a name="extract_filters.py"/>

<a href=“scripts/extract_filters.py”><h4>extract_filters.py</h4></a>

Extract filters from the network, calculate their stats and then based on that divide them into those with high and low average weight values.

| Arguments | Type | Description |
| --- | --- | --- |
| --model | string | The path where the network file is located |
| --thresh | value that can be converted to float | The threshold above which the average filter weight is considered high [Default: 1e-5] |
