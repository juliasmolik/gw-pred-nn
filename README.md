<img align="left" src="./logo.png" width="25%"> <h1 style="margin-left:300px;"> <br> <br> Genome wide prediction of regulatory <br> elements using neural networks

#### Interpretation of filters from convolutional layers of neural networks based on DNA motifs

In recent years, the advent of high-throughput sequencing technologies has generated vast amounts of genomic data, enabling the development of computational approaches to predict regulatory elements. The Basset tool ([Kelley et al., Genome Research 2016](https://genome.cshlp.org/content/26/7/990)), published in the last decade, uses convolutional neural networks to learn the functional activity of DNA sequences. It predicts cell type-specific chromatin openness (DNA sequence accessibility) in 164 cell types. The architecture of the neural network from the Basset tool has been modified and adapted to the problem of identifying the regulatory activity of regions of the genome in one type of tissue (human brain tumor tissue), as well as assigning functions to these regions (promoter or enhancer). The changed models are the result of [Marlena Osipowicz's master's thesis](https://github.com/marnifora/magisterka).

This repository is the result of Julia Smolik's master's thesis carried out under the supervision of Dr. Magdalena Machnicka at the Faculty of Mathematics, Informatics, and Mechanics of the University of Warsaw in Poland. The aim of the work was:
1. Detailed analysis of filters included in convolutional neural networks predicting the activity of regulatory areas, as well as identification of the similarity of the filters to DNA sequence motifs recognized by transcription factors.
2. Analysis of the connections of motifs encoded in the filters with DNA structural features described by a set of parameters ([Zhou et al., NAR 2013](https://academic.oup.com/nar/article/41/W1/W56/1105326)).

The results of all the analyzes are part of a publication on genome-wide prediction of regulatory elements using neural networks (manuscript in preparation).

---------------------------------------------------------------------------------------------------
### Documentation

---------------------------------------------------------------------------------------------------
### License
[MIT](https://choosealicense.com/licenses/mit/)

---------------------------------------------------------------------------------------------------
### Authors
[Julia Smolik](https://github.com/juliasmolik) <br>
[Marlena Osipowicz](https://github.com/marnifora) <br>
Magdalena Machnicka <br>
Bartek Wilczyński
