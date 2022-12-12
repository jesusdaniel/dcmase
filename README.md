# Community detection in multilayer degree corrected stochastic block models via DC-MASE

The degree-corrected multiple adjacency spectral embedding (DC-MASE) obtains a joint embedding from multilayer network data to estimate community memberships in the multilayer degree-corrected stochastic blockmodel (DC-SBM).

## Overview

Given a collection of *L* adjacency matrices representing graphs with $n$ aligned vertices (e.g. multilayer networks), DC-MASE obtains a joint embedding matrix with n rows, where each row corresponds to a vertex in the networks. After this, communities are obtained by applying a clustering procedure (such as K-means) to partition the latent positions into K groups in order to estimate the memberships in the multilayer degree-corrected stochastic block model. 

The joint embedding is calculated by performing a separate *adjacency spectral embedding*  (ASE) for each graph, which consists in computing the eigendecomposition of each adjacency matrix, with a possible eigenvalue scaling, followed by a row-normalization step (such as dividing the rows by its L2 norm), and then performs a joint singular value decomposition of the concatenated row-normalized ASEs. A pictorial representation is presented below.

![dc-mase embedding](https://github.com/jesusdaniel/dcmase/blob/main/img/DC-MASE2.png?raw=true)

# R Code

The R code in this repository implements DC-MASE and other methods to perform community detection in multilayer networks. To use this code, download all the content from the [R/](https://github.com/jesusdaniel/dcmase/tree/main/R) folder. Some examples and simulation experiments from the paper are contained in [this folder](https://github.com/jesusdaniel/dcmase/tree/main/Experiments).


# Data
The networks encode weighted edges representing the monthly number of flights between pairs of US airports. These data were obtained from the [T-100 Segment (US Carriers Only)](https://www.transtats.bts.gov/Fields.asp?gnoyr_VQ=GEE) database US Bureau of Transportation Statistics (BTS). A post-processed version of this dataset is included in this repository. The DC-MASE algorithm is illustrated using these data.

<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/USmap-K4-dcmase.png" height="400" />
<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/Bmatrices-K4-dcmase.png" height="400" />
<img src="https://github.com/jesusdaniel/dcmase/blob/main/Figures/degreecorrections-K4-dcmase.png" height="400" />


# References
Agterberg, J., Lubberts, Z., Arroyo, J., Joint Spectral Clustering in Multilayer Degree-Corrected Stochastic Blockmodels, [![arXiv shield](https://img.shields.io/badge/arXiv-2212.05053-red.svg?style=flat)](https://arxiv.org/abs/2212.05053)
