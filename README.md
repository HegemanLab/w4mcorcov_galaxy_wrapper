[![DOI](https://zenodo.org/badge/106058128.svg)](https://zenodo.org/badge/latestdoi/106058128)

[Repository 'w4mcorcov' in Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository?repository_id=96046af0e175c57d)


# w4mcorcov\_galaxy\_wrapper

This Galaxy tool presents an emulation of the SIMCA® OPLS-DA® S-PLOT® along with several plots from [the Bioconductor `ropls` package](https://dx.doi.org/10.18129/B9.bioc.ropls).
(OPLS-DA®, SIMCA®, and S-PLOT® are registered trademarks of the Umetrics company, [http://umetrics.com/about-us/trademarks](http://umetrics.com/about-us/trademarks).)
This tool plots the correlation of each feature with the OPLS-DA projection of multivariate data for two classes of data against the covariance of the same.

The original plot which this tool emulates is described in:
> Wiklund, Susanne and Johansson, Erik and Sjöström, Lina and Mellerowicz, Ewa J. and Edlund, Ulf and Shockcor, John P. and Gottfries, Johan and Moritz, Thomas and Trygg, Johan (2008). Visualization of GC/TOF-MS-Based Metabolomics Data for Identification of Biochemically Interesting Compounds Using OPLS Class Models. In Analytical Chemistry, 80 (1), pp. 115–122. doi:10.1021/ac0713510

# notes and thoughts

1. For lines from the corcov dataset for a single contrast (i.e., with a common combination of factorLevel1 and factorLevel2), there is a linear relationship between correlation and loadp; correlation has the advantage that it is more directly interpretable regarding how consistently a feature varies with the predictor.
