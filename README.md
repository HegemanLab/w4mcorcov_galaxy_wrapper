[![DOI](https://zenodo.org/badge/106058128.svg)](https://zenodo.org/badge/latestdoi/106058128)

[Repository 'w4mcorcov' in Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository?repository_id=96046af0e175c57d)


# OPLS-DA Contrasts - a Galaxy tool

The OPLS-DA Contrasts Galaxy tool presents an emulation of the SIMCA® OPLS-DA® S-PLOT®.

For each pair of treatment levels, this tool plots:
- an emulation of the S-PLOT®, i.e., for each feature, the correlation of the treatment level with the "parallel component" of the OPLS-DA® projection of multivariate data against the covariance of the same.  This plots the consistency of the effect of each feature on *the variation **between** the two treatment levels* against the degree to which the feature influences that variation.
- two OPLS-DA® plots (scores and cross-validation) from [the Bioconductor `ropls` package, doi:10.18129/B9.bioc.ropls](https://dx.doi.org/10.18129/B9.bioc.ropls).    
- an analog to the S-PLOT® for orthogonal features, i.e., for each feature, the correlation of the treatment level with the "orthogonal component" of the OPLS-DA® projection of multivariate data against the covariance of the same.  This plots the consistency of the effect of each feature on the *variation **within** the two treatment levels* against the degree to which the feature influences that variation.

Optionally, the features plotted may be limited to those found to be significant in a univariate test.

The original OPLS-DA® S-PLOT® which this tool emulates is described in:
> Wiklund *et al.* (2008). **Visualization of GC/TOF-MS-Based Metabolomics Data for Identification of Biochemically Interesting Compounds Using OPLS Class Models.** In *Analytical Chemistry, 80 (1), pp. 115–122.* [doi:10.1021/ac0713510](https://dx.doi.org/10.1021/ac0713510)

# Notes and thoughts

1. OPLS-DA®, SIMCA®, and S-PLOT® are registered trademarks of the Umetrics company, [http://umetrics.com/about-us/trademarks](http://umetrics.com/about-us/trademarks).
2. For lines from the corcov dataset for a single contrast (i.e., with a common combination of factorLevel1 and factorLevel2), there is a linear relationship between correlation and loadp; correlation has the advantage that it is more directly interpretable regarding how consistently a feature varies with the predictor.
