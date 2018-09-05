[![DOI](https://zenodo.org/badge/106058128.svg)](https://zenodo.org/badge/latestdoi/106058128) Latest public release

[![Build Status](https://travis-ci.org/HegemanLab/w4mcorcov_galaxy_wrapper.svg?branch=master)](https://travis-ci.org/HegemanLab/w4mcorcov_galaxy_wrapper) Current build status for master branch on GitHub

[Repository 'w4mcorcov' in Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository?repository_id=96046af0e175c57d)


# OPLS-DA Contrasts - a Galaxy tool

The OPLS-DA Contrasts Galaxy tool presents an emulation of the SIMCA® OPLS-DA® S-PLOT®.

For each pair of treatment levels, this tool plots:
- an emulation of the S-PLOT®, i.e., for each feature, the correlation of the treatment level with the "parallel component" of the OPLS-DA® projection of multivariate data against the covariance of the same.  In other words, this is a plot of the consistency of the effect of each feature on *the variation **between** the two treatment levels* against the degree to which the feature influences that variation.
- two OPLS-DA® plots (scores and cross-validation) from [the Bioconductor `ropls` package, doi:10.18129/B9.bioc.ropls](https://dx.doi.org/10.18129/B9.bioc.ropls).    
- an analog to the S-PLOT® for orthogonal features, i.e., for each feature, the correlation of the treatment level with the "orthogonal component" of the OPLS-DA® projection of multivariate data against the covariance of the same.  In other words, this is a plot of the consistency of the effect of each feature on the *variation **within** the two treatment levels* against the degree to which the feature influences that variation.

Optionally, the features plotted may be limited to those found to be significant in a univariate test.

The original OPLS-DA® S-PLOT® which this tool emulates is described in:
> Wiklund *et al.* (2008). **Visualization of GC/TOF-MS-Based Metabolomics Data for Identification of Biochemically Interesting Compounds Using OPLS Class Models.** In *Analytical Chemistry, 80 (1), pp. 115–122.* [doi:10.1021/ac0713510](https://dx.doi.org/10.1021/ac0713510)

# Notes

1. OPLS-DA®, SIMCA®, and S-PLOT® are registered trademarks of the Umetrics company, [http://umetrics.com/about-us/trademarks](http://umetrics.com/about-us/trademarks).
2. Available in the Galaxy Tool Shed at [https://toolshed.g2.bx.psu.edu/repository?repository_id=96046af0e175c57d](https://toolshed.g2.bx.psu.edu/repository?repository_id=96046af0e175c57d)

# Release notes

0.98.15

- Fix `unselected levels included in other` bug [https://github.com/HegemanLab/w4mclassfilter/issues/8](https://github.com/HegemanLab/w4mclassfilter/issues/8).

0.98.14

- Fix `read_data_frame` bug inherited by copy-and-paste from [https://github.com/HegemanLab/w4mclassfilter/issues/1](https://github.com/HegemanLab/w4mclassfilter/issues/1).
- Added missing support for underscores in class levels, fixing issue [https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/issues/7](https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/issues/7).

0.98.13

- Travis CI support, rearrange directories, 

0.98.12

- merge: iuc standards pull requests [3](https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/pull/3) and [4](https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/pull/4).

0.98.11

- bug fix: Readdress issue 2 - features now are pareto-scaled *and centered* per Wikland *op cit.*.

0.98.10

- new feature: C-plots of VIP versus correlation or relative covariance.
- bug fix: Handle issue 2 - features now are only pareto-scaled per Wikland *op cit.*.

0.98.9

- bug fix: Handle issue 1 - handle features removed by ropls because variance is less than 2.2e-16.

0.98.8

- new feature: Replace loadings plot with correlation-versus-covariance plot for orthogonal features, i.e., the consistency of features influencing within-treatment variation (which is linearly related to the loading of the orthogonal projection) versus consistency.  This eliminates the need for the parameter to suppress labels for features with extreme orthogonal loadings

0.98.7

- bug fix: Handle case of a treatment level with only one sample.

0.98.6

- bug fix: Set 'crossvalI' param (of R function 'ropls::opls') to the number of samples when the there are fewer than seven samples.

0.98.5

- bug fix: Fit feature-labels within clipping region of cor-vs.cov plot
- new feature: optionally (and by default) suppress labels for features with extreme orthogonal loadings

0.98.3

- Add support for two-level factors
- Add adjusted mz and rt to output tables
- Allow explicitly setting the number of features with extreme loadings to be labelled on the correlation vs. covariance plot
- Add loadings to corcov table

0.98.2

- first release

