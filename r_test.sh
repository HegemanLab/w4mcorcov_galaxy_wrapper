#!/bin/bash
__tool_directory__=.
# sampleMetadata	k10
sampleMetadata_in=test-data/input_sampleMetadata.tsv
# variableMetadata	k10_kruskal_k2.k1_sig	k10_kruskal_k3.k1_sig	k10_kruskal_k4.k1_sig	k10_kruskal_k3.k2_sig	k10_kruskal_k4.k2_sig	k10_kruskal_k4.k3_sig
variableMetadata_in=test-data/input_variableMetadata.tsv
# rows labeled as first column of variableMetadata; columns, as first column of sampleMetadata
dataMatrix_in=test-data/input_dataMatrix.tsv
# Name of the statistical test that was run in Univariate to produce the variable metadata file - it is a portion of the column names in that file
tesC=kruskal
# Name of the column of the sample metadata table corresponding to the qualitative variable used to define the contrasts - it is a portion of the column names in the variable metadata file
facC=k10
# Retain only pairwise-significant features" help="When true, analyze only features that differ significantly for the pair of levels being contrasted; when false, include any feature that varies significantly across all levels.
pairSigFeatOnly=TRUE
labelFeatures=TRUE
# comma-separated level-names (or comma-less regular expressions to match level-names) to consider in analysis; must match at least two levels; may include wild cards or regular expressions
levCSV=k[12],k[3-4]
# how to specify levels generically
matchingC=regex
#   pdf1: contrast grid
#     * with S-PLOTs in the upper triangle
#     * PLS score-plots in the lower triangle
#     * with level-labels along the diagonal
contrast_detail=test-data/output_contrast_detail.pdf
#   tsv1: cor and cov dataframe with colums:
#     * feature-ID
#     * factor-level 1
#     * factor-level 2, lexically greater than level 1
#     * Wiklund_2008 correlation
#     * Wiklund_2008 covariance
#     * Galindo_Prieto_2014 VIP for predictive components, VIP[4,p]
#     * Galindo_Prieto_2014 VIP for orthogonal components, VIP[4,o]
#     * Salient level, i.e., for the feature, the class-level having the greatest median intensity
#     * Salient robust coefficient of variation, i.e., for the feature, the mean absolute deviation of the intensity for the salient level divided by the median intensity for the salient level
#     * Raw salience, i.e., for the feature, the median of the class-level having the greatest intensity divided by the mean of the medians for all class-levels.
#     * (When filtering on significance of univariate tests) Significance of test of null hypothesis that there is no difference between the two classes, i.e, the pair-wise test.


contrast_corcov=test-data/output_contrast_corcov.tsv

# Run the script
bash -c " cd $__tool_directory__; \
  Rscript w4mcorcov_wrapper.R \
  dataMatrix_in '$dataMatrix_in' \
  sampleMetadata_in '$sampleMetadata_in' \
  variableMetadata_in '$variableMetadata_in' \
  tesC '$tesC' \
  facC '$facC' \
  pairSigFeatOnly '$pairSigFeatOnly' \
  levCSV '$levCSV' \
  matchingC '$matchingC' \
  contrast_detail '$contrast_detail' \
  contrast_corcov '$contrast_corcov' \
  labelFeatures '$labelFeatures' \
  "

# Repeat the test with pairSigFeatOnly FALSE
pairSigFeatOnly=FALSE
contrast_detail=test-data/output_contrast_detail_all.pdf
contrast_corcov=test-data/output_contrast_corcov_all.tsv

# Run the script
bash -c " cd $__tool_directory__; \
  Rscript w4mcorcov_wrapper.R \
  dataMatrix_in '$dataMatrix_in' \
  sampleMetadata_in '$sampleMetadata_in' \
  variableMetadata_in '$variableMetadata_in' \
  tesC '$tesC' \
  facC '$facC' \
  pairSigFeatOnly '$pairSigFeatOnly' \
  levCSV '$levCSV' \
  matchingC '$matchingC' \
  contrast_detail '$contrast_detail' \
  contrast_corcov '$contrast_corcov' \
  labelFeatures '$labelFeatures' \
  "

# Repeat the test with test none
tesC=none
labelFeatures=FALSE
contrast_detail=test-data/output_contrast_detail_global.pdf
contrast_corcov=test-data/output_contrast_corcov_global.tsv

# Run the script
bash -c " cd $__tool_directory__; \
  Rscript w4mcorcov_wrapper.R \
  dataMatrix_in '$dataMatrix_in' \
  sampleMetadata_in '$sampleMetadata_in' \
  variableMetadata_in '$variableMetadata_in' \
  tesC '$tesC' \
  facC '$facC' \
  pairSigFeatOnly '$pairSigFeatOnly' \
  levCSV '$levCSV' \
  matchingC '$matchingC' \
  contrast_detail '$contrast_detail' \
  contrast_corcov '$contrast_corcov' \
  labelFeatures '$labelFeatures' \
  "

