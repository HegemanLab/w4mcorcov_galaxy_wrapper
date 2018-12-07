w4msalience <- function(
  data_matrix   # a matrix of intensities; features as rows, and samples as columns
, sample_class  # a vector of sample class-levels; length(sample_class) == ncol(data_matrix)
, failure_action = stop
) {
  library(stats)
  # begin sanity checks
  if ( !is.vector(sample_class) || !( is.character(sample_class) || is.factor(sample_class) ) ) {
    failure_action("w4msalience:  Expected sample_class to be a vector of characters of factor-levels")
    return (NULL)
  }
  if ( !is.matrix(data_matrix) && !is.data.frame(data_matrix) ) {
    failure_action("w4msalience:  Expected data_matrix to be a matrix (or data.frame) of numeric")
    return (NULL)
  }
  # transpose data_matrix so that columns are the features
  t_data_matrix <- t(data_matrix)
  if ( !is.matrix(t_data_matrix) || !is.numeric(t_data_matrix) ) {
    failure_action("w4msalience:  Expected data_matrix to be a matrix (or data.frame) of numeric")
    return (NULL)
  }
  n_features <- ncol(t_data_matrix)
  n_features_plus_1 <- 1 + n_features
  features   <- colnames(t_data_matrix)
  n_samples  <- nrow(t_data_matrix)
  if ( length(sample_class) != n_samples ) {
    strF(data_matrix)
    strF(sample_class)
    failure_action(
      sprintf(
        "w4msalience:  The data_matrix has %d samples but sample_class has %d"
      , n_samples
      , length(sample_class)
      )
    )
    return (NULL)
  }
  # end sanity checks

  # "For each feature, 'select sample_class, median(intensity) from feature group by sample_class'."
  # The first column(s) of the result of aggregate has the classifier value(s) specified in the 'by' list.
  medianOfFeatureBySampleClassLevel <- aggregate(
      x = as.data.frame(t_data_matrix)
    , by = list(sample_class)
    , FUN = "median"
    )

  # "For each feature, 'select sample_class, max(intensity) from feature group by sample_class'."
  maxOfFeatureBySampleClassLevel <- aggregate(
      x = as.data.frame(t_data_matrix)
    , by = list(sample_class)
    , FUN = "max"
    )

  # "For each feature, 'select sample_class, rcv(intensity) from feature group by sample_class'."
  #   cv is less robust; deviation from normality degrades performance
  #     cv(x) == sd(x) / mean(x)
  #   rcv is a "robust" coefficient of variation, expressed as a proportion
  #     rcv(x) == mad(x) / median(x)
  madOfFeatureBySampleClassLevel <- aggregate(
      x = as.data.frame(t_data_matrix)
    , by = list(sample_class)
    , FUN = "mad"
  )
  # "robust coefficient of variation", i.e.,
  #    mad(feature-intensity for class-level max_level) / median(feature-intensity for class-level max_level)
  rcvOfFeatureBySampleClassLevel <- as.matrix(
    madOfFeatureBySampleClassLevel[ , 2:n_features_plus_1] /
      medianOfFeatureBySampleClassLevel[,2:n_features_plus_1]
  )
  rcvOfFeatureBySampleClassLevel[is.nan(rcvOfFeatureBySampleClassLevel)] <-
    max(9999,max(rcvOfFeatureBySampleClassLevel, na.rm = TRUE))

  relative_salient_distance <- sapply(
    X = 1:n_features
  , FUN = function(x) {
      my_df <- data.frame(
        max = maxOfFeatureBySampleClassLevel[ , 1 + x, drop = TRUE]
      , median = medianOfFeatureBySampleClassLevel[ , 1 + x, drop = TRUE]
      , mad = madOfFeatureBySampleClassLevel[ , 1 + x, drop = TRUE]
      )
      my_df <- my_df[ order(my_df$max, decreasing = TRUE), ][  1:2, ]
      result <-
        ( my_df$median[1] - my_df$median[2] ) /
        sqrt( my_df$mad[1] * my_df$mad[2] )
      if (is.infinite(result) || is.nan(result))
        result <- 0
      return(result)
    }
  )


  # "For each feature, 'select max(max_feature_intensity) from feature'."
  maxApplyMaxOfFeatureBySampleClassLevel <- sapply(
      X = 1:n_features
    , FUN = function(i) {
        match(
          max(maxOfFeatureBySampleClassLevel[, i + 1])
        , maxOfFeatureBySampleClassLevel[, i + 1]
        )
      }
    )

  # "For each feature, 'select mean(median_feature_intensity) from feature'."
  meanApplyMedianOfFeatureBySampleClassLevel <- sapply(
      X = 1:n_features
    , FUN = function(i) mean(medianOfFeatureBySampleClassLevel[, i + 1])
    )

  # Compute the 'salience' for each feature, i.e., how salient the intensity of a feature
  #   is for one particular class-level relative to the intensity across all class-levels.
  salience_df <- data.frame(
    # the feature name
    feature = features
    # the name (or factor-level) of the class-level with the highest median intensity for the feature
  , max_level = medianOfFeatureBySampleClassLevel[maxApplyMaxOfFeatureBySampleClassLevel,1]
    # the median intensity for the feature and the level max_level
  , max_median = sapply(
        X = 1:n_features
      , FUN = function(i) {
          maxOfFeatureBySampleClassLevel[
            maxApplyMaxOfFeatureBySampleClassLevel[i]
          , 1 + i
          ]
        }
    )
    # the distance between the maximum intensities for the feature at the two highest levels
  , relative_salient_distance = relative_salient_distance
    # the coefficient of variation (expressed as a proportion) for the intensity for the feature and the level max_level
  , max_rcv = sapply(
        X = 1:n_features
      , FUN = function(i) {
          rcvOfFeatureBySampleClassLevel[
            maxApplyMaxOfFeatureBySampleClassLevel[i]
          , i
          ]
        }
    )
    # the mean of the medians of intensity for all class-levels for the feature
  , mean_median = meanApplyMedianOfFeatureBySampleClassLevel
    # don't coerce strings to factors (this is a parameter for the data.frame constructor, not a column of the data.frame)
  , stringsAsFactors = FALSE
  )
  # raw salience is the ratio of the most-prominent level to the mean of all levels for the feature
  salience_df$salience <- sapply(
      X = 1:nrow(salience_df)
    , FUN = function(i) {
        with(
          salience_df[i,]
        , if (mean_median > 0) max_median / mean_median else 0
        )
      }
    )
  # "robust coefficient of variation"
  salience_df$salient_rcv <- salience_df$max_rcv

  return (salience_df)
}

