# center with 'colMeans()'
# ref: http://gastonsanchez.com/visually-enforced/how-to/2014/01/15/Center-data-in-R/
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

#### OPLS-DA
algoC <- "nipals"

cor_vs_cov <- function(matrix_x, ropls_x) {
  x_class <- class(ropls_x)
  if ( !( as.character(x_class) == "opls" ) ) { # || !( attr(class(x_class),"package") == "ropls" ) ) {
    stop( "cor_vs_cov: Expected ropls_x to be of class ropls::opls but instead it was of class ", as.character(x_class) )
  }
  result <- list()
  # suppLs$algoC - Character: algorithm used - "svd" for singular value decomposition; "nipals" for NIPALS
  if ( ropls_x@suppLs$algoC == "svd") {
    # Equations (1) and (2) from Wiklund 2008, doi:10.1021/ac0713510
    # scoreMN - Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
    score_matrix <- ropls_x@scoreMN
    score_matrix_transposed <- t(score_matrix)
    cov_divisor <- nrow(matrix_x) - 1
    result$covariance <- sapply(
      X = 1:ncol(matrix_x)
    , FUN = function(x) score_matrix_transposed %*% matrix_x[,x] / cov_divisor
    )
    score_sd <- sapply(
      X = 1:ncol(score_matrix)
      , FUN = function(x) sd(score_matrix[,x])
    )
    # xSdVn - Numerical vector: variable standard deviations of the 'x' matrix
    xSdVn <- ropls_x@xSdVn
    result$correlation <- sapply(
      X = 1:ncol(matrix_x)
    , FUN = function(x) {
        ( score_matrix_transposed / score_sd ) %*% ( matrix_x[,x] / (xSdVn[x] * cov_divisor) )
    }
    )
  } else {
    # Equations (1) and (2) from *Supplement to* Wiklund 2008, doi:10.1021/ac0713510
    mag <- function(one_dimensional) sqrt(sum(one_dimensional * one_dimensional))
    mag_xi <- sapply(X = 1:ncol(matrix_x), FUN = function(x) mag(matrix_x[,x]))
    score_matrix <- ropls_x@scoreMN
    score_matrix_transposed <- t(score_matrix)
    score_matrix_magnitude <- mag(score_matrix)
    result$covariance <- score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * score_matrix_magnitude )
    result$correlation <- score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * mag_xi )
  }
  return (result)
}

# S-PLOT and OPLS reference: Wiklund_2008 doi:10.1021/ac0713510
corcov_calc <- function(calc_env, failure_action = stop) {
  if ( ! is.environment(calc_env) ) {
    failure_action("corcov_calc: fatal error - 'calc_env' is not an environment")
    return ( FALSE )
  }

  # extract command-line parameters from the environment
  vrbl_metadata <- calc_env$vrbl_metadata
  smpl_metadata <- calc_env$smpl_metadata
  data_matrix <- calc_env$data_matrix
  pairSigFeatOnly <- calc_env$pairSigFeatOnly
  facC <- calc_env$facC
  tesC <- calc_env$tesC
  # extract the levels from the environment
  levCSV <- calc_env$levCSV
  # matchingC is one of { "none", "wildcard", "regex" }
  matchingC <- calc_env$matchingC

  # transform wildcard to regex
  if (matchingC == "wildcard") {
    levCSV <- gsub("[.]", "[.]", levCSV)
    levCSV <- utils::glob2rx(levCSV, trim.tail = FALSE)
  }
  # function to determine whether level is a member of levCSV
  isLevelSelected <- function(lvl) {
    matchFun <- if (matchingC == "regex") grepl else `==`
    return(
      0 < sum(sapply(X = strsplit(x = levCSV, split = ",", fixed = TRUE)[[1]], FUN = matchFun, lvl))
    )
  }

  # transpose matrix because ropls matrix is the transpose of XCMS matrix
  # Wiklund_2008 centers and pareto-scales data before OPLS-DA S-plot
  # center
  cdm <- center_colmeans(t(data_matrix))
  # pareto-scale
  my_scale <- sqrt(apply(cdm, 2, sd, na.rm=TRUE))
  scdm <- sweep(cdm, 2, my_scale, "/")

  # pattern to match column names like k10_kruskal_k4.k3_sig
  col_pattern <- sprintf('^%s_%s_(.*)[.](.*)_sig$', facC, tesC)
  # column name like k10_kruskal_sig
  intersample_sig_col  <- sprintf('%s_%s_sig', facC, tesC)
  # get the facC column from sampleMetadata, dropping to one dimension
  smpl_metadata_facC <- smpl_metadata[,facC]

  # allocate a slot in the environment for the contrast_list, each element of which will be a data.frame with this structure:
  #   - feature ID
  #   - value1
  #   - value2
  #   - Wiklund_2008 correlation
  #   - Wiklund_2008 covariance
  #   - Wiklund_2008 VIP
  calc_env$contrast_list <- list()
  # for each column name, extract the parts of the name matched by 'col_pattern', if any
  the_colnames <- colnames(vrbl_metadata)
  col_matches <- lapply(
    X = the_colnames,
    FUN = function(x) {
      regmatches( x, regexec(col_pattern, x) )[[1]]
    }
  )
  # process columns matching the pattern
  for (i in 1:length(col_matches)) {
    # for each potential match of the pattern
    col_match <- col_matches[[i]]
    if (length(col_match) > 0) {
      # it's an actual match; extract the pieces, e.g., k10_kruskal_k4.k3_sig
      vrbl_metadata_col <- col_match[1]               # ^^^^^^^^^^^^^^^^^^^^^  # Column name
      fctr_lvl_1 <- col_match[2]                      #             ^^         # Factor-level 1
      fctr_lvl_2 <- col_match[3]                      #                ^^      # Factor-level 2
      # only process this column if both factors are members of lvlCSV
      is_match <- isLevelSelected(fctr_lvl_1) && isLevelSelected(fctr_lvl_2)
      # TODO delete next line displaying currently-processed column
      cat(sprintf("%s %s %s %s\n", vrbl_metadata_col, fctr_lvl_1, fctr_lvl_2, is_match))
      # choose only samples with one of the two factors for this column
      chosen_samples <- smpl_metadata_facC %in% c(fctr_lvl_1, fctr_lvl_2)
      # transpose matrix because ropls matrix is the transpose of XCMS matrix
      # extract only the significantly-varying features and the chosen samples
      fully_significant   <- 1 == vrbl_metadata[,vrbl_metadata_col] * vrbl_metadata[,intersample_sig_col]
      overall_significant <- 1 == vrbl_metadata[,intersample_sig_col]
      col_selector <- if ( pairSigFeatOnly ) fully_significant else overall_significant 
      my_matrix <- scdm[ chosen_samples, col_selector, drop = FALSE ]
      # ropls::strF(my_matrix)
      # predictor has exactly two levels
      predictor <- smpl_metadata_facC[chosen_samples]
      file.pdfC <- sprintf("%s.pdf",vrbl_metadata_col)
      print(file.pdfC)
      if (is_match && ncol(my_matrix) > 1 && length(unique(predictor))> 1) {
        my_oplsda <- opls(
            my_matrix
          , predictor
          , algoC = algoC
          , predI = 1
          , orthoI = 1
          , printL = FALSE
          , plotL = FALSE
          , parDevNewL =  TRUE
          , file.pdfC = file.pdfC
          )
        my_oplsda_plot <- function(oplsda) {
          force(oplsda)
          return (
            function() {
              # typeVc <- c("correlation",      # 1
              #             "outlier",          # 2
              #             "overview",         # 3
              #             "permutation",      # 4
              #             "predict-train",    # 5
              #             "predict-test",     # 6
              #             "summary",          # 7 = c(2,3,4,9)
              #             "x-loading",        # 8
              #             "x-score",          # 9
              #             "x-variance",       # 10
              #             "xy-score",         # 11
              #             "xy-weight"         # 12
              #            )                    # [c(4,3,8,9)]      
              # for (i in c(4,3,8,9)) {
              # }
              plot(
                oplsda
              , typeVc = c("correlation",      # 1
                           "outlier",          # 2
                           "overview",         # 3
                           "permutation",      # 4
                           "predict-train",    # 5
                           "predict-test",     # 6
                           "summary",          # 7 = c(2,3,4,9)
                           "x-loading",        # 8
                           "x-score",          # 9
                           "x-variance",       # 10
                           "xy-score",         # 11
                           "xy-weight"         # 12
                          )[c(4,3,8,9)]               # [7]      
              , parCexN = 0.4
              , parDevNewL = TRUE
              , parLayL = FALSE
              , parEllipsesL = TRUE
              , file.pdfC = file.pdfC
              )
            }
          )
        }
        my_plot <- my_oplsda_plot(my_oplsda)
        my_plot()
        my_cor_vs_cov <- cor_vs_cov(
          matrix_x        = my_matrix
          , ropls_x       = my_oplsda
          )
        with(
          my_cor_vs_cov
        , {
            min_x <- min(covariance)
            max_x <- max(covariance)
            lim_x <- max(sapply(X=c(min_x, max_x), FUN=abs))
            plot(y = correlation, x = covariance, type="p", xlim=c(-lim_x, lim_x), ylim=c(-1,+1))
            low_x <- -0.7 * lim_x
            high_x <- 0.7 * lim_x
            text(x = low_x, y = -0.15, labels = fctr_lvl_1)
            text(x = high_x, y = 0.15, labels = fctr_lvl_2)
            #text(y = correlation, x = covariance, labels = names(covariance))
          }
        )
      } else {
        my_oplsda <- NULL
        cat("NO PLOT\n")
      }
    }
  }

  return ( TRUE )
}

# Calculate data for correlation-versus-covariance plot
#   Adapted from:
#     Wiklund_2008 doi:10.1021/ac0713510
#     Galindo_Prieto_2014 doi:10.1002/cem.2627
#     https://github.com/HegemanLab/extra_tools/blob/master/generic_PCA.R
cor_vs_cov <- function(matrix_x, ropls_x) {
  x_class <- class(ropls_x)
  if ( !( as.character(x_class) == "opls" ) ) { # || !( attr(class(x_class),"package") == "ropls" ) ) 
    stop( "cor_vs_cov: Expected ropls_x to be of class ropls::opls but instead it was of class ", as.character(x_class) )
  }
  result <- list()
  # suppLs$algoC - Character: algorithm used - "svd" for singular value decomposition; "nipals" for NIPALS
  if ( ropls_x@suppLs$algoC == "nipals") {
    # Equations (1) and (2) from *Supplement to* Wiklund 2008, doi:10.1021/ac0713510
    mag <- function(one_dimensional) sqrt(sum(one_dimensional * one_dimensional))
    mag_xi <- sapply(X = 1:ncol(matrix_x), FUN = function(x) mag(matrix_x[,x]))
    score_matrix <- ropls_x@scoreMN
    score_matrix_transposed <- t(score_matrix)
    score_matrix_magnitude <- mag(score_matrix)
    result$covariance <- score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * score_matrix_magnitude )
    result$correlation <- score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * mag_xi )
  } else {
    # WARNING - untested code - I don't have test data to exercise this branch
    # Equations (1) and (2) from Wiklund 2008, doi:10.1021/ac0713510
    # scoreMN - Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
    score_matrix <- ropls_x@scoreMN
    score_matrix_transposed <- t(score_matrix)
    cov_divisor <- nrow(matrix_x) - 1
    result$covariance <- sapply(
      X = 1:ncol(matrix_x)
    , FUN = function(x) score_matrix_transposed %*% matrix_x[,x] / cov_divisor
    )
    score_sd <- sapply(
      X = 1:ncol(score_matrix)
      , FUN = function(x) sd(score_matrix[,x])
    )
    # xSdVn - Numerical vector: variable standard deviations of the 'x' matrix
    xSdVn <- ropls_x@xSdVn
    result$correlation <- sapply(
      X = 1:ncol(matrix_x)
    , FUN = function(x) {
        ( score_matrix_transposed / score_sd ) %*% ( matrix_x[,x] / (xSdVn[x] * cov_divisor) )
      }
    )
  }
  result$correlation <- result$correlation[1,,drop = TRUE]
  result$covariance  <- result$covariance[1,,drop = TRUE]
  # Variant 4 of Variable Influence on Projection for OPLS from Galindo_Prieto_2014
  #    Length = number of features; labels = feature identifiers.  (The same is true for $correlation and $covariance.)
  result$vip4p     <- ropls_x@vipVn
  result$vip4o     <- ropls_x@orthoVipVn
  # get the level names
  level_names      <- sort(levels(as.factor(ropls_x@suppLs$y)))
  feature_count    <- length(ropls_x@vipVn)
  result$level1    <- rep.int(x = level_names[1], times = feature_count)
  result$level2    <- rep.int(x = level_names[2], times = feature_count)
  #strF(result$covariance)
  #print(sprintf("sd(covariance) = %f; sd(correlation) = %f", sd(result$covariance), sd(result$correlation)))
  superresult <- list()
  superresult$tsv1 <- data.frame(
    featureID           = names(ropls_x@vipVn)
  , factorLevel1        = result$level1
  , factorLevel2        = result$level2
  , correlation         = result$correlation
  , covariance          = result$covariance
  , row.names           = NULL
  )
  print(superresult$tsv1)
  #print(sprintf("sd(superresult$tsv1$covariance) = %f; sd(superresult$tsv1$correlation) = %f", sd(superresult$tsv1$covariance), sd(superresult$tsv1$correlation)))
  #strF(superresult$tsv1[,1:3])
  #strF(superresult$tsv1[,4:5])
  # Include thise in case future consumers of this routine want to use it in currently unanticipated ways
  result$oplsda    <- ropls_x          
  result$predictor <- ropls_x@suppLs$y   # in case future consumers of this routine want to use it in currently unanticipated ways
  return (result)
}

# # Wiklund_2008 centers and pareto-scales data before OPLS-DA S-plot
# cdm <- center_colmeans(dataMatrix)
# my_scale <- sqrt(apply(cdm, 2, sd, na.rm=TRUE))
# scdm <- sweep(cdm, 2, my_scale, "/")
# 
# my_matrix <- scdm
# #my_matrix <- dataMatrix
# 
# my_pca <- opls(my_matrix)
# dev.off()
# with( 
#   cor_vs_cov(
#     matrix_x = my_matrix
#   , ropls_x = my_pca
#   )
# , plot(y = correlation, x = covariance)
# )
# 
# # to suppress the summary plot, invoke with plotL = FALSE
# my_oplsda <- opls(my_matrix, sampleMetadata[, "gender"], algoC = algoC, predI = 1, orthoI = 1)
# dev.off()
# my_cor_vs_cov <-
#   cor_vs_cov(
#     matrix_x = my_matrix
#     , ropls_x = my_oplsda
#   )
# with(
#   my_cor_vs_cov
# , {
#     plot(y = correlation, x = covariance, type="n")
#     text(y = correlation, x = covariance)
#   }
# )
# plot( y = -cor(my_matrix, sampleMetadata[, "gender"]=="M"), x = -cov(sampleMetadata[, "gender"]=="M", my_matrix), type="n" )
# text( y = -cor(my_matrix, sampleMetadata[, "gender"]=="M"), x = -cov(sampleMetadata[, "gender"]=="M", my_matrix) )
# plot(y = my_cor_vs_cov$correlation, x = cor(my_matrix, sampleMetadata[, "gender"]=="M"))
# plot(y = my_cor_vs_cov$covariance, x = cov(my_matrix, sampleMetadata[, "gender"]=="M"))

