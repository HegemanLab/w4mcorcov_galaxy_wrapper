# compute and output detail plots
do_detail_plot <- function(
  x_dataMatrix
, x_predictor
, x_is_match
, x_algorithm
, x_prefix
, x_show_labels
, x_progress = print
, x_env
, x_crossval_i
) {
  off <- function(x) if (x_show_labels == "0") 0 else x
  if ( x_is_match
      && ncol(x_dataMatrix) > 0
      && length(unique(x_predictor)) > 1
      && x_crossval_i < nrow(x_dataMatrix)
  ) {
    my_oplsda <- opls(
        x      = x_dataMatrix
      , y      = x_predictor
      , algoC  = x_algorithm
      , predI  = 1
      , orthoI = if (ncol(x_dataMatrix) > 1) 1 else 0
      , printL = FALSE
      , plotL  = FALSE
      , crossvalI = x_crossval_i
      , scaleC = "pareto" # data centered and pareto scaled here only. This line fixes issue #2.
      )
    # strip out variables having negligible variance
    x_dataMatrix <- x_dataMatrix[ , names(my_oplsda@vipVn), drop = FALSE]
    my_oplsda_suppLs_y_levels <- levels(as.factor(my_oplsda@suppLs$y))

    fctr_lvl_1 <- my_oplsda_suppLs_y_levels[1]
    fctr_lvl_2 <- my_oplsda_suppLs_y_levels[2]
    do_s_plot <- function(
      x_env
    , predictor_projection_x   = TRUE
    , cplot_x      = FALSE
    , cor_vs_cov_x = NULL
    ) {
      if (cplot_x) {
        cplot_y_correlation <- (x_env$cplot_y == "correlation")
        default_lim_x <- 10
      } else {
        default_lim_x <- 1.2
      }
      if (is.null(cor_vs_cov_x)) {
        my_cor_vs_cov <- cor_vs_cov(
            matrix_x   = x_dataMatrix
          , ropls_x    = my_oplsda
          , predictor_projection_x = predictor_projection_x
          , x_progress
          , x_env
          )
      } else {
        my_cor_vs_cov <- cor_vs_cov_x
      }

      if (is.null(my_cor_vs_cov) || sum(!is.na(my_cor_vs_cov$tsv1$covariance)) < 2) {
        if (is.null(cor_vs_cov_x)) {
          x_progress(
            sprintf("No cor_vs_cov data produced for level %s versus %s", fctr_lvl_1, fctr_lvl_2)
          )
        }
        plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
        text(x=1, y=1, labels="too few covariance data")
        return(my_cor_vs_cov)
      }
      with(
        my_cor_vs_cov
      , {
          min_x <- min(covariance, na.rm = TRUE)
          max_x <- max(covariance, na.rm = TRUE)
          lim_x <- max(sapply(X=c(min_x, max_x), FUN=abs))

          # Regarding using VIP as a guide to selecting a biomarker:
          #   "It is generally accepted that a variable should be selected if vj>1, [27–29],
          #   but a proper threshold between 0.83 and 1.21 can yield more relevant variables according to [28]."
          #   (Mehmood 2012 doi:10.1186/1748-7188-6-27)
          plus_cor <- correlation
          plus_cov <- covariance
          if (length(plus_cor) != 0 && length(plus_cor) != 0) {
            cex <- 0.65
            if (projection == 1) {
              # predictor-projection
              vipcp <- pmax(0, pmin(1, (vip4p - 0.83) / (1.21 - 0.83)))
              if (!cplot_x) {
                # S-plot predictor-projection
                my_xlab <- "covariance(feature,t1)"
                my_x <- plus_cov
                my_ylab <- "correlation(feature,t1)"
                my_y <- plus_cor
                # X,Y limits for S-PLOT
                my_xlim <- c( -lim_x, lim_x ) * (1.0 + off(0.3))
                my_ylim <- c( -1.0, 1.0 ) * (1.0 + off(0.2) )
              } else {
                # C-plot predictor-projection
                my_xlab <- "variable importance in predictor-projection"
                my_x <- vip4p
                if (cplot_y_correlation) {
                  my_ylab <- "correlation(feature,t1)"
                  my_y <- plus_cor
                  my_ylim <- c( -1.0, 1.0 ) * (1.0 + off(0.2) )
                } else {
                  my_ylab <- "covariance(feature,t1)"
                  my_y <- plus_cov
                  my_ylim <- max(abs(plus_cov))
                  my_ylim <- c( -my_ylim, my_ylim ) * (1.0 + off(0.2) )
                }
                # X,Y limits for C-plot
                lim_x <- max(my_x, na.rm = TRUE) * 1.1
                lim_x <- min(lim_x, default_lim_x)
                my_xlim <- c( 0, lim_x ) # + off(0.2) )
              }
              my_load_distal <- loadp
              my_load_proximal <- loado
              red  <- as.numeric(correlation > 0) * vipcp
              blue <- as.numeric(correlation < 0) * vipcp
              alpha <- 0.1 + 0.4 * vipcp
              red[is.na(red)] <- 0
              blue[is.na(blue)] <- 0
              alpha[is.na(alpha)] <- 0
              my_col <- rgb(blue = blue, red = red, green = 0, alpha = alpha)
              main_label <- sprintf("%s for level %s versus %s"
                                   , x_prefix, fctr_lvl_1, fctr_lvl_2)
            } else {
              # orthogonal projection
              vipco <- pmax(0, pmin(1, (vip4o - 0.83) / (1.21 - 0.83)))
              if (!cplot_x) {
                # S-PLOT orthogonal-projection
                my_xlab <- "covariance(feature,to1)"
                my_x <- -plus_cov
                # X,Y limits for S-PLOT
                my_xlim <- c( -lim_x, lim_x ) * (1.0 + off(0.3))
                my_ylab <- "correlation(feature,to1)"
                my_y <- plus_cor
                my_ylim <- c( -1.0, 1.0 ) * (1.0 + off(0.2) )
              } else {
                # C-plot orthogonal-projection
                my_xlab <- "variable importance in orthogonal projection"
                my_x <- vip4o
                # C-plot orthogonal projection
                # X,Y limits for C-plot
                lim_x <- max(my_x, na.rm = TRUE) * 1.1
                my_xlim <- c( 0, lim_x ) # + off(0.2) )
                if (cplot_y_correlation) {
                  my_ylab <- "correlation(feature,to1)"
                  my_y <- plus_cor
                  my_ylim <- c( -1.0, 1.0 ) * (1.0 + off(0.2) )
                } else {
                  my_ylab <- "covariance(feature,to1)"
                  my_y <- plus_cov
                  my_ylim <- max(abs(plus_cov))
                  my_ylim <- c( -my_ylim, my_ylim ) * (1.0 + off(0.2) )
                }
              }
              my_load_distal <- loado
              my_load_proximal <- loadp
              alpha <- 0.1 + 0.4 * vipco
              alpha[is.na(alpha)] <- 0
              my_col <- rgb(blue = 0, red = 0, green = 0, alpha = alpha)
              main_label <- sprintf(
                "Features influencing orthogonal projection for %s versus %s"
              , fctr_lvl_1, fctr_lvl_2)
            }
            main_cex <- min(1.0, 46.0/nchar(main_label))
            my_feature_label_slant <- -30 # slant feature labels 30 degrees downward
            my_pch <- sapply(X = cor_p_value, function(x) if (x < 0.01) 16 else if (x < 0.05) 17 else 18)
            if ( sum(is.infinite(my_xlim)) == 0 ) {
              plot(
                y = my_y
              , x = my_x
              , type = "p"
              , xlim = my_xlim
              , ylim = my_ylim
              , xlab = my_xlab
              , ylab = my_ylab
              , main = main_label
              , cex.main = main_cex
              , cex = cex
              , pch = my_pch
              , col = my_col
              )
              low_x <- -0.7 * lim_x
              high_x <- 0.7 * lim_x
              if (projection == 1 && !cplot_x) {
                text(x = low_x, y = -0.05, labels =  fctr_lvl_1, col = "blue")
                text(x = high_x, y = 0.05, labels =  fctr_lvl_2, col = "red")
              }
              if ( x_show_labels != "0" ) {
                names(my_load_distal) <- tsv1$featureID
                names(my_load_proximal) <- tsv1$featureID
                if ( x_show_labels == "ALL" ) {
                  n_labels <- length(my_load_distal)
                } else {
                  n_labels <- as.numeric(x_show_labels)
                }
                n_labels <- min( n_labels, (1 + length(my_load_distal)) / 2 )
                labels_to_show <- c(
                  names(head(sort(my_load_distal), n = n_labels))
                , names(tail(sort(my_load_distal), n = n_labels))
                )
                labels <- unname(
                  sapply(
                    X = tsv1$featureID
                  , FUN = function(x) if ( x %in% labels_to_show ) x else ""
                  )
                )
                x_text_offset <- 0.024
                y_text_off <- 0.017
                if (!cplot_x) {
                  # S-plot
                  y_text_offset <- if (projection == 1) -y_text_off else y_text_off
                } else {
                  # C-plot
                  y_text_offset <-
                    sapply(
                      X = (my_y > 0)
                    , FUN = function(x) {
                        if (x) y_text_off else -y_text_off
                      }
                    )
                }
                label_features <- function(x_arg, y_arg, labels_arg, slant_arg) {
                  if (length(labels_arg) > 0) {
                    unique_slant <- unique(slant_arg)
                    if (length(unique_slant) == 1) {
                      text(
                        y = y_arg
                      , x = x_arg + x_text_offset
                      , cex = 0.4
                      , labels = labels_arg
                      , col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5) # grey semi-transparent labels
                      , srt = slant_arg
                      , adj = 0   # left-justified
                      )
                    } else {
                      for (slant in unique_slant) {
                        text(
                          y = y_arg[slant_arg == slant]
                        , x = x_arg[slant_arg == slant] + x_text_offset
                        , cex = 0.4
                        , labels = labels_arg[slant_arg == slant]
                        , col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5) # grey semi-transparent labels
                        , srt = slant
                        , adj = 0   # left-justified
                        )
                      }
                    }
                  }
                }
                if (!cplot_x) {
                  my_slant <- (if (projection == 1) 1 else -1) * my_feature_label_slant
                } else {
                  my_slant <- sapply(
                                X = (my_y > 0)
                              , FUN = function(x) if (x) 2 else -2
                              ) * my_feature_label_slant
                }
                if (length(my_x) > 1) {
                  label_features(
                    x_arg      = my_x  [my_x > 0]
                  , y_arg      = my_y  [my_x > 0] - y_text_offset
                  , labels_arg = labels[my_x > 0]
                  , slant_arg = (if (!cplot_x) -my_slant else (my_slant))
                  )
                  if (!cplot_x) {
                    label_features(
                      x_arg      = my_x  [my_x < 0]
                    , y_arg      = my_y  [my_x < 0] + y_text_offset
                    , labels_arg = labels[my_x < 0]
                    , slant_arg = my_slant
                    )
                  }
                } else {
                  if (!cplot_x) {
                    my_slant <- (if (my_x > 1) -1 else 1) * my_slant
                    my_y_arg <- my_y + (if (my_x > 1) -1 else 1) * y_text_offset
                  } else {
                    my_slant <- (if (my_y > 1) -1 else 1) * my_slant
                    my_y_arg <- my_y + (if (my_y > 1) -1 else 1) * y_text_offset
                  }
                  label_features(
                    x_arg = my_x
                  , y_arg = my_y_arg
                  , labels_arg = labels
                  , slant_arg = my_slant
                  )
                } # end if (length(my_x) > 1)
              } # end if ( x_show_labels != "0" )
            } else {
              plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
              text(x=1, y=1, labels="no S-plot is possible")
            } # end if (sum(is.infinte(my_xlim))==0)
          } else {
            plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
            text(x=1, y=1, labels="no S-plot is possible")
          } # end if (length(plus_cor) != 0 && length(plus_cor) != 0)
        } # end action
      ) # end with( my_cor_vs_cov, action )
      return (my_cor_vs_cov)
    } # end function do_s_plot
    my_cor_vs_cov <- do_s_plot(
      x_env = x_env
    , predictor_projection_x = TRUE
    , cplot_x = FALSE
    )
    typevc <- c("correlation",      # 1
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
               )                    # [c(3,8,9)] # [c(4,3,8,9)]
    if ( length(my_oplsda@orthoVipVn) > 0 ) {
      my_typevc <- typevc[c(9,3,8)]
    } else {
      my_typevc <- c("(dummy)","overview","(dummy)")
    }
    my_ortho_cor_vs_cov_exists <- FALSE
    for (my_type in my_typevc) {
      if (my_type %in% typevc) {
        tryCatch({
            if ( my_type != "x-loading" ) {
               plot(
                 x            = my_oplsda
               , typeVc       = my_type
               , parCexN      = 0.4
               , parDevNewL   = FALSE
               , parLayL      = TRUE
               , parEllipsesL = TRUE
               )
               if (my_type == "overview") {
                 sub_label <- sprintf("%s versus %s", fctr_lvl_1, fctr_lvl_2)
                 title(sub = sub_label)
               }
            } else {
              my_ortho_cor_vs_cov <- do_s_plot(
                x_env = x_env
              , predictor_projection_x = FALSE
              , cplot_x = FALSE
              )
              my_ortho_cor_vs_cov_exists <- TRUE
            }
          }
        , error = function(e) {
          x_progress(
            sprintf(
              "factor level %s or %s may have only one sample - %s"
            , fctr_lvl_1
            , fctr_lvl_2
            , e$message
            )
          )
        })
      } else {
        plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
        text(x=1, y=1, labels="no orthogonal projection is possible")
      }
    }
    cplot_p <- x_env$cplot_p
    cplot_o <- x_env$cplot_o
    if (cplot_p || cplot_o) {
      if (cplot_p) {
        if (!is.null(my_cor_vs_cov)) {
          do_s_plot(
            x_env = x_env
          , predictor_projection_x = TRUE
          , cplot_x = TRUE
          , cor_vs_cov_x = my_cor_vs_cov
          )
        } else {
          plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
          text(x=1, y=1, labels="no predictor projection is possible")
        }
        did_plots <- 1
      } else {
        did_plots <- 0
      }
      if (cplot_o) {
        if (my_ortho_cor_vs_cov_exists) {
          do_s_plot(
            x_env = x_env
          , predictor_projection_x = FALSE
          , cplot_x = TRUE
          , cor_vs_cov_x = my_ortho_cor_vs_cov
          )
        } else {
          plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
          text(x=1, y=1, labels="no orthogonal projection is possible")
        }
        did_plots <- 1 + did_plots
      }
      if (did_plots == 1) {
        plot(x=1, y=1, xaxt="n", yaxt="n", xlab="", ylab="", type="n", fg = "white")
      }
    }
    return (my_cor_vs_cov)
  } else {
    return (NULL)
  }
}

# S-PLOT and OPLS reference: Wiklund_2008 doi:10.1021/ac0713510
corcov_calc <- function(
  calc_env
, failure_action      = stop
, progress_action     = function(x) { }
, corcov_tsv_action   = function(t) { }
, salience_tsv_action = function(t) { }
, extra_plots         = c()
) {
  if ( ! is.environment(calc_env) ) {
    failure_action("corcov_calc: fatal error - 'calc_env' is not an environment")
    return ( FALSE )
  }
  if ( ! is.function(corcov_tsv_action) ) {
    failure_action("corcov_calc: fatal error - 'corcov_tsv_action' is not a function")
    return ( FALSE )
  }
  if ( ! is.function(salience_tsv_action) ) {
    failure_action("salience_calc: fatal error - 'salience_tsv_action' is not a function")
    return ( FALSE )
  }

  # extract parameters from the environment
  vrbl_metadata <- calc_env$vrbl_metadata
  vrbl_metadata_names <- vrbl_metadata[, 1]
  smpl_metadata <- calc_env$smpl_metadata
  data_matrix <- calc_env$data_matrix
  pair_significant_features_only <- calc_env$pairSigFeatOnly
  facC <- calc_env$facC
  tesC <- calc_env$tesC
  # extract the levels from the environment
  originalLevCSV <- levCSV <- calc_env$levCSV
  # matchingC is one of { "none", "wildcard", "regex" }
  matchingC <- calc_env$matchingC
  labelFeatures <- calc_env$labelFeatures

  # arg/env checking
  if (!(facC %in% names(smpl_metadata))) {
    failure_action(
      sprintf("bad parameter!  Factor name '%s' not found in sampleMetadata"
             , facC))
    return ( FALSE )
  }

  mz             <- vrbl_metadata$mz
  names(mz)      <- vrbl_metadata$variableMetadata
  mz_lookup      <- function(feature) unname(mz[feature])

  rt             <- vrbl_metadata$rt
  names(rt)      <- vrbl_metadata$variableMetadata
  rt_lookup      <- function(feature) unname(rt[feature])

  # calculate salience_df as data.frame(feature, max_level, max_median, salient_rcv, mean_median, salience, salient_rcv)
  salience_df <- calc_env$salience_df <- w4msalience(
    data_matrix    = data_matrix
  , sample_class   = smpl_metadata[,facC]
  , failure_action = failure_action
  )
  salience_tsv_action({
    with (
      salience_df
    , {
      my_df <<- data.frame(
        featureID       = feature
      , salientLevel    = max_level
      , salientRCV      = salient_rcv
      , relativeSalientDistance = relative_salient_distance
      , salience        = salience
      , mz              = mz_lookup(feature)
      , rt              = rt_lookup(feature)
      )
    })
    my_df[order(-my_df$relativeSalientDistance),]
  })

  # transform wildcards to regexen
  if (matchingC == "wildcard") {
    # strsplit(x = "hello,wild,world", split = ",")
    levCSV <- gsub("[.]", "[.]", levCSV)
    levCSV <- strsplit(x = levCSV, split = ",")
    levCSV <- sapply(levCSV, utils::glob2rx, trim.tail = FALSE)
    levCSV <- paste(levCSV, collapse = ",")
  }
  # function to determine whether level is a member of levCSV
  isLevelSelected <- function(lvl) {
    matchFun <- if (matchingC != "none") grepl else `==`
    return(
      Reduce(
        f = "||"
      , x = sapply(
              X = strsplit(
                    x = levCSV
                  , split = ","
                  , fixed = TRUE
                  )[[1]]
            , FUN = matchFun
            , lvl
            )
      )
    )
  }

  # transpose matrix because ropls matrix is the transpose of XCMS matrix
  tdm <- t(data_matrix)
  # Wiklund_2008 centers and pareto-scales data before OPLS-DA S-plot
  # However, data should be neither centered nor pareto scaled here because ropls::opls does that; this fixes issue #2.

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


  did_plot <- FALSE
  if (tesC != "none") {
    # for each column name, extract the parts of the name matched by 'col_pattern', if any
    the_colnames <- colnames(vrbl_metadata)
    if ( ! Reduce( f = "||", x = grepl(tesC, the_colnames) ) ) {
      failure_action(
        sprintf(
          "bad parameter!  variableMetadata must contain results of W4M Univariate test '%s'."
        , tesC))
      return ( FALSE )
    }
    col_matches <- lapply(
      X = the_colnames,
      FUN = function(x) {
        regmatches( x, regexec(col_pattern, x) )[[1]]
      }
    )
    ## first contrast each selected level with all other levels combined into one "super-level" ##
    # process columns matching the pattern
    level_union <- c()
    for (i in 1:length(col_matches)) {
      col_match <- col_matches[[i]]
      if (length(col_match) > 0) {
        # it's an actual match; extract the pieces, e.g., k10_kruskal_k4.k3_sig
        vrbl_metadata_col <- col_match[1]               # ^^^^^^^^^^^^^^^^^^^^^  # Column name
        fctr_lvl_1        <- col_match[2]               #             ^^         # Factor-level 1
        fctr_lvl_2        <- col_match[3]               #                ^^      # Factor-level 2
        # only process this column if both factors are members of lvlCSV
        is_match <- isLevelSelected(fctr_lvl_1) && isLevelSelected(fctr_lvl_2)
        if (is_match) {
          level_union <- c(level_union, col_match[2], col_match[3])
        }
      }
    }
    level_union <- unique(sort(level_union))
    overall_significant <- 1 == (
      if (intersample_sig_col %in% colnames(vrbl_metadata)) {
        vrbl_metadata[,intersample_sig_col]
      } else {
        1
      }
    )
    if ( length(level_union) > 2 ) {
      chosen_samples <- smpl_metadata_facC %in% level_union
      chosen_facC <- as.character(smpl_metadata_facC[chosen_samples])
      col_selector <- vrbl_metadata_names[ overall_significant ]
      my_matrix <- tdm[ chosen_samples, col_selector, drop = FALSE ]
      plot_action <- function(fctr_lvl_1, fctr_lvl_2) {
        progress_action(
          sprintf("calculating/plotting contrast of %s vs. %s"
                 , fctr_lvl_1, fctr_lvl_2)
        )
        predictor <- sapply(
          X = chosen_facC
        , FUN = function(fac) if ( fac == fctr_lvl_1 ) fctr_lvl_1 else fctr_lvl_2
        )
        my_cor_cov <- do_detail_plot(
          x_dataMatrix  = my_matrix
        , x_predictor   = predictor
        , x_is_match    = TRUE
        , x_algorithm   = "nipals"
        , x_prefix      = if (pair_significant_features_only) {
                            "Significantly contrasting features"
                          } else {
                            "Significant features"
                          }
        , x_show_labels = labelFeatures
        , x_progress    = progress_action
        , x_crossval_i  = min(7, length(chosen_samples))
        , x_env         = calc_env
        )
        if ( is.null(my_cor_cov) ) {
          progress_action("NOTHING TO PLOT")
        } else {
          my_tsv <- my_cor_cov$tsv1
          my_tsv$mz <- mz_lookup(my_tsv$featureID)
          my_tsv$rt <- rt_lookup(my_tsv$featureID)
          my_tsv["level1Level2Sig"] <- vrbl_metadata[
            match(my_tsv$featureID, vrbl_metadata_names)
          , vrbl_metadata_col
          ]
          tsv <<- my_tsv
          corcov_tsv_action(tsv)
          did_plot <<- TRUE
        }
      }
      if ( length(level_union) != 2 ) {
        fctr_lvl_2 <- "other"
        for ( fctr_lvl_1 in level_union[1:length(level_union)] ) {
          plot_action(fctr_lvl_1, fctr_lvl_2)
        }
      } else {
        plot_action(fctr_lvl_1 = level_union[1], fctr_lvl_2 = level_union[2])
      }
    }

    if ( length(level_union) > 1 ) {
      ## next, contrast each selected level with each of the other levels individually ##
      # process columns matching the pattern
      for (i in 1:length(col_matches)) {
        # for each potential match of the pattern
        col_match <- col_matches[[i]]
        if (length(col_match) > 0) {
          # it's an actual match; extract the pieces, e.g., k10_kruskal_k4.k3_sig
          vrbl_metadata_col <- col_match[1]               # ^^^^^^^^^^^^^^^^^^^^^  # Column name
          fctr_lvl_1        <- col_match[2]               #             ^^         # Factor-level 1
          fctr_lvl_2        <- col_match[3]               #                ^^      # Factor-level 2
          # only process this column if both factors are members of lvlCSV
          is_match <- isLevelSelected(fctr_lvl_1) && isLevelSelected(fctr_lvl_2)
          if (is_match) {
            progress_action(
              sprintf("calculating/plotting contrast of %s vs. %s."
                     , fctr_lvl_1, fctr_lvl_2
              )
            )
            # choose only samples with one of the two factors for this column
            chosen_samples <- smpl_metadata_facC %in% c(fctr_lvl_1, fctr_lvl_2)
            predictor <- smpl_metadata_facC[chosen_samples]
            # extract only the significantly-varying features and the chosen samples
            fully_significant   <- 1 == vrbl_metadata[,vrbl_metadata_col] *
              ( if (intersample_sig_col %in% colnames(vrbl_metadata)) {
                  vrbl_metadata[,intersample_sig_col]
                } else {
                  1
                }
              )
            col_selector <- vrbl_metadata_names[
              if (pair_significant_features_only) fully_significant else overall_significant
            ]
            my_matrix <- tdm[ chosen_samples, col_selector, drop = FALSE ]
            my_cor_cov <- do_detail_plot(
              x_dataMatrix  = my_matrix
            , x_predictor   = predictor
            , x_is_match    = is_match
            , x_algorithm   = "nipals"
            , x_prefix      = if (pair_significant_features_only) {
                                "Significantly contrasting features"
                              } else {
                                "Significant features"
                              }
            , x_show_labels = labelFeatures
            , x_progress    = progress_action
            , x_crossval_i  = min(7, length(chosen_samples))
            , x_env         = calc_env
            )
            if ( is.null(my_cor_cov) ) {
              progress_action("NOTHING TO PLOT.")
            } else {
              tsv <- my_cor_cov$tsv1
              tsv$mz <- mz_lookup(tsv$featureID)
              tsv$rt <- rt_lookup(tsv$featureID)
              tsv["level1Level2Sig"] <- vrbl_metadata[
                match(tsv$featureID, vrbl_metadata_names)
              , vrbl_metadata_col
              ]
              corcov_tsv_action(tsv)
              did_plot <- TRUE
            }
          } else {
            progress_action(
              sprintf("skipping contrast of %s vs. %s."
                     , fctr_lvl_1, fctr_lvl_2
              )
            )
          }
        }
      }
    }
  } else {
    # tesC == "none"
    # find all the levels for factor facC in sampleMetadata
    level_union <- unique(sort(smpl_metadata_facC))
    # identify the selected levels for factor facC from sampleMetadata
    level_include <- sapply(X = level_union, FUN = isLevelSelected)
    # discard the non-selected levels for factor facC
    level_union <- level_union[level_include]
    if ( length(level_union) > 1 ) {
      if ( length(level_union) > 2 ) {
        ## pass 1 - contrast each selected level with all other levels combined into one "super-level" ##
        completed <- c()
        lapply(
          X = level_union
        , FUN = function(x) {
            fctr_lvl_1 <- x[1]
            fctr_lvl_2 <- {
              if ( fctr_lvl_1 %in% completed )
                return("DUMMY")
              completed <<- c(completed, fctr_lvl_1)
              setdiff(level_union, fctr_lvl_1)
            }
            chosen_samples <- smpl_metadata_facC %in% c(fctr_lvl_1, fctr_lvl_2)
            fctr_lvl_2 <- "other"
            if (length(unique(chosen_samples)) < 1) {
              progress_action(
                sprintf("Skipping contrast of %s vs. %s; there are no chosen samples."
                , fctr_lvl_1, fctr_lvl_2)
              )
            } else {
              chosen_facC <- as.character(smpl_metadata_facC[chosen_samples])
              predictor <- sapply(
                X = chosen_facC
              , FUN = function(fac) {
                  if ( fac == fctr_lvl_1 ) fctr_lvl_1 else fctr_lvl_2
                }
              )
              my_matrix <- tdm[ chosen_samples, , drop = FALSE ]
              # only process this column if both factors are members of lvlCSV
              is_match <- isLevelSelected(fctr_lvl_1)
              if (is_match) {
                progress_action(
                  sprintf("Calculating/plotting contrast of %s vs. %s"
                  , fctr_lvl_1, fctr_lvl_2)
                )
                my_cor_cov <- do_detail_plot(
                  x_dataMatrix  = my_matrix
                , x_predictor   = predictor
                , x_is_match    = is_match
                , x_algorithm   = "nipals"
                , x_prefix      = "Features"
                , x_show_labels = labelFeatures
                , x_progress    = progress_action
                , x_crossval_i  = min(7, length(chosen_samples))
                , x_env         = calc_env
                )
                if ( is.null(my_cor_cov) ) {
                  progress_action("NOTHING TO PLOT...")
                } else {
                  tsv <- my_cor_cov$tsv1
                  tsv$mz <- mz_lookup(tsv$featureID)
                  tsv$rt <- rt_lookup(tsv$featureID)
                  corcov_tsv_action(tsv)
                  did_plot <<- TRUE
                }
              } else {
              }
            }
            "dummy" # need to return a value; otherwise combn fails with an error
          }
        )
      }
      ## pass 2 - contrast each selected level with each of the other levels individually ##
      completed <- c()
      utils::combn(
        x = level_union
      , m = 2
      , FUN = function(x) {
          fctr_lvl_1 <- x[1]
          fctr_lvl_2 <- x[2]
          chosen_samples <- smpl_metadata_facC %in% c(fctr_lvl_1, fctr_lvl_2)
          if (length(unique(chosen_samples)) < 1) {
            progress_action(
              sprintf("Skipping contrast of %s vs. %s. - There are no chosen samples."
                     , fctr_lvl_1, fctr_lvl_2
              )
            )
          } else {
            chosen_facC <- as.character(smpl_metadata_facC[chosen_samples])
            predictor <- chosen_facC
            my_matrix <- tdm[ chosen_samples, , drop = FALSE ]
            # only process this column if both factors are members of lvlCSV
            is_match <- isLevelSelected(fctr_lvl_1) && isLevelSelected(fctr_lvl_2)
            if (is_match) {
              progress_action(
                sprintf("Calculating/plotting contrast of %s vs. %s."
                       , fctr_lvl_1, fctr_lvl_2)
                )
              my_cor_cov <- do_detail_plot(
                x_dataMatrix  = my_matrix
              , x_predictor   = predictor
              , x_is_match    = is_match
              , x_algorithm   = "nipals"
              , x_prefix      = "Features"
              , x_show_labels = labelFeatures
              , x_progress    = progress_action
              , x_crossval_i  = min(7, length(chosen_samples))
              , x_env         = calc_env
              )
              if ( is.null(my_cor_cov) ) {
                progress_action("NOTHING TO PLOT.....")
              } else {
                tsv <- my_cor_cov$tsv1
                tsv$mz <- mz_lookup(tsv$featureID)
                tsv$rt <- rt_lookup(tsv$featureID)
                corcov_tsv_action(tsv)
                did_plot <<- TRUE
              }
            } else {
              progress_action(
                sprintf("Skipping contrast of %s vs. %s."
                       , fctr_lvl_1, fctr_lvl_2
                )
              )
            }
          }
          "dummy" # need to return a value; otherwise combn fails with an error
        }
      )
    } else {
      progress_action("NOTHING TO PLOT......")
    }
  }
  if (!did_plot) {
    failure_action(
      sprintf(
        "Did not plot. Does sampleMetadata have at least two levels of factor '%s' matching '%s' ?"
      , facC, originalLevCSV))
    return ( FALSE )
  }
  return ( TRUE )
}

# Calculate data for correlation-versus-covariance plot
#   Adapted from:
#     Wiklund_2008 doi:10.1021/ac0713510
#     Galindo_Prieto_2014 doi:10.1002/cem.2627
#     https://github.com/HegemanLab/extra_tools/blob/master/generic_PCA.R
cor_vs_cov <- function(
  matrix_x
, ropls_x
, predictor_projection_x = TRUE
, x_progress = print
, x_env
) {
  tryCatch({
      return(
        cor_vs_cov_try( matrix_x, ropls_x, predictor_projection_x, x_progress, x_env)
      )
    }
  , error = function(e) {
    x_progress(
      sprintf(
        "cor_vs_cov fatal error - %s"
      , as.character(e) # e$message
      )
    )
    return ( NULL )
  })
}

cor_vs_cov_try <- function(
  matrix_x                      # rows are samples; columns, features
, ropls_x                       # an instance of ropls::opls
, predictor_projection_x = TRUE # TRUE for predictor projection; FALSE for orthogonal projection
, x_progress = print            # function to produce progress and error messages
, x_env
) {
  my_matrix_x <- matrix_x
  my_matrix_x[my_matrix_x==0] <- NA
  fdr_features <- x_env$fdr_features

  x_class <- class(ropls_x)
  if ( !( as.character(x_class) == "opls" ) ) {
    stop(
      paste(
        "cor_vs_cov: Expected ropls_x to be of class ropls::opls but instead it was of class "
      , as.character(x_class)
      )
    )
  }
  if ( !ropls_x@suppLs$algoC == "nipals" ) {
    # suppLs$algoC - Character: algorithm used - "svd" for singular value decomposition; "nipals" for NIPALS
    stop(
      paste(
        "cor_vs_cov: Expected ropls::opls instance to have been computed by the NIPALS algorithm rather than "
      , ropls_x@suppLs$algoC
      )
    )
  }
  result <- list()
  result$projection <- projection <- if (predictor_projection_x) 1 else 2

  # I used equations (1) and (2) from Wiklund 2008, doi:10.1021/ac0713510
  #   (and not from the supplement despite the statement that, for the NIPALS algorithm,
  #   the equations from the supplement should be used) because of the definition of the
  #   Pearson/Galton coefficient of correlation is defined as
  #   $$
  #      \rho_{X,Y}= \frac{\operatorname{cov}(X,Y)}{\sigma_X \sigma_Y}
  #   $$
  #   as described (among other places) on Wikipedia at
  #     https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_population
  # The equations in the supplement said to use, for the predictive component t1,
  #      \rho_{t1,X_i}= \frac{\operatorname{cov}(t1,X_i)}{(\operatorname{mag}(t1))(\operatorname{mag}(X_i))}
  # but the results that I got were dramatically different from published results for S-PLOTs;
  # perhaps my data are not centered exactly the same way that theirs were.
  # The correlations calculated here are in agreement with those calculated with the code from
  #   page 22 of https://cran.r-project.org/web/packages/muma/muma.pdf


  # count the features/variables (one column for each sample)
  # count the features/variables (one column for each sample)
  n_features <- ncol(my_matrix_x)
  all_n_features <- x_env$fdr_features
  if (length(grep("^[0-9][0-9]*$", all_n_features)) > 0) {
    all_n_features <- as.integer(all_n_features)
  } else {
    all_n_features <- n_features
  }
  fdr_n_features <- max(n_features, all_n_features)
  # print("n_features")
  # print(n_features)

  # count the samples/observations (one row for each sample)
  n_observations <- nrow(my_matrix_x)

  # choose whether to plot the predictive score vector or orthogonal score vector
  if (predictor_projection_x)
     score_vector <- ropls_x@scoreMN
  else
     score_vector <- ropls_x@orthoScoreMN

  # compute the covariance of each feature with the score vector
  result$covariance <-
    sapply(
      X = 1:n_features
    , FUN = function(i) {
        cov(score_vector, my_matrix_x[ , i, drop = TRUE], use = "pairwise.complete.obs")
      }
    )
  # access covariance by feature name
  names(result$covariance) <- colnames(my_matrix_x)

  # compute the correlation of each feature with the score vector
  result$correlation <-
    sapply(
      X = 1:n_features
    , FUN = function(i) {
        cor(score_vector, my_matrix_x[ , i, drop = TRUE], use = "pairwise.complete.obs")
      }
    )
  # access correlation by feature name
  names(result$correlation) <- colnames(my_matrix_x)

  # eliminate NAs in either correlation or covariance
  nas <- is.na(result$correlation) | is.na(result$covariance)
  featureID          <- names(ropls_x@vipVn)
  featureID          <- featureID[!nas]
  result$correlation <- result$correlation[!nas]
  result$covariance  <- result$covariance[!nas]
  n_features <- length(featureID)

  correl_pci <- lapply(
    X = 1:n_features
  , FUN = function(i) {
      correl.ci(
        r = result$correlation[i]
      , n_obs = n_observations
      , n_vars = fdr_n_features
      )
    }
  )
  result$p_value_raw <- sapply(
    X = 1:n_features
  , FUN = function(i) correl_pci[[i]]$p.value.raw
  )
  result$p_value_raw[is.na(result$p_value_raw)] <- 1.0
  result$p_value <- sapply(
    X = 1:n_features
  , FUN = function(i) correl_pci[[i]]$p.value
  )
  result$p_value[is.na(result$p_value)] <- 1.0
  result$ci_lower <- sapply(
    X = 1:n_features
  , FUN = function(i) correl_pci[[i]]$CI["lower"]
  )
  result$ci_upper <- sapply(
    X = 1:n_features
  , FUN = function(i) correl_pci[[i]]$CI["upper"]
  )


  # extract "variant 4 of Variable Influence on Projection for OPLS" (see Galindo_Prieto_2014, DOI 10.1002/cem.2627)
  #    Length = number of features; labels = feature identifiers.  (The same is true for $correlation and $covariance.)
  result$vip4p     <- as.numeric(ropls_x@vipVn)[!nas]
  result$vip4o     <- as.numeric(ropls_x@orthoVipVn)[!nas]
  if (length(result$vip4o) == 0) result$vip4o <- NA
  # extract the loadings
  result$loadp     <- as.numeric(ropls_x@loadingMN)[!nas]
  result$loado     <- as.numeric(ropls_x@orthoLoadingMN)[!nas]
  # get the level names
  level_names      <- sort(levels(as.factor(ropls_x@suppLs$y)))
  fctr_lvl_1       <- level_names[1]
  fctr_lvl_2       <- level_names[2]
  result$level1    <- rep.int(x = fctr_lvl_1, times = n_features)
  result$level2    <- rep.int(x = fctr_lvl_2, times = n_features)
  greaterLevel <- sapply(
    X = result$correlation
  , FUN = function(my_corr) {
      tryCatch({
          if ( is.na(my_corr) || ! is.numeric( my_corr ) ) {
            NA
          } else {
            if ( my_corr < 0 ) fctr_lvl_1 else fctr_lvl_2
          }
        }
      , error = function(e) {
          print(my_corr)
          x_progress(
            sprintf(
              "cor_vs_cov -> sapply:  error - substituting NA - %s"
            , as.character(e)
            )
          )
          NA
        }
      )
    }
  )

  # build a data frame to hold the content for the tab-separated values file
  tsv1 <- data.frame(
    featureID     = featureID
  , factorLevel1  = result$level1
  , factorLevel2  = result$level2
  , greaterLevel  = greaterLevel
  , projection    = result$projection
  , correlation   = result$correlation
  , covariance    = result$covariance
  , vip4p         = result$vip4p
  , vip4o         = result$vip4o
  , loadp         = result$loadp
  , loado         = result$loado
  , cor_p_val_raw = result$p_value_raw
  , cor_p_value   = result$p_value
  , cor_ci_lower  = result$ci_lower
  , cor_ci_upper  = result$ci_upper
  )
  rownames(tsv1) <- tsv1$featureID

  # build the superresult, i.e., the result returned by this function
  superresult <- list()
  superresult$projection <- result$projection
  superresult$covariance <- result$covariance
  superresult$correlation <- result$correlation
  superresult$vip4p <- result$vip4p
  superresult$vip4o <- result$vip4o
  superresult$loadp <- result$loadp
  superresult$loado <- result$loado
  superresult$cor_p_value <- tsv1$cor_p_value
  superresult$details <- result

  # remove any rows having NA for covariance or correlation
  tsv1 <- tsv1[!is.na(tsv1$correlation),]
  tsv1 <- tsv1[!is.na(tsv1$covariance),]
  superresult$tsv1 <- tsv1

  # # I did not include these but left them commentd out in case future
  # #   consumers of this routine want to use it in currently unanticipated ways
  # result$superresult <- superresult
  # result$oplsda    <- ropls_x
  # result$predictor <- ropls_x@suppLs$y

  return (superresult)
}

# Code for correl.ci was adapted from correl function from:
#   @book{
#     Tsagris_2018,
#     author = {Tsagris, Michail},
#     year = {2018},
#     link = {https://www.researchgate.net/publication/324363311_Multivariate_data_analysis_in_R},
#     title = {Multivariate data analysis in R}
#   }
# which follows
#   https://en.wikipedia.org/wiki/Fisher_transformation#Definition

correl.ci <- function(r, n_obs, n_vars, a = 0.05, rho = 0) {
  ## r is the calculated correlation coefficient for n_obs pairs of observations of one variable
  ## a is the significance level
  ## rho is the hypothesised correlation
  zh0 <- atanh(rho) # 0.5*log((1+rho)/(1-rho)), i.e., Fisher's z-transformation for Ho
  zh1 <- atanh(r)   # 0.5*log((1+r)/(1-r)), i.e., Fisher's z-transformation for H1
  se <- (1 - r^2)/sqrt(n_obs - 3) ## standard error for Fisher's z-transformation of Ho
  test <- (zh1 - zh0)/se ### test statistic
  pvalue <- 2*(1 - pnorm(abs(test))) ## p-value
  pvalue.adj <- p.adjust(p = pvalue, method = "BY", n = n_vars)
  z_L <- zh1 - qnorm(1 - a/2)*se
  z_h <- zh1 + qnorm(1 - a/2)*se
  fish_l <- tanh(z_L) # (exp(2*z_l)-1)/(exp(2*z_l)+1), i.e., lower confidence limit
  fish_h <- tanh(z_h) # (exp(2*z_h)-1)/(exp(2*z_h)+1), i.e., upper confidence limit
  ci <- c(fish_l, fish_h)
  names(ci) <- c("lower", "upper")
  list(
    correlation = r
  , p.value.raw = pvalue
  , p.value = pvalue.adj
  , CI = ci
  )
}

# vim: sw=2 ts=2 et ai :
