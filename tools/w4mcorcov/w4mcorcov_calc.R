# center with 'colMeans()' - ref: http://gastonsanchez.com/visually-enforced/how-to/2014/01/15/Center-data-in-R/
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

#### OPLS-DA
algoC <- "nipals"

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
      && length(unique(x_predictor))> 1
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
    x_dataMatrix <- x_dataMatrix[,names(my_oplsda@vipVn), drop = FALSE]
    my_oplsda_suppLs_y_levels <- levels(as.factor(my_oplsda@suppLs$y))
    # x_progress(strF(x_dataMatrix))
    # x_progress(strF(my_oplsda))
    #x_progress(head(my_oplsda_suppLs_y_levels))
    #x_progress(unique(my_oplsda_suppLs_y_levels))
    fctr_lvl_1 <- my_oplsda_suppLs_y_levels[1]
    fctr_lvl_2 <- my_oplsda_suppLs_y_levels[2]
    do_s_plot <- function(
        x_env
      , predictor_projection_x   = TRUE
      , cplot_x      = FALSE
      , cor_vs_cov_x = NULL
      )
    {
      #print(ls(x_env))               # "cplot_y" etc
      #print(str(x_env$cplot_y))      # chr "covariance"
      if (cplot_x) {
        #print(x_env$cplot_y)         # "covariance"
        cplot_y_correlation <- (x_env$cplot_y == "correlation")
        #print(cplot_y_correlation)   # FALSE
      }
      if (is.null(cor_vs_cov_x)) {
        my_cor_vs_cov <- cor_vs_cov(
            matrix_x   = x_dataMatrix
          , ropls_x    = my_oplsda
          , predictor_projection_x = predictor_projection_x
          , x_progress
          )
      } else {
        my_cor_vs_cov <- cor_vs_cov_x
      }
      # print("str(my_cor_vs_cov)")
      # str(my_cor_vs_cov)
      if (is.null(my_cor_vs_cov) || sum(!is.na(my_cor_vs_cov$tsv1$covariance)) < 2) {
        if (is.null(cor_vs_cov_x)) {
          x_progress("No cor_vs_cov data produced")
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
          covariance <- covariance / lim_x
          lim_x <- 1.2
          # "It is generally accepted that a variable should be selected if vj>1, [27–29],
          #   but a proper threshold between 0.83 and 1.21 can yield more relevant variables according to [28]."
          #   (Mehmood 2012 doi:10.1186/1748-7188-6-27)
          plus_cor <- correlation
          plus_cov <- covariance
          cex <- 0.65
          which_projection <- if (projection == 1) "t1" else "o1"
          which_loading <- if (projection == 1) "parallel" else "orthogonal"
          if (projection == 1) { # predictor-projection
            vipcp <- pmax(0, pmin(1,(vip4p-0.83)/(1.21-0.83)))
            if (!cplot_x) { # S-plot predictor-projection
              my_xlab <- "relative covariance(feature,t1)"
              my_x <- plus_cov
              my_ylab <- "correlation(feature,t1)"
              my_y <- plus_cor
            } else { # C-plot predictor-projection
              my_xlab <- "variable importance in predictor-projection"
              my_x <- vip4p
              if (cplot_y_correlation) {
                my_ylab <- "correlation(feature,t1)"
                my_y <- plus_cor
              } else {
                my_ylab <- "relative covariance(feature,t1)"
                my_y <- plus_cov
              }
            }
            if (cplot_x) {
              lim_x <- max(my_x, na.rm = TRUE) * 1.1
              my_xlim <- c( 0, lim_x + off(0.2) )
            } else {
              my_xlim <- c( -lim_x - off(0.2), lim_x + off(0.2) )
            }
            my_ylim <- c( -1.0   - off(0.2), 1.0   + off(0.2) )
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
          } else { # orthogonal projection
            vipco <- pmax(0, pmin(1,(vip4o-0.83)/(1.21-0.83)))
            if (!cplot_x) {
              my_xlab <- "relative covariance(feature,to1)"
              my_x <- -plus_cov
            } else {
              my_xlab <- "variable importance in orthogonal projection"
              my_x <- vip4o
            }
            if (!cplot_x) { # S-plot orthogonal projection
              my_xlim <- c( -lim_x - off(0.2), lim_x + off(0.2) )
              my_ylab <- "correlation(feature,to1)"
              my_y <- plus_cor
            } else { # C-plot orthogonal projection
              lim_x <- max(my_x, na.rm = TRUE) * 1.1
              my_xlim <- c( 0, lim_x + off(0.2) )
              if (cplot_y_correlation) {
                my_ylab <- "correlation(feature,to1)"
                my_y <- plus_cor
              } else {
                my_ylab <- "relative covariance(feature,to1)"
                my_y <- plus_cov
              }
            }
            my_ylim <- c( -1.0   - off(0.2), 1.0   + off(0.2) )
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
          , pch = 16
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
              names(head(sort(my_load_distal),n = n_labels))
            , names(tail(sort(my_load_distal),n = n_labels))
            )
            labels <- unname(
              sapply(
                X = tsv1$featureID
              , FUN = function(x) if( x %in% labels_to_show ) x else ""
              )
            )
            x_text_offset <- 0.024
            y_text_off <- 0.017
            if (!cplot_x) { # S-plot
              y_text_offset <- if (projection == 1) -y_text_off else y_text_off
            } else { # C-plot
              y_text_offset <-
                sapply(
                  X = (my_y > 0)
                , FUN = function(x) { if (x) y_text_off else -y_text_off }
                )
            }
            label_features <- function(x_arg, y_arg, labels_arg, slant_arg) {
              # print("str(x_arg)")
              # print(str(x_arg))
              # print("str(y_arg)")
              # print(str(y_arg))
              # print("str(labels_arg)")
              # print(str(labels_arg))
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
                my_y_arg = my_y + (if (my_x > 1) -1 else 1) * y_text_offset
              } else {
                my_slant <- (if (my_y > 1) -1 else 1) * my_slant
                my_y_arg = my_y + (if (my_y > 1) -1 else 1) * y_text_offset
              }
              label_features(
                x_arg = my_x
              , y_arg = my_y_arg
              , labels_arg = labels
              , slant_arg = my_slant
              )
            }
          }
        }
      )
      return (my_cor_vs_cov)
    }
    my_cor_vs_cov <- do_s_plot(
      x_env = x_env
    , predictor_projection_x = TRUE
    , cplot_x = FALSE
    )
    typeVc <- c("correlation",      # 1
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
      my_typevc <- typeVc[c(9,3,8)]
    } else {
      my_typevc <- c("(dummy)","overview","(dummy)")
    }
    my_ortho_cor_vs_cov_exists <- FALSE
    for (my_type in my_typevc) {
      if (my_type %in% typeVc) {
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
        }, error = function(e) {
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
        do_s_plot(
          x_env = x_env
        , predictor_projection_x = TRUE
        , cplot_x = TRUE
        , cor_vs_cov_x = my_cor_vs_cov
        )
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
  vrbl_metadata_names <- vrbl_metadata[,1]
  smpl_metadata <- calc_env$smpl_metadata
  data_matrix <- calc_env$data_matrix
  pairSigFeatOnly <- calc_env$pairSigFeatOnly
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

  # calculate salience_df as data.frame(feature, max_level, max_median, max_rcv, mean_median, salience, salient_rcv)
  salience_df <- calc_env$salience_df <- w4msalience(
    data_matrix    = data_matrix
  , sample_class   = smpl_metadata[,facC]
  , failure_action = failure_action
  )
  salience_tsv_action({
    my_df <- data.frame(
      featureID    = salience_df$feature
    , salientLevel = salience_df$max_level
    , salientRCV   = salience_df$salient_rcv
    , salience     = salience_df$salience
    , mz           = mz_lookup(salience_df$feature)
    , rt           = rt_lookup(salience_df$feature)
    )
    my_df[order(-my_df$salience),]
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
        , x_algorithm   = algoC
        , x_prefix      = if (pairSigFeatOnly) {
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
              if ( pairSigFeatOnly ) fully_significant else overall_significant
            ]
            my_matrix <- tdm[ chosen_samples, col_selector, drop = FALSE ]
            my_cor_cov <- do_detail_plot(
              x_dataMatrix  = my_matrix
            , x_predictor   = predictor
            , x_is_match    = is_match
            , x_algorithm   = algoC
            , x_prefix      = if (pairSigFeatOnly) {
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
  } else { # tesC == "none"
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
              # strF(completed)
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
                , x_algorithm   = algoC
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
              , x_algorithm   = algoC
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
        "bad parameter!  sampleMetadata must have at least two levels of factor '%s' matching '%s'"
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
) {
  tryCatch({
    return(
      cor_vs_cov_try( matrix_x, ropls_x, predictor_projection_x, x_progress)
    )
  }, error = function(e) {
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
  matrix_x
, ropls_x
, predictor_projection_x = TRUE
, x_progress = print
) {
  x_class <- class(ropls_x)
  if ( !( as.character(x_class) == "opls" ) ) { # || !( attr(class(x_class),"package") == "ropls" ) )
    stop(
      paste(
        "cor_vs_cov: Expected ropls_x to be of class ropls::opls but instead it was of class "
      , as.character(x_class)
      )
    )
  }
  result <- list()
  result$projection <- projection <- if (predictor_projection_x) 1 else 2
  # suppLs$algoC - Character: algorithm used - "svd" for singular value decomposition; "nipals" for NIPALS
  if ( ropls_x@suppLs$algoC == "nipals") {
    # Equations (1) and (2) from *Supplement to* Wiklund 2008, doi:10.1021/ac0713510
    mag <- function(one_dimensional) sqrt(sum(one_dimensional * one_dimensional))
    mag_xi <- sapply(X = 1:ncol(matrix_x), FUN = function(x) mag(matrix_x[,x]))
    if (predictor_projection_x)
       score_matrix <- ropls_x@scoreMN
    else
       score_matrix <- ropls_x@orthoScoreMN
    score_matrix_transposed <- t(score_matrix)
    score_matrix_magnitude <- mag(score_matrix)
    result$covariance <-
      score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * score_matrix_magnitude )
    result$correlation <-
      score_matrix_transposed %*% matrix_x / ( score_matrix_magnitude * mag_xi )
  } else {
    # WARNING - untested code - I don't have test data to exercise this branch
    # Equations (1) and (2) from Wiklund 2008, doi:10.1021/ac0713510
    # scoreMN - Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
    if (predictor_projection_x)
       score_matrix <- ropls_x@scoreMN
    else
       score_matrix <- ropls_x@orthoScoreMN
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
  result$correlation <- result$correlation[ 1, , drop = TRUE ]
  result$covariance  <- result$covariance [ 1, , drop = TRUE ]

  # Variant 4 of Variable Influence on Projection for OPLS from Galindo_Prieto_2014
  #    Length = number of features; labels = feature identifiers.  (The same is true for $correlation and $covariance.)
  result$vip4p     <- as.numeric(ropls_x@vipVn)
  result$vip4o     <- as.numeric(ropls_x@orthoVipVn)
  result$loadp     <- as.numeric(ropls_x@loadingMN)
  result$loado     <- as.numeric(ropls_x@orthoLoadingMN)
  # get the level names
  level_names      <- sort(levels(as.factor(ropls_x@suppLs$y)))
  fctr_lvl_1       <- level_names[1]
  fctr_lvl_2       <- level_names[2]
  feature_count    <- length(ropls_x@vipVn)
  result$level1    <- rep.int(x = fctr_lvl_1, times = feature_count)
  result$level2    <- rep.int(x = fctr_lvl_2, times = feature_count)
  superresult <- list()
  if (length(result$vip4o) == 0) result$vip4o <- NA
  greaterLevel <- sapply(
    X = result$correlation
  , FUN = function(my_corr)
      tryCatch({
          if ( is.nan( my_corr ) ) {
            print("my_corr is NaN")
            NA 
          } else {
            if ( my_corr < 0 ) fctr_lvl_1 else fctr_lvl_2
          }
        }, error = function(e) {
          x_progress(
            sprintf(
              "cor_vs_cov -> sapply:  error - substituting NA - %s"
            , as.character(e)
            )
          )
          NA
        })
  )

  # begin fixes for https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/issues/1
  featureID          <- names(ropls_x@vipVn)
  greaterLevel       <- greaterLevel[featureID]
  result$correlation <- result$correlation[featureID]
  result$covariance  <- result$covariance[featureID]
  # end fixes for https://github.com/HegemanLab/w4mcorcov_galaxy_wrapper/issues/1

  tsv1 <- data.frame(
    featureID           = featureID
  , factorLevel1        = result$level1
  , factorLevel2        = result$level2
  , greaterLevel        = greaterLevel
  , projection          = result$projection
  , correlation         = result$correlation
  , covariance          = result$covariance
  , vip4p               = result$vip4p
  , vip4o               = result$vip4o
  , loadp               = result$loadp
  , loado               = result$loado
  , row.names           = NULL
  )
  tsv1 <- tsv1[!is.na(tsv1$correlation),]
  tsv1 <- tsv1[!is.na(tsv1$covariance),]
  superresult$tsv1 <- tsv1
  rownames(superresult$tsv1) <- tsv1$featureID
  superresult$projection <- result$projection
  superresult$covariance <- result$covariance
  superresult$correlation <- result$correlation
  superresult$vip4p <- result$vip4p
  superresult$vip4o <- result$vip4o
  superresult$loadp <- result$loadp
  superresult$loado <- result$loado
  superresult$details <- result
  result$superresult <- superresult
  # Include thise in case future consumers of this routine want to use it in currently unanticipated ways
  result$oplsda    <- ropls_x
  result$predictor <- ropls_x@suppLs$y   # in case future consumers of this routine want to use it in currently unanticipated ways
  return (superresult)
}

# vim: sw=2 ts=2 et :
