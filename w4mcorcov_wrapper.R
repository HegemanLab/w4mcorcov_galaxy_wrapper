#!/usr/bin/env Rscript

# This script assumes that it is being executed in a current working directory containing the following files:
#   - w4mcorcov_lib.R
#   - w4mcorcov_input.R
#   - w4mcorcov_calc.R

## constants
##----------

modNamC <- "w4mcorcov" ## module name

topEnvC <- environment()
nl <- "\n"

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## subroutines
##----------

source("w4mcorcov_lib.R")
source("w4mcorcov_util.R")
source("w4mcorcov_input.R")
source("w4mcorcov_salience.R")
source("w4mcorcov_calc.R")
source("w4mcorcov_output.R")

## log file
##---------

my_log <- function(x, ...) { cat(paste(iso8601.znow(), " ", x, ..., nl, sep=""))}
my_fatal <- function(x, ...) { 
  my_log("ERROR: ", x, ...)
  quit(save = "no", status = 11, runLast = TRUE)
}

my_log("Start of the '", modNamC, "' Galaxy module call: ")


########
# MAIN #
########

argVc <- unlist(parseCommandArgs(evaluate=FALSE))
my_env <- new.env()

##------------------------------
## Initializing
##------------------------------

## arguments
##----------

# files

my_env$dataMatrix_in        <- as.character(argVc["dataMatrix_in"])
my_env$sampleMetadata_in    <- as.character(argVc["sampleMetadata_in"])
my_env$variableMetadata_in  <- as.character(argVc["variableMetadata_in"])
my_env$contrast_detail      <- as.character(argVc["contrast_detail"])
my_env$contrast_corcov      <- as.character(argVc["contrast_corcov"])
my_env$contrast_salience    <- as.character(argVc["contrast_salience"])
# print(sprintf("contrast_salience: %s", my_env$contrast_salience))

# other parameters

my_env$tesC               <- as.character(argVc["tesC"])
my_env$facC               <- as.character(argVc["facC"])
my_env$pairSigFeatOnly    <- as.logical(argVc["pairSigFeatOnly"])
my_env$levCSV             <- as.character(argVc["levCSV"])
my_env$matchingC          <- as.character(argVc["matchingC"])
my_env$labelFeatures      <- as.character(argVc["labelFeatures"]) # number of features to label at each extreme of the loadings or 'ALL'
my_env$labelOrthoFeatures <- as.logical(argVc["labelOrthoFeatures"])

label_features <- my_env$labelFeatures
labelfeatures_check <- TRUE
if ( is.na(label_features) ) {
  labelfeatures_check <- FALSE
} else if ( is.null(label_features) ) {
  labelfeatures_check <- FALSE
} else if (label_features != "ALL") {
  if ( is.na(as.numeric(label_features)) )
    labelfeatures_check <- FALSE
  else if ( as.numeric(label_features) < 0 )
    labelfeatures_check <- FALSE
}
if ( !labelfeatures_check ) {
  my_log("invalid argument: labelFeatures")
  print(label_features)
  quit(save = "no", status = 10, runLast = TRUE)
}

tsv_action_factory <- function(file, colnames, append) {
  return (
    function(tsv) {
      write.table(
        x = tsv
      , file = file
      , sep = "\t"
      , quote = FALSE
      , row.names = FALSE
      , col.names = colnames
      , append = append
      )
    }
  )
}

corcov_tsv_colnames <- TRUE
corcov_tsv_append   <- FALSE
corcov_tsv_action <- function(tsv) {
  tsv_action_factory(
    file     = my_env$contrast_corcov
  , colnames = corcov_tsv_colnames
  , append   = corcov_tsv_append
  )(tsv)
  corcov_tsv_colnames <<- FALSE
  corcov_tsv_append   <<- TRUE
}

salience_tsv_colnames <- TRUE
salience_tsv_append   <- FALSE
salience_tsv_action <- function(tsv) {
  tsv_action_factory(
    file     = my_env$contrast_salience
  , colnames = salience_tsv_colnames
  , append   = salience_tsv_append
  )(tsv)
  salience_tsv_colnames <<- FALSE
  salience_tsv_append   <<- TRUE
}

my_log( "--------------------------  Reading input data  --------------------------")

# read_inputs is defined in w4mcorcov_input.R
my_result <- read_inputs(input_env = my_env, failure_action = my_log)

if ( is.logical(my_result) && my_result) {
  my_log( "--------------------------  Beginning data processing  --------------------------")

  # receiver for result of the call to corcov_calc
  my_result <- NULL

  # compute and plot the correlation_vs_covariance details plot
  #   The parameter settings here are generally taken from bioconductor ropls::plot.opls source.
  marVn <- c(4.6, 4.1, 2.6, 1.6)
  old_par <- par(
    font      = 2         # bold font face
  , font.axis = 2         # bold font face for axis
  , font.lab  = 2         # bold font face for x and y labels
  , lwd       = 2         # line-width - interpretation is device spcific
  , mar       = marVn     # margins
  , pch       = 18        # black diamond plot-character, see help for graphics::points
  # , mfrow     = c(2,2)    # two rows by two columns
  , pty       = "s"       # force plots to be square
  )
  plot2pdf(
    file.name = my_env$contrast_detail
  , width  = 8
  , height = 8
  , plot.function = function() {
      # plot layout four plots per page
      layout(matrix(1:4, byrow = TRUE, nrow = 2))
      my_result <<- corcov_calc(
          calc_env            = my_env
        , failure_action      = my_fatal
        , progress_action     = my_log
        , corcov_tsv_action   = corcov_tsv_action
        , salience_tsv_action = salience_tsv_action
        )
    }
  )
  par(old_par)
  
  my_log( "--------------------------  Finished data processing  --------------------------")
}

my_log( "End of the '", modNamC, "' Galaxy module call")

if (is.logical(my_result) && my_result) {
  quit(save = "no", status = 0, runLast = TRUE)
} else {
  my_log("failure :(")
  quit(save = "no", status = 10, runLast = TRUE)
}
