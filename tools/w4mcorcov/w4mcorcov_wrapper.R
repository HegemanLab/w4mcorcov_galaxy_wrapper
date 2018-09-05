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

# from: https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}

script.dir <-  LocationOfThisScript()

source(paste(script.dir, "w4mcorcov_lib.R", sep="/")) 
source(paste(script.dir, "w4mcorcov_util.R", sep="/")) 
source(paste(script.dir, "w4mcorcov_input.R", sep="/")) 
source(paste(script.dir, "w4mcorcov_salience.R", sep="/")) 
source(paste(script.dir, "w4mcorcov_calc.R", sep="/")) 
source(paste(script.dir, "w4mcorcov_output.R", sep="/")) 

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

errorPrint(sessionInfo())

argVc <- unlist(parseCommandArgs(evaluate=FALSE))
errorCat("\n\n---\n\nArguments that were passed to R are as follows:\n")
errorPrint(argVc)

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

# other parameters

my_env$tesC               <- as.character(argVc["tesC"])
my_env$facC               <- as.character(argVc["facC"])
my_env$pairSigFeatOnly    <- as.logical(argVc["pairSigFeatOnly"])
my_env$levCSV             <- as.character(argVc["levCSV"])
my_env$matchingC          <- as.character(argVc["matchingC"])
my_env$labelFeatures      <- as.character(argVc["labelFeatures"]) # number of features to label at each extreme of the loadings or 'ALL'
my_env$cplot_o            <- as.logical(argVc["cplot_o"]) # TRUE if orthogonal C-plot is requested
my_env$cplot_p            <- as.logical(argVc["cplot_p"]) # TRUE if parallel C-plot is requested
my_env$cplot_y            <- as.character(argVc["cplot_y"]) # Choice of covariance/correlation for Y-axis on C-plot

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
  if ( my_env$cplot_p || my_env$cplot_o ) {
    old_par <- par(
      font        = 2         # bold font face
    , font.axis   = 2         # bold font face for axis
    , font.lab    = 2         # bold font face for x and y labels
    , lwd         = 2         # line-width - interpretation is device spcific
    , pch         = 18        # black diamond plot-character, see help for graphics::points
    , pty         = "m"       # do not force plots to be square
    , no.readonly = TRUE      # only save writable parameters
    )
    pdf_height <- 12
    pdf_width  <- 8
    my_layout <- function() {
      # lay out 2 columns by 3 rows with extra width at the margin of individual plots
      layout(
        matrix(
          # blank row  plot 1 & 2  blank row  plot 3 & 4  blank row  plot 5 & 6 blank row
          c(0,0,0,0,0, 0,1,0,2,0,  0,0,0,0,0, 0,3,0,4,0,  0,0,0,0,0, 0,5,0,6,0, 0,0,0,0,0)
        , nrow = 7
        , ncol = 5
        , byrow = TRUE
        )
        # slim columns 1, 3, and 5
      , widths  = c(0.1, 0.9, 0.1, 0.9, 0.1)
        # slim rows 1, 3, 5, and 7
      , heights = c(0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.1)
      )
    }
  } else {
    old_par <- par(
      font        = 2         # bold font face
    , font.axis   = 2         # bold font face for axis
    , font.lab    = 2         # bold font face for x and y labels
    , lwd         = 2         # line-width - interpretation is device spcific
    , pch         = 18        # black diamond plot-character, see help for graphics::points
    , pty         = "m"       # do not force plots to be square
    , no.readonly = TRUE      # only save writable parameters
    )
    pdf_height <- 8
    pdf_width  <- 8
    my_layout <- function() {
      # lay out 2 columns by 2 rows with extra width at the margin of individual plots
      layout(
        matrix(
          # blank row  plot 1 & 2  blank row  plot 3 & 4  blank row
          c(0,0,0,0,0, 0,1,0,2,0,  0,0,0,0,0, 0,3,0,4,0,  0,0,0,0,0)
        , nrow = 5
        , ncol = 5
        , byrow = TRUE
        )
        # slim columns 1, 3, and 5
      , widths  = c(0.1, 0.9, 0.1, 0.9, 0.1)
        # slim rows 1, 3, and 5
      , heights = c(0.1, 0.9, 0.1, 0.9, 0.1)
      )
    }
  }
  plot2pdf(
    file.name = my_env$contrast_detail
  , width  = pdf_width
  , height = pdf_height
  , plot.function = function() {
      # plot layout four or six plots per page
      my_layout()
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
