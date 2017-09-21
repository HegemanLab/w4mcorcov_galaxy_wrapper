#!/usr/bin/env Rscript

# This script assumes that it is being executed in a current working directory containing the following files:
#   - w4mcorcov_lib.R
#   - w4mcorcov_input.R
#   - w4mcorcov_calc.R

## constants
##----------

modNamC <- "w4mcorcov" ## module name

topEnvC <- environment()
flgC <- "\n"

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## log file
##---------

my_print <- function(x, ...) { cat(paste(x, ..., sep=""))}

my_print("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), flgC)

## subroutines
##----------

source("w4mcorcov_lib.R")
source("w4mcorcov_util.R")
source("w4mcorcov_input.R")
source("w4mcorcov_calc.R")
source("w4mcorcov_output.R")


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
my_env$contrast_grid        <- as.character(argVc["contrast_grid"])
my_env$contrast_detail      <- as.character(argVc["contrast_detail"])
my_env$contrast_corcov      <- as.character(argVc["contrast_corcov"])
my_env$variableMetadata_out <- as.character(argVc["variableMetadata_out"])

# other parameters

my_env$tesC            <- as.character(argVc["tesC"])
my_env$facC            <- as.character(argVc["facC"])
my_env$pairSigFeatOnly <- as.character(argVc["pairSigFeatOnly"])
my_env$levCSV          <- as.character(argVc["levCSV"])
my_env$matchingC       <- as.character(argVc["matchingC"])

# read_inputs is defined in w4mcorcov_input.R
my_result <- read_inputs(input_env = my_env, failure_action = my_print)

if ( is.logical(my_result) && my_result) {
  print("Data Matrix")
  ropls::strF(my_env$data_matrix)
  print("Variable Metadata")
  ropls::strF(my_env$vrbl_metadata)
  print("Sample Metadata")
  ropls::strF(my_env$smpl_metadata)

  my_result <- corcov_calc(calc_env = my_env, failure_action = my_print)
}

if (is.logical(my_result) && my_result) {
  cat("success :)\n")
} else {
  cat("failure :(\n")
}
