# read_data_frame - read a w4m data frame, with error handling
#   e.g., data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
read_data_frame <- function(file_path, kind_string, failure_action = failure_action) {
  my.env <- new.env()
  my.env$success <- FALSE
  my.env$msg <- sprintf("no message reading %s", kind_string)
  tryCatch(
    expr = {
      my.env$data    <- utils::read.delim( fill = FALSE, file = file_path )
      my.env$success <- TRUE
    }
  , error = function(e) {
     my.env$ msg <- sprintf("%s read failed", kind_string)
    }
  )
  if (!my.env$success) {
    failure_action(my.env$msg)
    return ( FALSE )
  }
  return (my.env)
}

# read one of three XCMS data elements: dataMatrix, sampleMetadata, variableMetadata
# returns respectively: matrix, data.frame, data.frame, or FALSE if there is a failure
read_xcms_data_element <- function(xcms_data_in, xcms_data_type, failure_action = stop) {
  # note that 'stop' effectively means 'throw'; if 'warning' and 'message' are caught, they mean 'throw' as well
  my_failure_action <- function(...) { failure_action("read_xcms_data_element: ", ...) }
  # xcms_data_type must be in c("sampleMetadata", "variableMetadata", "dataMatrix")
  if ( ! is.character(xcms_data_type) ) {
    my_failure_action(sprintf("bad parameter xcms_data_type '%s'", deparse(xcms_data_type)))
    return ( FALSE )
  }
  if ( 1 != length(xcms_data_type)
       || ! ( xcms_data_type %in% c("sampleMetadata", "variableMetadata", "dataMatrix") ) 
  ) {
    my_failure_action( sprintf("bad parameter xcms_data_type '%s'", xcms_data_type) )
    return ( FALSE )
  }
  if ( is.character(xcms_data_in) ) {
    # case: xcms_data_in is a path to a file
    xcms_data_input_env <- read_data_frame( xcms_data_in, sprintf("%s input", xcms_data_type) )
    if (!xcms_data_input_env$success) {
      my_failure_action(xcms_data_input_env$msg)
      return ( FALSE )
    }
    return ( xcms_data_input_env$data )
  } else {
    # case: xcms_data_in is invalid
    my_failure_action( sprintf("xcms_data_in has unexpected type %s", typeof(xcms_data_in)) )
    return ( FALSE )
  }
}

read_inputs <- function(input_env, failure_action = print) {
  if ( ! is.environment(input_env) ) {
    failure_action("read_inputs: fatal error - 'input_env' is not an environment")
    return ( FALSE )
  }

  if (!is.null(sampleMetadata_in <- input_env$sampleMetadata_in)) {
    # ---
    # read in the sample metadata
    read_data_result <- tryCatchFunc(
      expr = {
        read_xcms_data_element(xcms_data_in = sampleMetadata_in, xcms_data_type = "sampleMetadata")
      }
    )
    if ( read_data_result$success ) {
      smpl_metadata <- read_data_result$value
    } else {
      failure_action(read_data_result$msg)
      return ( FALSE )
    }

    # extract rownames
    rownames(smpl_metadata) <- smpl_metadata[,1]
    
    input_env$smpl_metadata <- smpl_metadata
    # ...
  } else {
    failure_action("read_inputs: fatal error - 'sampleMetadata_in' is missing from 'input_env'")
    return ( FALSE )
  }

  if (!is.null(variableMetadata_in <- input_env$variableMetadata_in)) {
    # ---
    # read in the variable metadata
    read_data_result <- tryCatchFunc(
      expr = {
        read_xcms_data_element(xcms_data_in = variableMetadata_in, xcms_data_type = "variableMetadata")
      }
    )
    if ( read_data_result$success ) {
      vrbl_metadata <- read_data_result$value
    } else {
      failure_action(read_data_result$msg)
      return (FALSE)
    }
    

    # extract rownames (using make.names to handle degenerate feature names)
    err.env <- new.env()
    err.env$success <- FALSE
    err.env$msg <- "no message setting vrbl_metadata rownames"
    tryCatch(
      expr = {
        rownames(vrbl_metadata) <- make.names( vrbl_metadata[,1], unique = TRUE )
        vrbl_metadata[,1] <- rownames(vrbl_metadata)
        err.env$success     <- TRUE
      }
    , error = function(e) {
       err.env$ msg <- sprintf("failed to set rownames for vrbl_metadata read because '%s'", e$message) 
      }
    )
    if (!err.env$success) {
      failure_action(err.env$msg)
      return ( FALSE )
    }

    input_env$vrbl_metadata <- vrbl_metadata
    # ...
  } else {
    failure_action("read_inputs: fatal error - 'variableMetadata_in' is missing from 'input_env'")
    return ( FALSE )
  }

  if (!is.null(dataMatrix_in <- input_env$dataMatrix_in)) {
    # ---
    # read in the data matrix
    read_data_result <- tryCatchFunc(
      expr = {
        read_xcms_data_element(xcms_data_in = dataMatrix_in, xcms_data_type = "dataMatrix")
      }
    )
    if ( read_data_result$success ) {
      data_matrix <- read_data_result$value
    } else {
      failure_action(read_data_result$msg)
      return (FALSE)
    }

    if ( ! is.matrix(data_matrix) ) {
      # extract rownames (using make.names to handle degenerate feature names)
      err.env <- new.env()
      err.env$success <- FALSE
      err.env$msg <- "no message setting data_matrix rownames"
      tryCatch(
        expr = {
          rownames(data_matrix) <- make.names( data_matrix[,1], unique = TRUE )
          err.env$success     <- TRUE
        }
      , error = function(e) {
         err.env$msg <- sprintf("failed to set rownames for data_matrix read because '%s'", e$message) 
        }
      )
      if (!err.env$success) {
        failure_action(err.env$msg)
        return ( FALSE )
      }

      # remove rownames column
      data_matrix <- data_matrix[,2:ncol(data_matrix)]

      # convert data_matrix to matrix from data.frame
      data_matrix <- as.matrix(data_matrix)
    }

    input_env$data_matrix <- data_matrix
    # ...
  } else {
    failure_action("read_inputs: fatal error - 'dataMatrix_in' is missing from 'input_env'")
    return ( FALSE )
  }

  return ( TRUE )
}

