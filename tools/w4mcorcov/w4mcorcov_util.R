# tryCatchFunc wraps an expression that produces a value if it does not stop:
#   tryCatchFunc produces a list
#   On success of expr(), tryCatchFunc produces
#     list(success TRUE, value = expr(), msg = "")
#   On failure of expr(), tryCatchFunc produces
#     list(success = FALSE, value = NA, msg = "the error message")
tryCatchFunc <- function(expr) {
  # format error for logging
  format_error <- function(e) {
    paste(c("Error { message:", e$message, ", call:", e$call, "}"), collapse = " ")
  }
  retval <- NULL
  tryCatch(
    expr = {
      retval <- ( list( success = TRUE, value = eval(expr = expr), msg = "" ) )
    }
    , error = function(e) {
      retval <<- list( success = FALSE, value = NA, msg = format_error(e) )
    }
  )
  return (retval)
}

errorSink <- function(which_function, ...) {
  var_args <- "..."
  tryCatch(
    var_args <<- (deparse(..., width.cutoff = 60))
  , error = function(e) {print(e$message)}
  )
  if (var_args == "...")
    return
  # format error for logging
  format_error <- function(e) {
    sprintf(
      "Error\n{  message: %s\n, arguments: %s\n}\n"
    , e$message
    , Reduce(f = paste, x = var_args)
    )
  }
  format_warning <- function(e) {
    sprintf(
      "Warning\n{  message: %s\n, arguments: %s\n}\n"
    , e$message
    , Reduce(f = paste, x = var_args)
    )
  }
  sink_number <- sink.number()
  sink(stderr())
  tryCatch(
    var_args <- (deparse(..., width.cutoff = 60))
  , expr = {
      retval <- which_function(...)
    }
    , error = function(e) cat(format_error(e), file = stderr())
    , warning = function(w) cat(format_warning(w), file = stderr())
  )
  while (sink.number() > sink_number) {
    sink()
  }
}
errorPrint <- function(...) {
  errorSink(which_function = print, ...)
}
errorCat <- function(...) {
  errorSink(which_function = cat, ..., "\n")
}


# # pseudo-inverse - computational inverse of non-square matrix a
# p.i <- function(a) {
#   solve(t(a) %*% a) %*% t(a)
# } 

# vim: sw=2 ts=2 et ai :
