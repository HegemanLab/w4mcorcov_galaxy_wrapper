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


# # pseudo-inverse - computational inverse of non-square matrix a
# p.i <- function(a) {
#   solve(t(a) %*% a) %*% t(a)
# } 


