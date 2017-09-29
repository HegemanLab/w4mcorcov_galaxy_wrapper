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

# turn off all plotting devices
dev.off.all <- function() {
  while (!is.null(dev.list())) { dev.off() }
}
  
# capture plot and write to PDF; then close any devices opened in the process
plot2pdf <- function(
  file.name
, plot.function
, width = 12
, height = 12
) {
  # capture plot and write to PDF
  cur.dev <- dev.list()
  filename <- file.name
  pdf(file = filename, width = width, height = height)
  plot.function()
  # close any devices opened in the process
  dev.off()
  if (is.null(cur.dev)) {
      dev.off.all()
  } else {
      while ( length(dev.list()) > length(cur.dev) ) { dev.off() }
  }
}

# print and capture plot and write to PDF; then close any devices opened in the process
#   This is needed for ggplot which does not print the plot when invoked within a function.
print2pdf <- function(
  file.name
, plot.function
, width = 12
, height = 12
) {
  plot2pdf(
    file.name = file.name
  , width = width
  , height = height
  , plot.function = function() {
      print(plot.function())
    }
  )
}

# iso8601.filename.fragment <- function()
# {
#   strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y%m%d-%H%M%S")
# }
#
# pdf.name <- function(name)
# {
#   paste0(name, "_", iso8601.filename.fragment(), ".pdf")
# }
#
# tsv.name <- function(name)
# {
#   paste0(name, "_", iso8601.filename.fragment(), ".tsv")
# }
#
# # pseudo-inverse - computational inverse non-square matrix a
# p.i <- function(a) {
#   solve(t(a) %*% a) %*% t(a)
# } 


