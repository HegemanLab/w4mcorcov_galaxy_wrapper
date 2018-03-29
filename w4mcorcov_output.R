
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

iso8601.znow <- function()
{
  strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

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

