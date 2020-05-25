#' Call the fastp binary with additional arguments.
#'
#' Adapted from the Rbowtie package
#'
#' @keywords internal
#'
#' @param bin The name of the fastp binary
#' @param execute Logical scalar, whether to execute the command. If FALSE,
#'   return a string with the shell command.
#' @param args A character string containing the arguments that will be passed
#'   to the binary.
#'
#' @return If \code{execute} is TRUE, returns the console output of running the
#'   fastp command. If \code{execute} is FALSE, returns the shell command.
#'
.fastpBin <- function(bin="fastp", args="", execute=TRUE) {
  if (is.null(args) || args=="") {
    stop("The fastp binary needs to be called with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  bin <- match.arg(bin)
  call <- paste(shQuote(file.path(system.file(package="Rfastp"), bin)), args)
  if (!execute) {
    return(call)
  }
  output <- system(call, intern=TRUE)
  return(output)
}


#' Call the gencore binary with additional arguments.
#'
#' Adapted from the gencore package
#'
#' @keywords internal
#'
#' @param bin The name of the gencore binary
#' @param execute Logical scalar, whether to execute the command. If FALSE,
#'   return a string with the shell command.
#' @param args A character string containing the arguments that will be passed
#'   to the binary.
#'
#' @return If \code{execute} is TRUE, returns the console output of running the
#'   gencore command. If \code{execute} is FALSE, returns the shell command.
#'
.gencoreBin <- function(bin="gencore", args="", execute=TRUE) {
  if (is.null(args) || args=="") {
    stop("The fastp binary needs to be called with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  bin <- match.arg(bin)
  call <- paste(shQuote(file.path(system.file(package="Rfastp"), bin)), args)
  if (!execute) {
    return(call)
  }
  output <- system(call, intern=TRUE)
  return(output)
}


