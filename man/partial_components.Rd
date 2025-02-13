% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quoll.R
\name{partial_components}
\alias{partial_components}
\title{Render components of a partial inverse}
\usage{
partial_components(
  partial,
  fmt,
  type = c("inv", "adj", "det"),
  determinant = TRUE
)
}
\arguments{
\item{partial}{a partial inverse object}

\item{fmt}{a list of sprintf formatters}

\item{type}{whether the inverse, adjoint or determinant is generated}

\item{determinant}{should uneccessary factors be retained}
}
\value{
a list of character vectors that define
\itemize{
\item \code{cyc} assignments to cycle weight variables
\item \code{det} assignments to determinant variables
\item \code{el} assignments to elements of the inverse or adjoint matrix (if type!="det")
}
}
\description{
Render components of an partial inverse
}
\details{
This is the workhose function that generates expressions to
evaluate the elements of an inverse or adjoint matrix.

The \code{fmt} argument is a list of functions and strings that
construct syntactic elements of the generated code:
\itemize{
\item \code{asub(i,j)} function to generate the i,j element of the weight matrix
\item \code{rsub(i,j)} function to generate the i,j element of the result matrix
\item \code{csub(k)} function to generate the k-th cycle weight variable
\item \code{dsub(k)} function to generate the k-th determinant variable
\item \code{cassign(l,r)} vectorized function to generate cycle weight assignments (l=r)
\item \code{dassign(l,r)} vectorized function to generate determinant assignments (l=r)
\item \code{rassign(l,r)} vectorized function to generate result assignments (l=r)
\item \code{ratio(n,d)} function to generate the ratio operation (n/d)
\item \code{mul} string for the multiplication operator
}
}
\seealso{
\code{\link[=partial_latex]{partial_latex()}}, \code{\link[=partial_R]{partial_R()}}, \code{\link[=partial_cpp]{partial_cpp()}}
}
