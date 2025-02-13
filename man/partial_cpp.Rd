% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quoll.R
\name{partial_cpp}
\alias{partial_cpp}
\title{Render partial object as R code}
\usage{
partial_cpp(partial, type = c("inv", "adj"), determinant = TRUE, fname = type)
}
\arguments{
\item{partial}{a partial inverse object}

\item{type}{whether the inverse or adjoint is generated}

\item{determinant}{should uneccessary factors of the determinant be retained}

\item{fname}{the function name}
}
\value{
a character string of R code.
}
\description{
Render a \code{partial} object as Rcpp code to calculate the partial
inverse or adjoint matrix.
}
\details{
This function translates a \code{partial} object generated by
\code{\link[=partial]{partial()}} to Rcpp code to evaluate the elements of the inverse
or adjoint matrix.  In the generated code,
\itemize{
\item the argument \code{A} is a matrix of (unsigned) edge weights,
\item the \code{C} variables are cycle weights,
\item the \code{D} variables are determinants of strongly connected components.
}

If \code{determinant=FALSE}, any \code{D} variables not required for the
subsequent calculation areomitted.  If \code{determinant=TRUE}, no
variables are omitted and the product of the \code{D} variables is the
determinant.
}
\seealso{
\code{\link[=partial]{partial()}}
}
