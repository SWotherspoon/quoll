% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quoll.R
\name{graph_diagram}
\alias{graph_diagram}
\title{Graph diagram}
\usage{
graph_diagram(
  A,
  rankdir = c("BT", "TB", "LR", "RL"),
  vertex.col = colorspec_vector(A),
  edge.col = colorspec_matrix(A)
)
}
\arguments{
\item{A}{adjacency matrix of the undirected graph}

\item{rankdir}{the layout direction of the graph}

\item{vertex.col}{a colorspec vector of vertex colors}

\item{edge.col}{a colorspec matrix of edge colors}
}
\value{
a dot description of the graph
}
\description{
Create a diagram of an undirected graph
}
\details{
Creates a dot description of a graph, with specified vertices
highlighted in orange, to be rendered with \code{\link[DiagrammeR:grViz]{DiagrammeR::grViz()}}.
The adjacency matrix is assumed to be symmetric.
}
\examples{
# Adjacency matrix of a directed graph
A <- matrix(c(1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,1),4,4)
graph_diagram(A)
}
