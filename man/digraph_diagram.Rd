% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quoll.R
\name{digraph_diagram}
\alias{digraph_diagram}
\title{Directed graph diagram}
\usage{
digraph_diagram(
  A,
  rankdir = c("BT", "TB", "LR", "RL"),
  vertex.col = colorspec_vector(A),
  edge.col = colorspec_matrix(A),
  hfrac = 0.99,
  split.edges = FALSE
)
}
\arguments{
\item{A}{adjacency matrix of the graph}

\item{rankdir}{the layout direction of the graph}

\item{vertex.col}{a colorspec vector of vertex colors}

\item{edge.col}{a colorspec matrix of edge colors}

\item{hfrac}{the fraction of a bidirectional edge to highlight}

\item{split.edges}{should bidirectional edges be drawn as two
unidirectional edges?}
}
\value{
a dot description of the graph
}
\description{
Create a diagram of a directed graph
}
\details{
Creates a dot description of a directed graph, with specified
vertices highlighted in orange, to be rendered with
\code{\link[DiagrammeR:grViz]{DiagrammeR::grViz()}}.
}
\examples{
# Adjacency matrix of a directed graph
A <- matrix(c(1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1),4,4)
digraph_diagram(A)
}
\seealso{
\code{\link[=graph_diagram]{graph_diagram()}} for undirected graphs

\code{\link[=community_diagram]{community_diagram()}} for signed community matrices
}
