##' Create a diagram of an undirected graph
##'
##' Creates a dot description of a graph, with specified vertices
##' highlighted in orange, to be rendered with [DiagrammeR::grViz()].
##' The adjacency matrix is assumed to be symmetric.
##'
##' @title Graph diagram
##' @param A adjacency matrix of the undirected graph
##' @param rankdir the layout direction of the graph
##' @param vertex.col a colorspec vector of vertex colors
##' @param edge.col a colorspec matrix of edge colors
##' @return a dot description of the graph
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,1),4,4)
##' graph_diagram(A)
##' @export
graph_diagram <- function(A,rankdir=c("BT","TB","LR","RL"),
                          vertex.col=colorspec_vector(A),
                          edge.col=colorspec_matrix(A)) {

  rankdir <- match.arg(rankdir)
  nms <- if(is.null(rownames(A))) seq_len(nrow(A)) else rownames(A)

  pal <- attr(vertex.col,"palette")
  col <- pal[vertex.col]
  str <- paste(
    "graph G {",
    sprintf("  graph [rankdir=%s,overlap=true,fontsize=8]",rankdir),
    "  node [shape=circle,style=filled]",
    paste(sprintf("  %d [label='%s',color='%s',fillcolor='%s']",seq_len(nrow(A)),nms,col,col),collapse="\n"),
    "  edge []",
    sep="\n")

  pal <- attr(edge.col,"palette")

  ## Non loop edges
  ij <- which(A!=0L,arr.ind=TRUE)
  ij <- ij[ij[,1L] < ij[,2L],,drop=FALSE]
  if(nrow(ij)>0L) {
    col <- pal[edge.col[ij]]
    str <- paste(str,paste(
      sprintf("  %d -- %d [color='%s']",ij[,1],ij[,2],col),
      collapse="\n"),
      sep="\n")
  }

  ## Self loops
  ii <- which(diag(A)!=0L)
  if(length(ii)>0L) {
    col <- pal[edge.col[cbind(ii,ii)]]
    str <- paste(str,paste(
      sprintf("  %d -- %d [color='%s']",ii,ii,col),
      collapse="\n"),
      sep="\n")
  }

  str <- paste(str,"}",sep="\n")
  str
}


##' Create a diagram of a directed graph
##'
##' Creates a dot description of a directed graph, with specified
##' vertices highlighted in orange, to be rendered with
##' [DiagrammeR::grViz()].
##'
##' @title Directed graph diagram
##' @param A adjacency matrix of the graph
##' @param rankdir the layout direction of the graph
##' @param vertex.col a colorspec vector of vertex colors
##' @param edge.col a colorspec matrix of edge colors
##' @param hfrac the fraction of a bidirectional edge to highlight
##' @param split.edges should bidirectional edges be drawn as two
##'   unidirectional edges?
##' @return a dot description of the graph
##' @seealso [graph_diagram()] for undirected graphs
##' @seealso [community_diagram()] for signed community matrices
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1),4,4)
##' digraph_diagram(A)
##' @export
digraph_diagram <- function(A,rankdir=c("BT","TB","LR","RL"),
                            vertex.col=colorspec_vector(A),
                            edge.col=colorspec_matrix(A),
                            hfrac=0.99,split.edges=FALSE) {

  rankdir <- match.arg(rankdir)
  nms <- if(is.null(rownames(A))) seq_len(nrow(A)) else rownames(A)

  pal <- attr(vertex.col,"palette")
  col <- pal[vertex.col]
  str <- paste(
    "digraph G {",
    sprintf("  graph [rankdir=%s,overlap=true,fontsize=8]",rankdir),
    "  node [shape=circle,style=filled]",
    paste(sprintf("  %d [label='%s',color='%s',fillcolor='%s']",seq_len(nrow(A)),nms,col,col),collapse="\n"),
    "  edge [arrowsize=0.5]",
    sep="\n")

  pal <- attr(edge.col,"palette")

  ## Bidirectional edges
  if(!split.edges) {
    ij <- which((A!=0L & t(A)!=0L), arr.ind=TRUE)
    ij <- ij[ij[,1] < ij[,2],,drop=FALSE]
    if(nrow(ij)>0L) {
      colb <- edge.col[ij]
      colf <- edge.col[ij[,2:1,drop=FALSE]]
      split <- ifelse(colb > 2 & colf > 2,0.5,1-hfrac)
      col <- ifelse(colf==colb,pal[colf],ifelse(colf>colb,
                    sprintf("%s:%s;%0.2f",pal[colf],pal[colb],split),
                    sprintf("%s;%0.2f:%s",pal[colf],split,pal[colb])))
      str <- paste(str,paste(
        sprintf("  %d -> %d [dir=both,color='%s']",ij[,1],ij[,2],col),
        collapse="\n"),
        sep="\n")
    }
  }

  if(!split.edges) {
    ## Unidirectional (non loop) edges
    ij <- which((A!=0L & t(A)==0L), arr.ind=TRUE)
  } else {
    ## Non loop edges
    ij <- which(A!=0L, arr.ind=TRUE)
    ij <- ij[ij[,1]!=ij[,2],,drop=FALSE]
  }
  if(nrow(ij)>0L) {
    col <- pal[edge.col[ij]]
    str <- paste(str,paste(
      sprintf("  %d -> %d [dir=forward,color='%s']",ij[,1],ij[,2],col),
      collapse="\n"),
      sep="\n")
  }

  ## Self loops
  ii <- which(diag(A)!=0L)
  if(length(ii)>0L) {
    col <- pal[edge.col[cbind(ii,ii)]]
    str <- paste(str,paste(
      sprintf("  %d -> %d [dir=forward,color='%s']",ii,ii,col),
      collapse="\n"),
      sep="\n")
  }

  str <- paste(str,"}",sep="\n")
  str
}

##' Create a diagram of the directed graph representation of a
##' network model.
##'
##' Creates a dot description of the directed graph representation of
##' a network model, to be rendered with [DiagrammeR::grViz()].
##'
##' @title Network diagram
##' @param S a signed community matrix
##' @param rankdir the layout direction of the graph
##' @param vertex.col a colorspec vector of vertex colors
##' @param edge.col a colorspec matrix of edge colors
##' @param hfrac the fraction of a bidirectional edge to highlight
##' @param split.edges should bidirectional edges be drawn as two
##'   unidirectional edges?
##' @return a dot description of the graph
##' @export
community_diagram <- function(S,rankdir=c("BT","TB","LR","RL"),
                              vertex.col=colorspec_vector(S),
                              edge.col=colorspec_matrix(S),
                              hfrac=0.99,split.edges=FALSE) {

  rankdir <- match.arg(rankdir)
  nms <- if(is.null(rownames(S))) seq_len(nrow(S)) else rownames(S)

  pal <- attr(vertex.col,"palette")
  col <- pal[vertex.col]
  str <- paste(
    "digraph G {",
    sprintf("  graph [rankdir=%s,overlap=true,fontsize=8]",rankdir),
    "  node [shape=circle,style=filled]",
    paste(sprintf("  %d [label='%s',color='%s',fillcolor='%s']",seq_len(nrow(S)),nms,col,col),collapse="\n"),
    "  edge [arrowsize=0.5]",
    sep="\n")

  pal <- attr(edge.col,"palette")

  if(!split.edges) {
    ## Bidirectional edges
    ij <- which((S!=0L & t(S)!=0L), arr.ind=TRUE)
    ij <- ij[ij[,1] > ij[,2],,drop=FALSE]
    if(nrow(ij)>0L) {
      colb <- edge.col[ij]
      colf <- edge.col[ij[,2:1,drop=FALSE]]
      split <- ifelse(colb > 2 & colf > 2,0.5,1-hfrac)
      col <- ifelse(colf==colb,pal[colf],ifelse(colf>colb,
                    sprintf("%s:%s;%0.2f",pal[colf],pal[colb],split),
                    sprintf("%s;%0.2f:%s",pal[colf],split,pal[colb])))
      str <- paste(str,paste(
        sprintf("  %d -> %d [dir=both,arrowtail=%s,arrowhead=%s,color='%s']",ij[,2],ij[,1],
                ifelse(S[ij[,2:1,drop=FALSE]]==1L,"vee","dot"),
                ifelse(S[ij]==1L,"vee","dot"),
                col),
        collapse="\n"),
        sep="\n")
    }
  }

  if(!split.edges) {
    ## Unidirectional (non loop) edges
    ij <- which((S!=0L & t(S)==0L), arr.ind=TRUE)
  } else {
    ## Non loop edges
    ij <- which(S!=0L, arr.ind=TRUE)
    ij <- ij[ij[,1]!=ij[,2],,drop=FALSE]
  }
  if(nrow(ij)>0L) {
    col <- pal[edge.col[ij]]
    str <- paste(str,paste(
      sprintf("  %d -> %d [dir=forward,arrowhead=%s,color='%s']",ij[,2],ij[,1],
              ifelse(S[ij]==1L,"vee","dot"),
              col),
      collapse="\n"),
      sep="\n")
  }

  ## Self loops
  ii <- which(diag(S)!=0L)
  if(length(ii)>0L) {
    col <- pal[edge.col[cbind(ii,ii)]]
    str <- paste(str,paste(
      sprintf("  %d -> %d [dir=forward,arrowhead=%s,color='%s']",ii,ii,
                ifelse(S[cbind(ii,ii)]==1L,"vee","dot"),
                col),
      collapse="\n"),
      sep="\n")
  }

  str <- paste(str,"}",sep="\n")
  str
}


##' Specify the colouring of vertices and edges in a graph or
##' digraph diagram.
##'
##' In the graph diagrams created by [graph_diagram()],
##' [digraph_diagram()], and [community_diagram()], vertex colours are
##' specified by a vector of color indices into a palette of vertex
##' colours, and edge colours are specified by a matrix of color
##' indices into a palette of edge colours.
##'
##' The `colorspec_vector` and `colorspec_matrix` functions create
##' vectors and matrices to specify edge and vertex colours.  By
##' default, all vertices are given the first colour in the vertex
##' palette, and all edges are given the first colour in the edge
##' palette.  The user can then modify these as desired.
##'
##' The `vertex_colorspec` and `edge_colorspec` functions create
##' colorspec vectors and matrices from a symbolic representation of
##' vertex and edge colors. For `vertex_colorspec`, a vector of vertex
##' names is given a named argument that determines the colour applied
##' to those vertices.  For `edge_colorspec`, a block of edge
##' specifications is given as a named argument that determines the
##' colour applied to those edges.  The symbolic representation only
##' identifies the (directed) edges to be coloured - it cannot create
##' or change edges.
##'
##' @title Graph color features
##' @param A an adjacency matrix
##' @param palette a vector of colors
##' @param vertex.col the default vertex colour
##' @param edge.col the default edge colour
##' @param loop.col the default loop color
##' @param ... symbolic vertex and edge coloring directives
##' @return a colorspec vector or matrix
##' @examples
##' # Define a directed graph
##' A <- digraph({
##'   vertices(A,B,C,D,E)
##'   A %->% B %->% D
##'   A %->% C %->% D
##'   B %<->% E
##' })
##' # Specify vertex colors
##' g <- digraph_diagram(A,
##'   vertex.col=vertex_colspec(A,orange=c(A,D)),
##'   edge.col=edge_colspec(A,tomato={A %->% B %->% D},
##'                           orangered={A %->% C %->% D}))
##' cat(g)
##' @rdname colorspec
##' @export
colorspec_vector <- function(A,palette=c("aliceblue")) {
  h <- rep(1,nrow(A))
  attr(h,"palette") <- palette
  h
}

##' @rdname colorspec
##' @export
colorspec_matrix <- function(A,palette=c("gray80","gray90")) {
  H <- matrix(1L,nrow(A),ncol(A),dimnames=dimnames(A))
  diag(H) <- 2L
  attr(H,"palette") <- palette
  H
}


##' @rdname colorspec
##' @export
vertex_colspec <- function(A,vertex.col="aliceblue",...) {
  map <- setNames(as.list(seq_len(nrow(A))),rownames(A))
  spec <- substitute(list(...))
  h <- colorspec_vector(A,palette=c(vertex.col,names(spec)[-1L]))
  for(k in seq_len(length(spec)-1L)) {
    h[eval(spec[[k+1]],map)] <- k+1L
  }
  h
}

##' @rdname colorspec
##' @export
edge_colspec <- function(A,edge.col="gray80",loop.col="gray90",...) {

  `%-o%` <- `%->%` <- function(a,b) {
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    H[b,a] <<- col
    b
  }
  
  `%o-%` <-`%<-%` <- function(a,b) {
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    H[a,b] <<- col
    b
  }
  
  `%<->%` <- `%o-o%` <- `%o->%` <- `%<-o%` <- function(a,b) {
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    H[b,a] <<- col
    H[a,b] <<- col
    b
  }
  
  map <- setNames(as.list(seq_len(nrow(A))),rownames(A))
  spec <- substitute(list(...))
  H <- colorspec_matrix(A,palette=c(edge.col,loop.col,names(spec)[-1L]))
  for(k in seq_len(length(spec)-1L)) {
    col <- k+2L
    eval(spec[[k+1]])
  }
  H
}


##' Symbolically specify a directed graph.
##'
##'
##' Generate an adjacency matrix from a sequence of symbolic
##' directives that describe a directed graph.  The directives are:
##'
##' * `vertices(v1,v2,...)` defines the vertices of the graph.
##' * `a %->% b` defines a unidirectional edge from `a` to `b`
##' * `a %<-% b` defines a unidirectional edge from `b` to `a`
##' * `a %<->% b` is equivalent to `a %<-% b` and `a %->% b`
##'
##' @title Define an adjacency matrix
##' @param defn a sequence of statements describing the digraph
##' @return an adjacency matrix
##' @importFrom stats setNames
##' @examples
##' # Define a directed graph
##' A <- digraph({
##'   vertices(A,B,C,D,E)
##'   A %->% B %->% D
##'   A %->% C %->% D
##'   B %<-% E
##' })
##' A
##' @export
digraph <- function(defn) {

  vertices  <- function(...) {
    if(!is.null(map)) stop("vertices already defined")
    nms <- as.character(match.call()[-1])
    n <- length(nms)
    map <<- setNames(as.list(seq_len(n)),nms)
    A <<- matrix(0L,n,n,dimnames=list(nms,nms))
  }
  
  `%->%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    A[b,a] <<- 1L
    b
  }
  
  `%<-%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    A[a,b] <<- 1L
    b
  }
  
  `%<->%` <-function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    A[b,a] <<- 1L
    A[a,b] <<- 1L
    b
  }
  
  
  A <- 0L
  map <- NULL
  eval(substitute(defn))
  A
}

##' Symbolically specify a directed graph representing a network
##' model.
##'
##'
##' Generate a signed adjacency matrix from a sequence of symbolic
##' directives that describe the network model.  The directives are:
##'
##' * `vertices(v1,v2,...)` defines the vertices of the graph.
##' * `loop(v1,v2,...)` defines self loops for the specified
##'   vertices.
##' * `a %->% b` defines a unidirectional edge where `a` has a
##'   positive impact on `b`
##' * `a %-o% b` defines a unidirectional edge where `a` has a
##'   negative impact on `b`
##' * `a %<-% b` defines a unidirectional edge where `b` has a
##'   positive impact on `a`
##' * `a %o-% b` defines a unidirectional edge where `b` has a
##'   negative impact on `a`
##' * `a %o->% b` is equivalent to `a %0-% b` and `a %->% b`
##' * `a %<-o% b` is equivalent to `a %<-% b` and `a %-o% b`
##' * `a %<->% b` is equivalent to `a %<-% b` and `a %->% b`
##' * `a %<->% b` is equivalent to `a %o-% b` and `a %-o% b`
##'
##'
##' @title Define a signed community matrix
##' @param defn a sequence of statements describing the network model
##' @return a signed adjacency matrix
##' @references Dambacher, J. M., Li, H. W., & Rossignol,
##'   P. A. (2002).  Relevance of community structure in assessing
##'   indeterminacy of ecological predictions. *Ecology*, 83(5),
##'   1372-1385.
##' @importFrom stats setNames
##' @examples
##' # Define a signed community matrix
##' S <- community({
##'   vertices(A,B,C,D,E)
##'   loop(A,B,C,D,E)
##'   A %->% B %->% D
##'   A %->% C %->% D
##'   B %o->% E
##' })
##' S
##' @export
community <- function(defn) {

  vertices  <- function(...) {
    if(!is.null(map)) stop("vertices already defined")
    nms <- as.character(match.call()[-1])
    n <- length(nms)
    map <<- setNames(as.list(seq_len(n)),nms)
    S <<- matrix(0L,n,n,dimnames=list(nms,nms))
  }
  loop <- function(...) {
    if(is.null(map)) eval(substitute(vertices(...)))
    for(v in eval(substitute(c(...)),map))
      S[v,v] <<- -1L
  }
  
  `%->%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- 1L
    b
  }
  
  `%-o%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- -1
    b
  }
  
  `%<-%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[a,b] <<- 1L
    b
  }
  
  `%o-%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[a,b] <<- -1L
    b
  }
  
  `%o->%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- 1L
    S[a,b] <<- -1L
    b
  }
  
  `%<-o%` <- function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- -1L
    S[a,b] <<- 1L
    b
  }
  
  `%<->%` <-function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- 1L
    S[a,b] <<- 1L
    b
    }
  
  `%o-o%` <-function(a,b) {
    if(is.null(map)) stop("no vertices defined")
    a <- eval(substitute(a),map)
    b <- eval(substitute(b),map)
    S[b,a] <<- -1L
    S[a,b] <<- -1L
    b
  }
  
  S <- 0L
  map <- NULL
  eval(substitute(defn))
  S
}




##' Find all paths between two vertices, ignoring self loops
##'
##' @title Paths in a directed graph
##' @param A adjacency matrix of a directed graph
##' @param from initial vertex of path
##' @param to final vertex of path
##' @return a list of paths, represented as a vector of vertices
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1),4,4)
##' # Find all paths from vertex 1 to vertex 4
##' find_paths(A,1,4)
##' @export
find_paths <- function(A,from,to) {

  if(from==to) return(list(from))
  paths <- list()

  dfs <- function(v,p) {
    ## Adjacent unvisited vertices
    vs <- A[,v]
    vs[p] <- 0L
    for(v in which(vs!=0L)) {
      if(v==to)
        paths[[length(paths)+1L]] <<- c(p,to)
      else
        dfs(v,c(p,v))
    }
  }
  dfs(from,from)
  paths
}


##' Find all cycles in a directed graph, optionally including self
##' loops.
##'
##' Returns a list of the cycles of a directed graph, optionally
##' including self loops.  Cycles are represented as vectors of
##' vertices, lowest numbered vertex first.
##'
##' @title Cycles of a directed graph
##' @param A adjacency matrix of a directed graph
##' @param loops should self loops be included
##' @return a list of cycles
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1),4,4)
##' # Find all cycles in the graph, excluding loops
##' directed_cycles(A,loops=FALSE)
##' @export
directed_cycles <- function(A,loops=TRUE) {

  cycles <- list()

  dfs <- function(v0,v,p) {
    ## Adjacent vertices
    vs <- A[,v]
    ## Remove vertices less than v0, visited vertices and loops
    vs[c(seq_len(v0-1L),p,v)] <- 0L
    ## Search unvisited vertices
    for(v in which(vs!=0L)) {
      if(v==v0)
        cycles[[length(cycles)+1L]] <<- c(v0,p)
      else
        dfs(v0,v,c(p,v))
      }
  }
  for(v0 in seq_len(nrow(A))) dfs(v0,v0,integer(0L))
  if(loops) cycles <- c(cycles,lapply(which(unname(diag(A)!=0L)),function(v) v))

  cycles
}


##' Create the intersection graph of a list of sets
##'
##' Creates an adjacency matrix representing the intersection graph of
##' a collection of sets.  In the intersection graph, vertices
##' represent sets and edges connect sets that have a non-empty
##' intersection.
##'
##' Two implementations are provided:
##' * intersection_graph() - works with sets of any type
##' * intersection_graph_int() - optimized for integer sets
##'
##' @title Intersection graph
##' @param sets a list of vectors, each representing a set
##' @param maxel the maximum element of any set
##' @return the adjacency matrix of the intersection graph
##' @examples
##' # List of sets
##' sets <- list(c(1,2,3), c(2,3,4), c(3,4,5))
##' # Create the intersection graph
##' intersection_graph(sets)
##' @export
intersection_graph <- function(sets) {

  A <- matrix(0L,length(sets),length(sets))
  for(i in seq_len(nrow(A))) {
    for(j in seq_len(i)) {
      if(any(sets[[i]] %in% sets[[j]]))
        A[j,i] <- A[i,j] <- 1L
    }
  }
  A
}


##' @rdname intersection_graph
##' @export
intersection_graph_int <- function(sets,maxel) {

  X <- matrix(0L,maxel,length(sets))
  for(i in seq_len(ncol(X))) X[sets[[i]],i] <- 1L
  ifelse(crossprod(X)>0,1,0)
}


##' Find the maximal independent sets of a graph
##'
##' Returns a list of the maximal independent sets of a subgraph of a
##' graph, represented as vectors of vertices. The indendependend sets
##' can be filtered to exclude set containing certain vertices, and to
##' include only sets containing a (single) specified vertex.
##'
##' @title Maximal independent sets
##' @param A adjacency matrix of a graph
##' @param exclude sets that include these vertices are excluded.
##' @param include only sets that include this vertex are included.
##' @param subgraph the subset of vertices to search
##' @return a list of maximal independent sets
##' @export
maximal_IS <- function(A,exclude=c(),include=c(),subgraph=1:nrow(A)) {

  n <- nrow(A)
  isets <- list()
  marked <- logical(n)
  marked[exclude] <- TRUE

  dfs <- function(set,marked) {
    ## All independent vertices
    vs <- which(set==0L)
    if(!length(vs)) {
      isets[[length(isets)+1L]] <<- which(set==1L)
    } else {
      vs <- vs[!marked[vs]]
      for(v in vs) {
        s <- set
        s[A[,v]!=0L] <- -1L
        s[v] <- 1L
        marked[v] <- TRUE
        dfs(s,marked)
      }
    }
  }
  
  set <- rep(-1L,n)
  set[subgraph] <- 0L
  if(length(include)>0L) {
    set[A[,include]!=0L] <- -1L
    set[include] <- 1L
  }
  dfs(set,marked)
  isets
}


##' Find the vertices reachable from a root vertex in a directed graph
##'
##' @title Reachable vertices
##' @param A an adjacency matrix of a directed graph
##' @param root the root vertex
##' @return a vector of vertices reachable from the root
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,0,0,1,0,0,0,0,1,1,0,0,1,1),4,4)
##' # Find all vertices reachable from vertex 1
##' reachable(A,1)
##' @export
reachable <- function(A,root) {
  reach <- logical(nrow(A))
  reach[root] <- TRUE
  repeat {
    connected <- rowSums(A[,reach,drop=FALSE]!=0L)!=0L
    if(!any(connected & !reach)) break
    reach[connected] <- TRUE
  }
  which(reach)
}


##' Find the connected components of a graph
##'
##' Returns a list of the connected components of a graph, represented
##' as vectors of vertices.
##'
##' @title Connected components
##' @param A an adjacency matrix of a graph
##' @return a list of vectors of the vertices in each connected component
##' @examples
##' # Adjacency matrix of a directed graph
##' A <- matrix(c(1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1),4,4)
##' # Find all connected components
##' connected_components(A)
##' @export
connected_components <- function(A) {
  components <- list()
  vs <- seq_len(nrow(A))
  while(length(vs)) {
    comp <- reachable(A,vs[1L])
    components[[length(components)+1L]] <- comp
    vs <- setdiff(vs,comp)
  }
  components
}


##' Generate elements of the inverse of a matrix
##'
##' Given `S=sign(X)`, `partial` generates a symbolic representation
##' for a specified subset of elements of the inverse or adjoint of `X`,
##' expressed in terms of the elements of the weight matrix `A=abs(X)`.
##' The functions [partial_latex()], [partial_R()] and [partial_cpp()]
##' render the symbolic representation returned by `partial` as LaTeX,
##' R or C++ code.
##'
##' `partialA` is an alternative implementation of `partial` that is
##' often faster when many elements of the inverse or adjoint are
##' required.  `partial` and `partialA` index cycles differently, and
##' so can produce different but equivalent expressions.
##'
##' @title Partial matrix inverse
##' @param S a signed adjacency matrix
##' @param i vector of row indices
##' @param j vector of column indices
##' @return a partial object for the partial inverse of a matrix with the
##' same sign structure as the signed adjacency matrix S.
##' @seealso [partial_latex()], [partial_R()], [partial_cpp()]
##' @examples
##' # Tridiagonal sign matrix
##' n <- 3
##' S <- diag(-1,n,n)+rbind(0,diag(1,n-1,n))+cbind(0,diag(1,n,n-1))
##' # Generate R function for all elements of the inverse
##' cat(partial_R(partial(S,row(S),col(S))))
##' # Generate latex for all elements of the inverse
##' cat(partial_latex(partial(S,row(S),col(S))))
##' @export
partial <- function(S,i=integer(0),j=integer(0)) {

  ## Find all cycles and loops
  cycs <- directed_cycles(S)
  ## Adjoin missing loops
  zs <- unname(which(diag(S)==0L))
  zero <- length(cycs) + seq_along(zs)
  cycs <- c(cycs,as.list(zs))
  ## Form cycle graph
  C <- intersection_graph_int(cycs,nrow(S))
  ## Find MIS for each connected component of C
  ccmps <- connected_components(C)
  isets <- lapply(ccmps,function(cmp) maximal_IS(C,exclude=zero,subgraph=cmp))
  ## Vertices in each strongly connected component of S
  vcmps <- lapply(ccmps,function(cmp) sort(unique(unlist(cycs[cmp]))))
  
  pterm <- function(path) {
    ## Determinant factors for the inverse and adjoint
    inv <- vapply(vcmps,function(cmp) any(path %in% cmp),logical(1L))
    adj <- which(!inv)
    inv <- which(inv)
    ## Cycle weights for the path from each component
    wgt <- lapply(ccmps[inv],function(cmp) {
      cmp <- cmp[vapply(cycs[cmp],function(cyc) !any(path %in% cyc),logical(1L))]
      maximal_IS(C,exclude=zero,subgraph=cmp)
    })
    list(path=path,wgt=wgt,inv=inv,adj=adj)
  }
  
  inverse <- function(i,j) {
    if(is.character(i)) i <- match(i,rownames(S))
    if(is.character(j)) j <- match(j,rownames(S))
    ts <- lapply(find_paths(S,j,i),pterm)
    list(i=i,j=j,terms=ts)
  }
  
  els <- .mapply(inverse,list(i,j),NULL)
  if(length(zero) > 0L) cycs <- cycs[-zero]
  structure(list(S=S,components=vcmps,cycles=cycs,determinants=isets,elements=els),class="partial")
}


##' @rdname partial
##' @export
partialA <- function(S,i=integer(0),j=integer(0)) {

  ## Find all cycles
  cycs <- directed_cycles(S,loops=FALSE)
  ## Adjoin loops
  loops <- unname(c(which(diag(S)!=0L),which(diag(S)==0L)))
  idx <- length(cycs) + order(loops)
  zero <- as.integer(diag(S)==0L)
  cycs <- c(cycs,as.list(loops))
  ## Form cycle graph
  C <- intersection_graph_int(cycs,nrow(S))
  ## Find MIS for each connected component of C
  ccmps <- connected_components(C)
  isets <- lapply(ccmps,function(cmp) maximal_IS(C,subgraph=cmp))
  ## Vertices in each strongly connected component of S
  vcmps <- lapply(ccmps,function(cmp) sort(unique(unlist(cycs[cmp]))))
  ## Index loops in isets
  index <- lapply(isets,function(isets) t(vapply(isets,function(iset) idx %in% iset,integer(nrow(S)))))

  pterm <- function(path) {
    ## Determinant factors for the inverse and adjoint
    inv <- vapply(vcmps,function(cmp) any(path %in% cmp),logical(1L))
    adj <- which(!inv)
    inv <- which(inv)
    ## Cycle weights for the path from each component
    wgt <- lapply(inv,function(k) {
      path <- path[path %in% vcmps[[k]]]
      if(length(path) == length(vcmps[[k]])) return(list(integer(0L)))
      p <- integer(nrow(S))
      p[path] <- 1L
      keep <- (index[[k]]%*%p==length(path)) & (index[[k]]%*%(zero*(1L-p))==0L)
      pcycs <- idx[path]
      lapply(isets[[k]][keep],function(iset) iset[!(iset %in% pcycs)])
    })
    list(path=path,wgt=wgt,inv=inv,adj=adj)
  }

  inverse <- function(i,j) {
    if(is.character(i)) i <- match(i,rownames(S))
    if(is.character(j)) j <- match(j,rownames(S))
    ts <- lapply(find_paths(S,j,i),pterm)
    list(i=i,j=j,terms=ts)
  }

  els <- .mapply(inverse,list(i,j),NULL)
  if(sum(zero) > 0L) {
    isets <- lapply(seq_along(isets),function(k) isets[[k]][index[[k]]%*%zero==0L])
    cycs <- cycs[seq_len(length(cycs)-sum(zero))]
  }
  structure(list(S=S,components=vcmps,cycles=cycs,determinants=isets,elements=els),class="partial")
}



## This version retains zero cycles.
partial1 <- function(S,i=integer(0),j=integer(0)) {

  ## Find all cycles
  cycs <- directed_cycles(S,loops=FALSE)
  ## Adjoin missing loops
  zero <- length(cycs) + unname(which(diag(S)==0L))
  cycs <- c(cycs,as.list(seq_len(nrow(S))))
  ## Form cycle graph
  C <- intersection_graph_int(cycs,nrow(S))
  ## Find MIS for each connected component of C
  ccmps <- connected_components(C)
  dets <- lapply(ccmps,function(cmp) maximal_IS(C,exclude=zero,subgraph=cmp))
  ## Vertices in each strongly connected component of S
  vcmps <- lapply(ccmps,function(cmp) sort(unique(unlist(cycs[cmp]))))

  pterm <- function(path) {
    ## Determinant factors for the inverse and adjoint
    inv <- vapply(vcmps,function(cmp) any(path %in% cmp),logical(1L))
    adj <- which(!inv)
    inv <- which(inv)
    ## Cycle weights for the path from each component
    wgt <- lapply(ccmps[inv],function(cmp) {
      cmp <- cmp[vapply(cycs[cmp],function(cyc) !any(path %in% cyc),logical(1L))]
      maximal_IS(C,exclude=zero,subgraph=cmp)
    })
    list(path=path,wgt=wgt,inv=inv,adj=adj)
  }

  inverse <- function(i,j) {
    if(is.character(i)) i <- match(i,rownames(S))
    if(is.character(j)) j <- match(j,rownames(S))
    ts <- lapply(find_paths(S,j,i),pterm)
    list(i=i,j=j,terms=ts)
  }

  els <- .mapply(inverse,list(i,j),NULL)
  structure(list(S=S,components=vcmps,cycles=cycs,determinants=dets,elements=els),class="partial")
}


## This version retains zero cycles.
partialA1 <- function(S,i=integer(0),j=integer(0)) {

  ## Find all cycles
  cycs <- directed_cycles(S,loops=FALSE)
  ## Adjoin loops
  idx <- length(cycs) + seq_len(nrow(S))
  zero <- as.integer(diag(S)==0L)
  cycs <- c(cycs,as.list(seq_len(nrow(S))))
  ## Form cycle graph
  C <- intersection_graph_int(cycs,nrow(S))
  ## Find MIS for each connected component of C
  ccmps <- connected_components(C)
  isets <- lapply(ccmps,function(cmp) maximal_IS(C,subgraph=cmp))
  ## Vertices in each strongly connected component of S
  vcmps <- lapply(ccmps,function(cmp) sort(unique(unlist(cycs[cmp]))))
  ## Index loops in isets
  index <- lapply(isets,function(isets) t(vapply(isets,function(iset) idx %in% iset,integer(nrow(S)))))

  pterm <- function(path) {
    ## Determinant factors for the inverse and adjoint
    inv <- vapply(vcmps,function(cmp) any(path %in% cmp),logical(1L))
    adj <- which(!inv)
    inv <- which(inv)
    ## Cycle weights for the path from each component
    wgt <- lapply(inv,function(k) {
      path <- path[path %in% vcmps[[k]]]
      if(length(path) == length(vcmps[[k]])) return(list(integer(0L)))
      p <- integer(nrow(S))
      p[path] <- 1L
      keep <- (index[[k]]%*%p==length(path)) & (index[[k]]%*%(zero*(1L-p))==0L)
      pcycs <- idx[path]
      lapply(isets[[k]][keep],function(iset) iset[!(iset %in% pcycs)])
    })
    list(path=path,wgt=wgt,inv=inv,adj=adj)
  }

  inverse <- function(i,j) {
    if(is.character(i)) i <- match(i,rownames(S))
    if(is.character(j)) j <- match(j,rownames(S))
    ts <- lapply(find_paths(S,j,i),pterm)
    list(i=i,j=j,terms=ts)
  }

  els <- .mapply(inverse,list(i,j),NULL)
  if(length(zero) > 0L)
    isets <- lapply(seq_along(isets),function(k) isets[[k]][index[[k]]%*%zero==0L])
  structure(list(S=S,components=vcmps,cycles=cycs,determinants=isets,elements=els),class="partial")
}



##' Print a partial inverse object
##'
##' @title Print partial inverse
##' @param x a partial inverse
##' @param ... ignored
##' @export
print.partial <- function(x,...) {
  cat(sprintf("Partial inverse of %d x %d matrix\n",nrow(x$S),ncol(x$S)))
}


##' Render components of an partial inverse
##'
##' This is the workhose function that generates expressions to
##' evaluate the elements of an inverse or adjoint matrix.
##'
##' The `fmt` argument is a list of functions and strings that
##' construct syntactic elements of the generated code:
##'
##' * `asub(i,j)` function to generate the i,j element of the weight matrix
##' * `rsub(i,j)` function to generate the i,j element of the result matrix
##' * `csub(k)` function to generate the k-th cycle weight variable
##' * `dsub(k)` function to generate the k-th determinant variable
##' * `cassign(l,r)` vectorized function to generate cycle weight assignments (l=r)
##' * `dassign(l,r)` vectorized function to generate determinant assignments (l=r)
##' * `rassign(l,r)` vectorized function to generate result assignments (l=r)
##' * `ratio(n,d)` function to generate the ratio operation (n/d)
##' * `mul` string for the multiplication operator
##'
##' @title Render components of a partial inverse
##' @param partial a partial inverse object
##' @param fmt a list of sprintf formatters
##' @param type whether the inverse, adjoint or determinant is generated
##' @param determinant should uneccessary factors be retained
##' @return a list of character vectors that define
##'   * `cyc` assignments to cycle weight variables
##'   * `det` assignments to determinant variables
##'   * `el` assignments to elements of the inverse or adjoint matrix (if type!="det")
##' @seealso [partial_latex()], [partial_R()], [partial_cpp()]
##' @export
partial_components <- function(partial,fmt,type=c("inv","adj","det"),determinant=TRUE) {
  
  if(!inherits(partial,"partial")) stop("partial must be a partial object")
  type <- match.arg(type)
  if(type=="det") determinant <- TRUE
  
  ## Only keep needed determinants and cycles
  if(determinant)
    keep.det <- seq_along(partial$determinants)
  else
    keep.det <- sort(unique(unlist(lapply(partial$elements,function(el) lapply(el$terms,`[`,type)))))
  keep.cyc <- sort(unique(c(
    unlist(partial$determinants[keep.det]),
    unlist(lapply(partial$elements,function(el) lapply(el$terms,`[`,"wgt"))))))
  
  ## Render cycle weights
  cnames <- fmt$csub(keep.cyc)
  cyc <- vapply(partial$cycles[keep.cyc],render_cycle,character(1L),S=partial$S,fmt=fmt)
  cyc <- setNames(fmt$cassign(cnames,cyc),cnames)
  
  ## Render determinants
  dnames <- fmt$dsub(keep.det)
  det <- vapply(partial$determinants[keep.det],render_Cexpr,character(1L),fmt=fmt)
  det <- setNames(fmt$dassign(dnames,det),dnames)
  
  ## Render elements
  if(type=="det") {
    list(cyc=cyc,det=det)
  } else {
    lhs <- vapply(partial$elements,function(el) fmt$rsub(el$i,el$j),character(1L))
    if(type=="inv")
      rhs <- vapply(partial$elements,render_inverse_element,character(1L),S=partial$S,fmt=fmt)
    else
      rhs <- vapply(partial$elements,render_adjoint_element,character(1L),S=partial$S,fmt=fmt)
    el <- setNames(fmt$rassign(lhs,rhs),lhs)
    list(cyc=cyc,det=det,el=el)
  }
}


##' Render a `partial` object as LaTeX expressions for the partial
##' inverse or adjoint matrix.
##'
##' This function translates a `partial` object generated by
##' [partial()] to LaTeX code for the elements of the inverse or
##' adjoint matrix.  In the generated LaTeX,
##'
##' * the `a` are (unsigned) edge weights,
##' * the `c` are cycle weights,
##' * the `d` are determinants of strongly connected components
##'
##' If `determinant=FALSE`, any `D` variables not required for the
##' subsequent calculation areomitted.  If `determinant=TRUE`, no
##' variables are omitted and the product of the `D` variables is the
##' determinant.
##'
##' @title Render partial object as LaTeX
##' @param partial a partial inverse object
##' @param type whether the inverse or adjoint is generated
##' @param determinant should uneccessary factors be retained
##' @param use.names should variable names be used as subscripts
##' @return a character string of LaTeX code.
##' @seealso [partial()], [highlight_ncycle()]
##' @export
partial_latex <- function(partial,type=c("inv","adj"),determinant=TRUE,use.names=FALSE) {

  type <- match.arg(type)
  
  ## Define syntactic features for LaTeX
  if(use.names) {
    nms <- rownames(partial$S)
    asub <- function(i,j) sprintf("a_{\\mathrm{\\scriptsize %s},\\mathrm{\\scriptsize %s}}",nms[i],nms[j])
    rsub <- function(i,j) sprintf("r_{\\mathrm{\\scriptsize %s},\\mathrm{\\scriptsize %s}}",nms[i],nms[j])
  } else {
    asub <- function(i,j) sprintf("a_{%d,%d}",i,j)
    rsub <- function(i,j) sprintf("r_{%d,%d}",i,j)
  }
  fmt <- list(asub=asub,rsub=rsub,
              csub=function(k) sprintf("c_{%d}",k),
              dsub=function(k) sprintf("d_{%d}",k),
              cassign=function(l,r) sprintf("  %s &= %s\\\\",l,r) ,
              dassign=function(l,r) sprintf("  %s &= %s\\\\",l,r),
              rassign=function(l,r) sprintf("  %s &= %s\\\\",l,r),
              ratio=function(n,d) sprintf("%s%s^{-1}",n,d),
              mul="",recip="%s^{-1}",ratio="%s%s^{-1}")
  
  ## Generate code components
  code <- partial_components(partial,fmt,type,determinant)

  paste("\\begin{aligned}",
        paste(code$cyc,collapse="\n"),
        paste(code$det,collapse="\n"),
        paste(code$el,collapse="\n"),
        "\\end{aligned}",
        sep="\n")
}


##' Render a `partial` object as R code to calculate the partial
##' inverse or adjoint matrix.
##'
##' This function translates a `partial` object generated by
##' [partial()] to R code to evaluate the elements of the inverse or
##' adjoint matrix.  In the generated code,
##'
##' * the argument `A` is a matrix of (unsigned) edge weights,
##' * the `C` variables are cycle weights,
##' * the `D` variables are determinants of strongly connected components.
##'
##' If `determinant=FALSE`, any `D` variables not required for the
##' subsequent calculation areomitted.  If `determinant=TRUE`, no
##' variables are omitted and the product of the `D` variables is the
##' determinant.
##'
##' @title Render partial object as R code
##' @param partial a partial inverse object
##' @param type whether the inverse or adjoint is generated
##' @param determinant should uneccessary factors be retained
##' @return a character string of R code.
##' @seealso [partial()]
##' @export
partial_R <- function(partial,type=c("inv","adj"),determinant=TRUE) {
  
  type <- match.arg(type)
  ## Define syntactic features for R
  fmt <- list(asub=function(i,j) sprintf("A[%d,%d]",i,j),
              rsub=function(i,j) sprintf("R[%d,%d]",i,j),
              csub=function(k) sprintf("C%d",k),
              dsub=function(k) sprintf("D%d",k),
              cassign=function(l,r) sprintf("  %s <- %s",l,r),
              dassign=function(l,r) sprintf("  %s <- %s",l,r),
              rassign=function(l,r) sprintf("  %s <- %s",l,r),
              ratio=function(n,d) sprintf("%s/%s",n,d),
              mul="*")
  
  ## Generate code components
  code <- partial_components(partial,fmt,type,determinant)

  paste("function(A){",
        if(type=="inv") "  ## Inverse matrix" else "  ## Adjoint matrix",
        sprintf("  R <- matrix(NA,%d,%d,dimnames=dimnames(A))",nrow(partial$S),ncol(partial$S)),
        "  ## Cycle weights",
        paste(code$cyc,collapse="\n"),
        "  ## Determinants",
        paste(code$det,collapse="\n"),
        "  ## Elements",
        paste(code$el,collapse="\n"),
        "  R",
        "}",
        sep="\n")
}


##' Render a `partial` object as Rcpp code to calculate the partial
##' inverse or adjoint matrix.
##'
##' This function translates a `partial` object generated by
##' [partial()] to Rcpp code to evaluate the elements of the inverse
##' or adjoint matrix.  In the generated code,
##'
##' * the argument `A` is a matrix of (unsigned) edge weights,
##' * the `C` variables are cycle weights,
##' * the `D` variables are determinants of strongly connected components.
##'
##' If `determinant=FALSE`, any `D` variables not required for the
##' subsequent calculation areomitted.  If `determinant=TRUE`, no
##' variables are omitted and the product of the `D` variables is the
##' determinant.
##'
##' @title Render partial object as R code
##' @param partial a partial inverse object
##' @param type whether the inverse or adjoint is generated
##' @param determinant should uneccessary factors of the determinant be retained
##' @param fname the function name
##' @return a character string of R code.
##' @seealso [partial()]
##' @export
partial_cpp <- function(partial,type=c("inv","adj"),determinant=TRUE,fname=type) {

  ## Define syntactic features for c++
  fmt <- list(asub=function(i,j) sprintf("A(%d,%d)",i-1,j-1),
              rsub=function(i,j) sprintf("R(%d,%d)",i-1,j-1),
              csub=function(k) sprintf("C%d",k),
              dsub=function(k) sprintf("D%d",k),
              cassign=function(l,r) sprintf("  const double %s = %s;",l,r),
              dassign=function(l,r) sprintf("  const double %s = %s;",l,r),
              rassign=function(l,r) sprintf("  %s = %s;",l,r),
              ratio=function(n,d) sprintf("%s/%s",n,d),
              mul="*")
  
  ## Generate code components
  code <- partial_components(partial,fmt,type,determinant)
  
  paste(sprintf("NumericMatrix %s(const NumericMatrix& A) {",fname),
        sprintf("  NumericMatrix R(%d, %d);",nrow(partial$S),ncol(partial$S)),
        "  std::fill(R.begin(), R.end(), NumericMatrix::get_na());",
        "  // Cycle weights",
        paste(code$cyc,collapse="\n") ,
        "  // Determinants",
        paste(code$det,collapse="\n"),
        "  // Elements",
        paste(code$el,collapse="\n"),
        "  return R;",
        "}",
        sep="\n")
}

##' Render a `partial` object as LaTeX expressions for the determinant.
##'
##' This function translates a `partial` object generated by
##' [partial()] to LaTeX code for the determinant.  In the generated LaTeX,
##'
##' * the `a` are (unsigned) edge weights,
##' * the `c` are cycle weights,
##' * the `d` are determinants of strongly connected components
##'
##' The determinant of the matrix (d_0) is the product of the determinants of
##' the strongly connected components.
##'
##' @title Render determinant as LaTeX
##' @param determinant a determinant object
##' @param use.names should variable names be used as subscripts
##' @return a character string of LaTeX code.
##' @seealso [partial()]
##' @export
determinant_latex <- function(determinant,use.names=FALSE) {

  ## Define syntactic features for LaTeX
  if(use.names) {
    nms <- rownames(determinant$S)
    asub <- function(i,j) sprintf("a_{\\mathrm{\\scriptsize %s},\\mathrm{\\scriptsize %s}}",nms[i],nms[j])
    rsub <- function(i,j) sprintf("r_{\\mathrm{\\scriptsize %s},\\mathrm{\\scriptsize %s}}",nms[i],nms[j])
  } else {
    asub <- function(i,j) sprintf("a_{%d,%d}",i,j)
    rsub <- function(i,j) sprintf("r_{%d,%d}",i,j)
  }
  fmt <- list(asub=asub,rsub=rsub,
              csub=function(k) sprintf("c_{%d}",k),
              dsub=function(k) sprintf("d_{%d}",k),
              cassign=function(l,r) sprintf("  %s &= %s\\\\",l,r),
              dassign=function(l,r) r,
              ratio=function(n,d) sprintf("%s%s^{-1}",n,d),
              mul="")
  
  ## Generate code components
  code <- partial_components(determinant,fmt,type="det")
  
  ## Generate determinant
  if(length(code$det)==1L) {
    code$det <- sprintf("  d &= %s\\\\",code$det)
  } else {
    code$det <- c(
      sprintf("  %s &= %s\\\\",names(code$det),code$det),
      sprintf("  d &= %s\\\\",paste(names(code$det),collapse="")))
  }
  
  paste("\\begin{aligned}",
        paste(code$cyc,collapse="\n"),
        paste(code$det,collapse="\n"),
        "\\end{aligned}",
        sep="\n")
}


##' Render a `partial` object as R code to calculate the determinant.
##'
##' This function translates a `partial` object generated by
##' [partial()] to R code to evaluate the determinant.  In the
##' generated code,
##'
##' * the argument `A` is a matrix of (unsigned) edge weights,
##' * the `C` variables are cycle weights,
##' * the `D` variables are determinants of strongly connected components.
##'
##' The determinant of the matrix is the product of the determinants
##' of the strongly connected components.
##'
##' @title Render determinant as R code
##' @param determinant a determinant object
##' @return a character string of R code.
##' @seealso [partial()]
##' @export
determinant_R <- function(determinant) {

  ## Define syntactic features for R
  fmt <- list(asub=function(i,j) sprintf("A[%d,%d]",i,j),
              rsub=function(i,j) sprintf("R[%d,%d]",i,j),
              csub=function(k) sprintf("C%d",k),
              dsub=function(k) sprintf("D%d",k),
              cassign=function(l,r) sprintf("  %s <- %s",l,r),
              dassign=function(l,r) r,
              ratio=function(n,d) sprintf("%s/%s",n,d),
              mul="*")
  
  ## Generate code components
  code <- partial_components(determinant,fmt,type="det")
  
  
  ## Generate determinant
  if(length(code$det)==1L) {
    code$det <- sprintf("  %s",code$det)
  } else {
    code$det <- c(
      sprintf("  %s <- %s",names(code$det),code$det),
      sprintf("  %s",paste(names(code$det),collapse="*")))
  }
  
  paste("function(A){",
        "  ## Cycle weights",
        paste(code$cyc,collapse="\n"),
        "  ## Determinants",
        paste(code$det,collapse="\n"),
        "}",
        sep="\n")
}




##' Characteristic polynomial of a matrix
##'
##' Given `S=sign(X)`, `cpoly` generates a symbolic representation of
##' the characteristic polynomial of the matrix `X`. The functions
##' [cpoly_latex()] and [cpoly_R()] render this symbolic
##' representation as LaTeX and R code.
##'
##' If `factor=TRUE`, the characteristic polynomial is factored into
##' separate components, one for each strongly connected component of
##' the corresponding graph.
##'
##' @title Characteristic polynomial
##' @param S a signed adjacency matrix
##' @param factor should the characteristic polynomial be factored
##' @return a `cpoly` object
##' @seealso [cpoly_latex()], [cpoly_R()]
##' @examples
##' # Tridiagonal sign matrix
##' n <- 3
##' S <- diag(-1,n,n)+rbind(0,diag(1,n-1,n))+cbind(0,diag(1,n,n-1))
##' # Generate R function for coeffs of the characteristic polynomial
##' cat(cpoly_R(cpoly(S)))
##' # Generate latex for the characteristic polynomial
##' cat(cpoly_latex(cpoly(S)))
##' @export
cpoly <- function(S,factor=TRUE) {

  ## Find cycles and adjoin index loops
  cycs <- directed_cycles(S)
  n <- length(cycs)
  cycs <- c(cycs,as.list(seq_len(nrow(S))))
  index <- n+seq_len(nrow(S))
  
  ## Cycle intersection graph
  C <- intersection_graph_int(cycs,nrow(S))
  
  
  if(factor) {
    ## Connected components
    ccmps <- connected_components(C)
    vcmps <- lapply(ccmps,function(cmp) sort(unique(unlist(cycs[cmp]))))
    
    ## Polynomial coefficients for each connected component
    coeffs <- lapply(seq_along(ccmps),function(k) {
      isets <- maximal_IS(C,subgraph=ccmps[[k]])
      ## Split by order and remove indices
      ord <- vapply(isets,function(iset) sum(iset %in% index),integer(1L))
      isets <- lapply(isets,setdiff,index)
      lapply(0:length(vcmps[[k]]),function(k) isets[ord==k])
    })
  } else {
    ## Treat as one connected component
    isets <- maximal_IS(C)
    ## Split by order and remove indices
    ord <- vapply(isets,function(iset) sum(iset %in% index),integer(1L))
    isets <- lapply(isets,setdiff,index)
    coeffs <- list(lapply(0:nrow(S),function(k) isets[ord==k]))
  }
  
  structure(list(S=S,cycles=cycs[-index],coeffs=coeffs),class="cpoly")
}


##' Print a characteristic polynomial object
##'
##' @title Print characteristic polynomial
##' @param x a characteristic polynomial
##' @param ... ignored
##' @export
print.cpoly <- function(x,...) {
  cat(sprintf("Characteristic polynomial of %d x %d matrix\n",nrow(x$S),ncol(x$S)))
}


##' Render components of a `cpoly` object
##'
##' This is the workhorse function that generates expressions to
##' evaluate the coeffcients of the characteristic polynomial.
##'
##' The `fmt` argument is a list of functions and strings that
##' construct syntactic elements of the generated code:
##'
##' * `asub(i,j)` function to generate the i,j element of the weight matrix
##' * `csub(k)` function to generate the k-th cycle weight variable
##' * `bsub(k)` function to generate the k-th polynomial coefficient variable
##' * `cassign(l,r)` vectorized function to generate cycle weight assignments (l=r)
##' * `bassign(l,r)` vectorized function to generate polynomial coefficient assignments (l=r)
##' * `mul` string for the multiplication operator
##'
##' @title Render components of a characteristic polynomial
##' @param cpoly a `cpoly` object
##' @param fmt a list of format functions and strings
##' @return a list of character vectors that define
##'   * `cyc` assignments to cycle weight variables
##'   * `cft` lists of assignments to polynomial coefficients, one for each component
##' @seealso [cpoly()]
##' @export
cpoly_components <- function(cpoly,fmt) {

  if(!inherits(cpoly,"cpoly")) stop("cpoly must be a cpoly object")

  render_coeffs <- function(coeffs,n0) {
    n <- length(coeffs)
    bnames <- fmt$bsub(seq_len(n)+n0-1L)
    cft <- vapply(seq_len(n),function(k) {
      cf <- render_Cexpr(coeffs[[k]],fmt)
      if((n-k)%%2L!=0L)
        cf <- if(length(coeffs[[k]])>1L) sprintf("-(%s)",cf) else sprintf("-%s",cf)
      cf
    },character(1L))
    setNames(fmt$bassign(bnames,cft),bnames)
  }
  
  ## Render cycle weights
  cnames <- fmt$csub(seq_along(cpoly$cycles))
  cyc <- vapply(cpoly$cycles,render_cycle,character(1L),S=cpoly$S,fmt=fmt)
  cyc <- setNames(fmt$cassign(cnames,cyc),cnames)
  
  ## Render coefficients
  n0 <- c(0,cumsum(lengths(cpoly$coeffs)[-length(cpoly$coeffs)]))
  cft <- .mapply(render_coeffs,list(cpoly$coeffs,n0),NULL)
  list(cyc=cyc,cft=cft)
}


##' Render a `cpoly` object as a LaTeX expression for the
##' characteristic polynomial.
##'
##' This function translates a `cpoly` object generated by [cpoly()]
##' to LaTeX expressions for the characteristic polynomial.  In the
##' generated LaTeX,
##'
##' * the `a` are (unsigned) edge weights,
##' * the `c` are cycle weights,
##' * the `b` are coefficients of the characteristic polynomial.
##'
##' @title Render a `cpoly` object as LaTeX
##' @param cpoly a `cpoly` object
##' @param use.names should variable names be used as subscripts
##' @return a character string of LaTeX code
##' @seealso [cpoly()]
##' @export
cpoly_latex <- function(cpoly,use.names=FALSE) {

  ## Define syntactic features for LaTeX
  if(use.names) {
    nms <- rownames(cpoly$S)
    asub <- function(i,j) sprintf("a_{\\mathrm{\\scriptsize %s},\\mathrm{\\scriptsize %s}}",nms[i],nms[j])
  } else {
    asub <- function(i,j) sprintf("a_{%d,%d}",i,j)
  }
  fmt <- list(asub=asub,
              csub=function(k) sprintf("c_{%d}",k),
              bsub=function(k) sprintf("b_{%d}",k),
              cassign=function(l,r) sprintf("  %s &= %s\\\\",l,r),
              bassign=function(l,r) sprintf("  %s &= %s\\\\",l,r),
              mul="")

  ## Generate code components
  code <- cpoly_components(cpoly,fmt)
  
  ## Generate polynomial expressions
  polynomial <- function(bs) {
    if(length(bs)==1L) return(names(bs))
    if(length(bs)==2L) return(sprintf("%s+%sx",names(bs)[1],names(bs)[2]))
    paste(c(names(bs)[1],
            sprintf("%sx",names(bs)[2]),
            sprintf("%sx^{%d}",names(bs)[-(1:2)],seq.int(2,length(bs)-1))),
          collapse="+")
  }
  poly <- sprintf("  p_{%d}(x) &= %s\\\\",seq_along(code$cft),vapply(code$cft,polynomial,character(1L)))
  
  ## Generate polynomial
  bnames <- names(code$cft)
  paste("\\begin{aligned}",
        paste(code$cyc,collapse="\n"),
        paste(unlist(code$cft),collapse="\n"),
        paste(poly,collapse="\n"),
        "\\end{aligned}",
        sep="\n")
}



##' Render a `cpoly` object as R code to evaluate the
##' coefficients of the characteristic polynomial.
##'
##' This function translates a `cpoly` object generated by [cpoly()]
##' to R code to evaluate the coefficients of the characteristic
##' polynomial.  In the generated R code,
##'
##' * the argument `A` is a matrix of (unsigned) edge weights,
##' * the `C` variables are cycle weights,
##' * the `B` variables are coefficients of the characteristic polynomial.
##'
##' @title Render cpoly object as R code
##' @param cpoly a `cpoly` object
##' @return a character string of R code
##' @seealso [cpoly()]
##' @export
cpoly_R <- function(cpoly) {

  ## Define syntactic features for R
  fmt <- list(asub=function(i,j) sprintf("A[%d,%d]",i,j),
              csub=function(k) sprintf("C%d",k),
              bsub=function(k) sprintf("B%d",k),
              cassign=function(l,r) sprintf("  %s <- %s",l,r),
              bassign=function(l,r) sprintf("  %s <- %s",l,r),
              mul="*")
  
  ## Generate code components
  code <- cpoly_components(cpoly,fmt)
  
  ## Generate coefficient vectors
  ## Generate polynomial expressions
  bvector <- function(bs)
    paste("c(",paste(names(bs),collapse=","),")",sep="")
  ret <- vapply(code$cft,bvector,character(1L))
  ret <- paste("  list(",paste(sprintf("p%d=%s",seq_along(ret),ret),collapse=",\n       "),")",sep="")
  
  
  paste("function(A){",
        "  ## Cycle weights",
        paste(code$cyc,collapse="\n"),
        "  ## Coefficients",
        paste(unlist(code$cft),collapse="\n"),
        ret,
        "}",
        sep="\n")
}




##' Highlight negative cycle weights in LaTeX
##'
##' The `highlight_negative` function highlights negative cycle
##' weights in LaTeX code generated by [partial_latex()].
##'
##' @title Highlight negative cycles in LaTeX
##' @param latex a character string of LaTeX code
##' @param col the highlight colour
##' @return a character string of LaTeX code
##' @seealso [partial_latex()]
##' @export
highlight_ncycle <- function(latex,col="red") {
  m <- regmatches(latex,gregexec("c_\\{(\\d+)\\} &= -",latex))[[1]]
  if(length(m)==0L) {
    latex
  } else {
    gsub(sprintf("(c_\\{(%s)\\})",paste(m[2L,],collapse="|")),
         sprintf("\\\\color{%s}{\\1}",col),
         latex)
  }
}


##' Evaluate a polynomial
##'
##' This function evaluates a polynomial with coefficients `cs` at
##' values `x`.  If `x` is a matrix, then `x` must be square and the
##' matrix polynomial is evaluated.
##'
##' @title Evaluate a polynomial
##' @param cs a vector of polynomial coefficients
##' @param x a vector or matrix of values
##' @return a vector or matrix of polynomial values
##' @export
eval_poly <- function(cs,x) {
  n <- length(cs)
  if(is.matrix(x)) {
    p <- diag(cs[n],nrow(x),nrow(x))
    for(k in rev(seq_len(n-1L))) {
      p <- x%*%p
      diag(p) <- diag(p)+cs[k]
    }
  } else {
    p <- rep(cs[n],length(x))
    for(k in rev(seq_len(n-1L))) p <- cs[k]+x*p
  }
  p
}






##' Utility functions for rendering partial inverse objects
##'
##' @title Render R utility functions
##' @param S the signed adjacency matrix
##' @param fmt a list of format functions and strings
##' @param path a vector of vertices
##' @param cycle a vector of vertices
##' @param cexpr a vector of cycle numbers
##' @param sign a sign character
##' @param a,b (unsigned) multiplicands
##' @param element an element of a partial object
##' @param term a term of an element
##' @param dets indices of component determinants
##' @return a character string
##' @rdname render_utils
render_path <- function(path,S,fmt) {
  if(length(path)==1L) return("-1")
  is <- path[-1]
  js <- path[-length(path)]
  s <- if((-1)^length(path)*prod(sign(S[cbind(is,js)]))>0L) "+" else "-"
  paste(s,paste(fmt$asub(is,js),collapse=fmt$mul),sep="")
}

##' @rdname render_utils
render_cycle <- function(cycle,S,fmt) {
  if(length(cycle)==0L) return("+1")
  is <- c(cycle[-1],cycle[1])
  js <- cycle
  s <- if((-1)^(length(cycle)+1L)*prod(sign(S[cbind(is,js)]))>0L) "" else "-"
  paste(s,paste(fmt$asub(is,js),collapse=fmt$mul),sep="")
}

##' @rdname render_utils
render_Cexpr <- function(cexpr,fmt,sign="+") {
  if(length(cexpr)==0L) return("0")
  s <- paste(
    vapply(cexpr,function(cset) {
      paste(sign,
            if(length(cset)==0L) "1" else paste(fmt$csub(cset),collapse=fmt$mul),
            sep="")
    },
    character(1L)),
    collapse="")
  if(substr(s,1L,1L)=="+") substr(s,2L,nchar(s)) else s
}

##' @rdname render_utils
render_product <- function(a,b,fmt) {
  if(a=="1") return(b)
  if(b=="1") return(a)
  paste(a,b,sep=fmt$mul)
}

##' @rdname render_utils
render_Dproduct <- function(dets,fmt) {
  if(length(dets)==0L) return("1")
  if(length(dets)==1L) return(fmt$dsub(dets))
  paste("(",
        paste(vapply(dets,function(k) fmt$dsub(k),character(1L)),collapse=fmt$mul),
        ")",sep="")
}


##' @rdname render_utils
render_wpath <- function(term,S,fmt) {
  ## Form numerator, coalescing signs into path
  sgn <- -1
  cyc <- "1"
  for(w in term$wgt) {
    cw <- render_Cexpr(w,fmt)
    ## Add parentheses if more than one product of cycles
    if(length(w)==1L) {
      if(substr(cyc,1L,1L)=="-") {
        sgn <- -sgn
        cw <- substr(cw,2L,nchar(cw))
      }
    } else {
      cw <- sprintf("(%s)",cw)
    }
    cyc <- render_product(cyc,cw,fmt)
  }
  ## Multiply by path
  path <- render_path(term$path,S,fmt)
  if(substr(path,1L,1L)=="-") sgn <- -sgn
  path <- substr(path,2L,nchar(path))
  list(sgn=if(sgn==-1L) "-" else "+",
       wgt=render_product(path,cyc,fmt))
}


##' @rdname render_utils
render_inverse_element <- function(element,S,fmt) {
  
  inverse_term <- function(term) {
    if(all(lengths(term$wgt)==0L)) return("0")
    sn <- render_wpath(term,S,fmt)
    dprod <- render_Dproduct(term$inv,fmt)
    if(dprod=="1") {
      paste(sn$sgn,sn$wgt,sep="")
    } else {
      paste(sn$sgn,fmt$ratio(sn$wgt,dprod),sep="")
    }
  }
  
  tms <- vapply(element$terms,inverse_term,character(1L))
  tms <- tms[tms!="0"]
  s <- if(length(tms)==0L) "0" else paste(tms,collapse="")
  if(substr(s,1L,1L)=="+") substr(s,2L,nchar(s)) else s
}

##' @rdname render_utils
render_adjoint_element <- function(element,S,fmt) {
  
  adjoint_term <- function(term) {
    if(all(lengths(term$wgt)==0L)) return("0")
    sn <- render_wpath(term,S,fmt)
    dprod <- render_Dproduct(term$adj,fmt)
    paste(sn$sgn,render_product(sn$wgt,dprod,fmt),sep="")
  }
  
  tms <- vapply(element$terms,adjoint_term,character(1L))
  tms <- tms[tms!="0"]
  s <- if(length(tms)==0L) "0" else paste(tms,collapse="")
  if(substr(s,1L,1L)=="+") substr(s,2L,nchar(s)) else s
}

