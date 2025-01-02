
# quoll

<!-- badges: start -->
<!-- badges: end -->

quoll uses graph theoretic techniques to derive symbolic representations of elements of the inverse
or adjoint of a matrix $X$, given only the sign pattern of $X$.  

The complexity of the generated symbolic representations is dictated by the number of cycles in the 
directed graph for has $X$ as its weighted adjacency matrix, and so `quoll` is only practical for
relatively small, sparse matrices.

## Installation

The current version of quoll can be installed from GitHub using the remotes package. 
```r
# install.packages("remotes")
remotes::install_github("SWotherspoon/quoll")
```

## Example

This is a basic example which shows you how to solve a common problem:

```r
library(quoll)
## Define a sign matrix
S <- matrix(c(-1,1,0,0,1,-1,1,0,0,1,-1,1,0,0,1,-1),4,4)
## Calculate the (1,2) and the (4,3) elements of the inverse
pinv <- partial(S,c(1,4),c(2,3))
## Render as LaTeX
cat(render_latex(pinv))
```

The gnerated expressions can be rendered as LaTeX, or as R, c++ or Fortran code.



