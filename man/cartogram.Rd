\name{cartogram}
\alias{cartogram}
\title{Fit a diffusion cartogram}
\description{
  This uses Mark Newman's C code  to fit the diffusion cartogram 
  developed by  Michael Gastner and Mark Newman.

  One calls this with a matrix (or a data.frame)
  containing the counts or densities 
  (e.g. population, number of incidences of a particular disease, etc.)
  for that location.
  By default, the code returns a list of two matrices (x and y) each with one more row and column as there are
  in \code{pop}.
  These two "grids" can then be used to predict or compute the transformed or projected value 
  of a point in the new cartogram using the \code{predict} method in this package.

  This is a reasonably straightforward interface to Mark Newman's code.  Therefore
  to understand what it does and how to use, read the documentation
  at \url{http://www-personal.umich.edu/~mejn/cart/doc/}, specifically
  "Using the cart program".

  Note  that the use of the word "cart" has nothing to do with classification and regression trees.
}
\usage{
cartogram(pop, zero = TRUE, blur = 0, sea = NA)
}
\arguments{
  \item{pop}{a matrix.  This is oriented as Xs along the columns and Y's going down each column. 
    In other words, rows in the R matrix correspond to the same y element in each of the X vectors.
   }
  \item{zero}{a logical value that indicates whether to expand the grid by 1 more row and column. This is
   what is done by the cart application by Mark Newman so we repeat it here to allow for obtaining comparable results.}
  \item{blur}{a non-negative real value that is used in the C code to add  Gaussian noise to the diffusion.}
  \item{sea}{either a missing value (NA) or a multiplier that is passed to 
   \code{\link{addBoundary}} to embed the pop matrix within a larger
   grid. The values of the outer  grid are constant and  the mean
   of the original \code{pop} matrix is used.}
}
\details{
 See the paper in the references.
}
\value{
 A list with two entries, both 
 matrices with as many rows and columns as \code{pop}
 or the enlarged grid if \code{sea} is not an \code{NA}.
}
\references{
"Diffusion-based method for producing density equalizing maps", 
 Michael T. Gastner and M. E. J. Newman, Proc. Natl. Acad. Sci. USA 101, 7499-7504 (2004)
 \url{http://www.pnas.org/cgi/content/abstract/101/20/7499}}

\author{Duncan Temple Lang created the interface for R}
\seealso{ 
 \code{\link{addBoundary}}
}
\examples{
  filename = system.file("data", "uspop.dat", package = "Rcartogram")
  pop = read.table(filename)
  grid = cartogram(pop)

  predict(grid, 3.1, 4.5)
 
  # note that this is equivalent to using the interp application in
  # the cart(ogram) distribution with inputs as 2.1 3.4
  # i.e. echo "2.1 3.4" | interp 1024 512 output.dat
  #        

    # now predict lots of values.
  o = system.time(predict(grid, 
                          runif(10000, 1, 1024), 
                          runif(10000, 1, 512)))

}
\keyword{hplot}


