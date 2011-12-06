\name{addBoundary}
\alias{addBoundary}
\alias{addBoundary,ANY-method}
\title{Embed matrix within larger matrix}
\description{
  This function is useful for taking the grid of 
  densities/counts as a matrix and embedding
   it in a larger matrix whose additional values
  are just the mean of our original matrix.
  This new matrix can then be handed to
  \code{\link{cartogram}} to fit the 
  diffusion cartogram. The extra "padding" or "sea"
  helps the algorithm converge.
}
\usage{
addBoundary(pop, sea = 2, land.mean = mean(unlist(pop)))
}
\arguments{
  \item{pop}{the matrix of densities/counts to be embedded.}
  \item{sea}{a number which is applied to the dimensions of the original matrix
    to determine how many rows and columns should be added to the original matrix.
    Twice as many rows and columns are added as this product,
    an equal number of rows above and below, and an equal  number of columns to the left and to the right.
  }
  \item{land.mean}{the scalar value that we use as the value for the
     cells we with which we surround the existing matrix. This is the
     value for all the cells in the sea region.}
}
\value{
 A matrix with dimension
  \code{ floor((1 + 2 * seq) * dim(pop)) }
 and whose "outer" values in the added rows is the mean of 
 the original matrix \code{pop}.
}
\references{
"Diffusion-based method for producing density equalizing maps", 
 Michael T. Gastner and M. E. J. Newman, Proc. Natl. Acad. Sci. USA 101, 7499-7504 (2004)
 \url{http://www.pnas.org/cgi/content/abstract/101/20/7499}}

\author{Duncan Temple Lang created the interface for R}

\seealso{ 
 \code{\link{cartogram}}
}
\examples{

 m = matrix(1, 5, 7)
 mm = col(m)

 ex = addBoundary(mm)
 dim(ex)
 image(ex)

 ex = addBoundary(mm, 1)
 dim(ex)
}
\keyword{hplot}
\keyword{manip}

