cartogram =
function(pop, zero = TRUE, blur = 0.0, sea = NA)
{
      # Force this to be a matrix (at least for now. A data frame is not contiguous)
   pop = as.matrix(pop)

   if(!is.na(sea)) 
         # Add some padding with the constant value of the mean of the original matrix.
      pop = addBoundary(pop, sea)


      # we reverse this here as we are treating the columns as y's in the C code.
   dim = rev(dim(pop))

if(TRUE)  {
   x = expand.grid(1:(ncol(pop) + zero), 1:(nrow(pop) + zero))
   y = x[,2]
   x = x[,1]
} else     
      # Currently ignored.
      # We were letting the caller pass their own grid.
   if(missing(y) && (is.matrix(x) || is.data.frame(x))) {
	y = x[, 2]
        x = x[, 1]
   }
   
  if(is.matrix(pop)) {
         # Could do lapply(seq(length = ncol(pop)), function(i) as.numeric(pop[, i])), but no need.
      popEls = pop
  } else if(is.data.frame(pop)) {
       # Turn the pop into a collection of column vectors.
       # At present, this won't happen. We appear to need the vectors to be contiguous
       # although the C code doesn't look as if it cares. But the fftw code might.
     popEls = lapply(pop, as.numeric)
  }

     # Zero based counting for C.
   x = as.numeric(x) - 1 
   y = as.numeric(y) - 1

   tmp = .Call("R_makecartogram", popEls, x, y, dim, as.numeric(blur), PACKAGE = "Rcartogram")
     # convert the results into matrices.
   ans = lapply(tmp, function(x) matrix(x, nrow(pop) + zero, ncol(pop) + zero, byrow = TRUE))
   ans = structure(ans,
	           names  = c("x", "y"),
	           class = "Cartogram")

    # add in information about the expanded sea if relevant.
   ans
}

mapToGrid =
function(x, y, grid, offsets = c(0, 0))
{
  size = dim(grid$grid)
  xs = seq(grid$x[1], grid$x[2], length = size[1])
  xs = xs + (xs[2] - xs[1])/2
  xt = sapply(x, function(val) which( val < xs)[1]) + offsets[1]
  ys = seq(grid$y[1], grid$y[2], length = size[2])
  ys = ys + (ys[2] - ys[1])/2
  yt = sapply(y, function(val) which( val < ys)[1])  + offsets[2]
  cbind(xt, yt)
}

image.ExpandedMatrix =
function(x, ...)
{
  image(seq(1, nrow(x)), seq(1, ncol(x)), x)
}

transform.Cartogram =
function(`_data`, x, y = NULL, ...)
{
  predict(`_data`, x, y, ...) 
}

predict.Cartogram =
function(object, x, y = NULL, ...)
{
  if(missing(y) || is.null(y)) {
    y = x[,2]
    x = x[,1]
  }

  if(length(x) != length(y)) {
    len = max(length(x), length(y))
    length(x) = len
    length(y) = len    
  }
     
  
  # Avoid problems with the same vector (tmp) being sent to C twice due to R's
  # copy-on-modify rules
  tmp_x = rep(as.numeric(NA), length(x))
  tmp_y = rep(as.numeric(NA), length(y))
  ans = list(x = tmp_x, y = tmp_y)
  
  .Call("R_predict", object, as.numeric(x),  as.numeric(y), ans, dim(object$x), PACKAGE = "Rcartogram")
  ans
}

if(FALSE)  # an in R version.
predict.Cartogram =
 #
 #  Essentials taken from Mark Newman's interp.c code.
 #
 #   This is a simple, inefficient version at present. 
 #
function(object, x, y = NULL, ...)
{
  if(missing(y)) {
    y = x[, 2]
    x = x[, 1]
  }
   
  ix = as.integer(x)
  iy = as.integer(y)
  dx = x - ix
  dy = y - iy

  # This could be done much more efficiently with some vectorized operations. 
  # and if we really care, it can be done very easily in C.

  sapply(object,  
          function(m) {
            sapply(seq(along = ix),
                    function(i)
                       pred(m, ix, iy, dx, dy, i))
            })
}

pred = 
function(m, ix, iy, dx, dy, i = 1) {
   # + need to go at the end of the line or the R parser
   # thinks these are separate expressions with the last 3 preceded by a +
   #  e.g.  1
   #        +2
   # which gives expression{ 1, 2}
 (1-dx[i])*(1-dy[i]) * m[, ix[i]][ iy[i] ] +
       dx[i] * (1-dy[i]) * m[, ix[i] + 1L][iy[i]]  + 
         (1-dx[i]) * dy[i]* m[, ix[i]][iy[i] + 1L]  + 
           dx[i] * dy[i]* m[, ix[i] + 1L][iy[i] + 1L]
}




setGeneric("addBoundary",
           function(pop, sea = 2, land.mean = mean(unlist(pop)))
             standardGeneric("addBoundary"))

#setMethod("addBoundary", c("CartogramGrid", mean = "missing"),
#         function(pop, sea = 2, mean = mean(unlist(pop)))
#            addBoundary(pop, sea, attr(pop, "mean"ff))
#         )

setMethod("addBoundary", "ANY",
  #
  #  Take a matrix and add padding around it in the form
  #  of  rows and columns so that the original matrix is contained
  #  within the center of this new matrix. The value for each of the
  #  pad cells we add is the mean of the original matrix.
  #
  #  sea is a factor that determines the number of rows and columns to add
  #  on each side of the original matrix.
  #   In other words, if we add 2 rows, then we add 2 above and 2 below
  #  for a total of 4.
  #
  # The return value has attributes that identify what rows and columns
  # were added and the rows and columns by which one can identify the original
  # submatrix.
function(pop, sea = 2, land.mean = mean(unlist(pop)))
{
      extra = as.integer( sea * dim(pop))

      pop = as.matrix(pop)  
      pad.top = matrix(land.mean, extra[1], ncol(pop) + 2*extra[2])
      pad.left = matrix(land.mean, nrow(pop), extra[2])
      structure(rbind(pad.top, 
                      cbind(pad.left, pop, pad.left),
                      pad.top),
                class = "ExpandedMatrix",
                extra = extra,
                orig = dim(pop),
                land = list(rows = seq(extra[1] + 1, length = nrow(pop)),
                            cols = seq(extra[2] + 1, length = ncol(pop))))
})
