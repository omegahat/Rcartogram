library(Rcartogram)

N = 200
N2 = N/2
circle = FALSE # TRUE
plotInExpanded = TRUE

 # Create a matrix with different regions.
dat = matrix(0, N, N)

 # These are bottom left, top-left, bottom right, top right, middle.
densities = c(bl = 10, tl = 20, br = 30, tr = 40, middle = 50) * 10

# bottom left
dat[1:N2, 1:N2] = densities[1]
# top left - first 1:N2 columns
dat[(N2+1):N, 1:N2] = densities[2]
# bottom right
dat[1:N2, (N2 + 1):N] = densities[3]
# top right
dat[(N2+1):N, (N2+1):N] = densities[4]
# 
if(circle) {
      # a circle in the middle
   dat[ sqrt((col(dat) - N2)^2 + (row(dat) - N2)^2) < .1*N ] = densities[5]
} else
      # Or a rectangular boundary
  dat[ abs(col(dat) - N2) + abs(row(dat) - N2) < 10 ] = densities[5]

# Plot the image.
#image(1:N, 1:N, t(dat))

showData =
function(dat, xlim = c(1, ncol(dat)) + c(-10, 10), ylim = c(1, nrow(dat)) + c(-10, 10))
{  
  cols = rainbow(length(densities))
  names(cols) = as.character(densities)
  plot(1:N, 1:N, type = "n", xlim = xlim, ylim = ylim)
  invisible(sapply(1:N, function(i) points(rep(i, N), 1:N, col = cols[as.character(dat[,i])], pch = 22)))
}

showData(dat)

image(1:N, 1:N, t(dat))
text(N/4, N/4, t(dat)[N2, N/4])
text(N/4, 3*N/4, t(dat)[N2, 3*N/4])
text(3*N/4, N/4, t(dat)[3*N/4, N/4])
text(3*N/4, 3*N/4, t(dat)[3*N/4, 3*N/4])
text(N2, N2, t(dat)[N2, N2])

text(c(N/4, N/4, 3*N/4, 3*N/4, N/2),
     c(N/4, 3*N/4, N/4, 3*N/4, N/2),
     densities
     )

# Look at the data, 
dat[c(1, 101, 200), c(1, 101, 200)][3:1,]


#######################################

# Now compute the region borders.

tmp = apply(dat[1:N2, 1:(N2 + 1)], 1, function(x) min(which(x != densities[1]))) - 1
boundary.bottomleft =
 cbind( c(1:(N2-1),
          tmp,
          tmp[N2]:1),
        c(rep(1, N2-1), 1:N2, rep(N2, tmp[N2])))
polygon(boundary.bottomleft, lwd = 2)


 # Flip this for the top left

boundary.topleft = boundary.bottomleft
boundary.topleft[,2] = N + 1 - boundary.topleft[,2] 

########################

# bottom right
tmp = N - tmp
boundary.bottomright =
 cbind( c(N:(N2+1),
          tmp,
          tmp[N2]:N),
        c(rep(1, N2), 1:N2, rep(N2, N - tmp[N2] + 1)))

boundary.topright = boundary.bottomright
boundary.topright[,2] = N + 1 - boundary.topright[,2] 

###############################################

boundaries = list(boundary.bottomleft, boundary.topleft, boundary.bottomright, boundary.topright)

invisible(sapply(boundaries, polygon, lwd = 2))

##########################################################
# This is where we do the transformation

 big = addBoundary(dat, .25)
  # Tells us how much we have added to each dimension.
 added = as.integer((dim(big) - dim(dat))/2)

   # show the results.
if(plotInExpanded) {
image(1:nrow(big), 1:ncol(big), t(big), xlab = "X", ylab = "Y")

text(c(N/4, N/4, 3*N/4, 3*N/4, N2) + added[1],
     c(N/4, 3*N/4, N/4, 3*N/4, N2) + added[2],
     c(t(dat)[N2, N/4],
       t(dat)[N2, 3*N/4],
       t(dat)[3*N/4, N/4],
       t(dat)[3*N/4, 3*N/4],
       t(dat)[N2, N2]))
}

    # For a 1500 x 1500 grid, on a Macbook pro (4Gb RAM, 2.6 Ghz),
    # about 140 seconds. On a Linux box with 32Gb RAM, 4 dual core 2.4Ghz chip
    # 210 seconds!
 cart = cartogram(big)

    # We have to move the boundaries relative to the edge of big.
new.boundaries = lapply(boundaries,
                          function(x) {
                                    # So this is our indexing within the origina grid but moved to the
                                    # where we embedded x earlier in the larger "sea"
                               predict(cart, x[,1] + added[1], x[,2] + added[2])
                          })

if(plotInExpanded) {
 invisible(sapply(new.boundaries, polygon, lwd = 4, lty = 3, border = "blue"))
} else
 invisible(sapply(new.boundaries, function(x) polygon(x$x - added[1], x$y - added[2], lwd = 4, lty = 3, border = "blue")))

