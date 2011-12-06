makeGrid =
  #
  # This takes a collection of population densities
  # and a box/rectangle (or a grid of x's and y's)
  # and creates a matrix with the population densities
  # spread appropriately through the box.
  # It uses the maps package to do this to determine the
  # polygons/boundaries for each 
  #
  # For each (x, y) point on the "grid", it figures out
  # which polygon the point is in and then determine the
  # population density.
  #
function(x, y, densities, dim = c(100, 100), land.mean = mean(densities))
{
  library(maps)

    # if x has 4 elements and we have no y,
    # break x up into  x and y
  if(missing(y) && length(x) == 4) {
     y = x[3:4]
     x = x[1:2]
  }

    # If x has two elements, then compute the
    # grid of x's, equally spaced
    # So we end up with a vector of x values and a vector
    # of y's and this gives us the (x, y) points for our grid.
  if(length(x) == 2)
     x = seq(min(x), max(x), length = dim[1])
  if(length(y) == 2)
     y = seq(min(y), max(y), length = dim[2])

    # Compute all the (x,y) points in our matrix.
  grid = expand.grid(x, y)

     # now figure out the name of the state in which
     # each of the x,y values is located.
  ids = map.where('state', grid)
  ids = gsub(":.*", "", ids) 

     
  names(densities) = tolower(names(densities))

   # Construct our matrix by getting the population densities
   # at each (x,y) pair based on what state it was in and looking
   # this up in 
  values = matrix(densities[ids], length(x), length(y), byrow = TRUE)

   # Where there are missing values, use the mean value.
  values[is.na(values)] = land.mean
  structure(list(grid = t(values), x = range(x), y = range(y),
                 mean = mean(densities)),
             class = "CartogramGrid")
#  t(values[nrow(values):1, ])
}

plot.CartogramGrid =
  #
  # Should we transpose the grid or not.
  #
function(x, ...)
  image(seq(x$x[1], x$x[2], length = ncol(x$grid)),
        seq(x$y[1], x$y[2], length = nrow(x$grid)),        
        t(x$grid), xlab = "latitude", ylab = "longitude")




statePopulations =
  #
  # Get the US state populations by reading them from Wikipedia.
  #
  #
function()
{
  pop = htmlParse("http://en.wikipedia.org/wiki/List_of_U.S._states_by_population", error = function(...){})
  els = pop["//table//tr/child::*[3] | //table//tr/child::*[4]" ]
  vals = sapply(els[ - (1:2) ], xmlValue)
  i = seq(1, by = 2, length = length(vals)/2)
  structure(as.numeric(gsub("^&0*([0-9]+)\\..*", "\\1", vals[i+1])), names = vals[i])
}

if(exists("run") && run) {
    # Read the state populations.
  gridSize = 400
  pop = statePopulations()
  library(maps)
  mm = map('usa', plot = FALSE)
  g = makeGrid(mm$range, , pop, c(gridSize, gridSize))

  library(Rcartogram)

  big = addBoundary(g$grid, 2, g$mean)
  image(1:nrow(big), 1:ncol(big), big)
  cart = cartogram(t(big))

  load("states.rda")
  new.polys = lapply(countyBoundaries$us,
                      function(pol) {
                        z = Rcartogram:::mapToGrid(pol[,1], pol[,2], g, c(800, 800))
                        z.new = predict(cart, z)
                      })

  tmp = do.call("rbind", lapply(new.polys, function(x) cbind(x$x, x$y)))
#  quartz()

  plot(0, type = "n", xlim = range(tmp[,1], na.rm = TRUE), ylim = range(tmp[,2], na.rm = TRUE))
  invisible(lapply(new.polys, polygon, border = "blue"))

  plot(0, type = "n", xlim = range(tmp[,1], na.rm = TRUE), ylim = range(tmp[,2], na.rm = TRUE))
  state.win = sapply(states, function(x) diff(colSums(x[,1:2]))  < 0)
  names(state.win) = gsub("\\..*$", "", names(state.win))
  names(state.win) = gsub("-", " ", names(state.win))

  i = match(tolower(names(new.polys)), names(state.win), 0)

  invisible(mapply(polygon, new.polys, ifelse(state.win[tolower(names(new.polys[i]))], "blue", "red")))

  invisible(polygon(do.call("rbind", lapply(new.polys[i], function(x) cbind(x$x, x$y))),
                    border = "blue",
                    ))
}


