library(Rcartogram)

data(uspop)

out = cartogram(uspop)

if(FALSE) {
   g = expand.grid(1:ncol(pop), 1:nrow(pop))
   out = cartogram(pop, g[,1], g[,2])
}

