library(Rcartogram)

m = matrix(0, 100, 100)
g = expand.grid(20:80, 20:80)
m[g[,1], g[,2]] = rpois(nrow(g)^2, 7)
# Not quite what we want. Need to 
# take account that we are adding this many
# values so it will change the mean
# 
m[m == 0] = mean(m)

cartogram(m)
