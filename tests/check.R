library(Rcartogram)

if(file.exists("~/cart-1.2.2/output.dat")) {
p = read.table("~/cart-1.2.2/output.dat")
x = matrix(p[[1]], 513, 1025, byrow = TRUE)
y = matrix(p[[2]], 513, 1025, byrow = TRUE)

data(uspop)
out = cartogram(uspop)

max(abs(out$x - x)/(out$x + .00000001))
max(abs(out$y - y)/(out$y + .00000001))
}

