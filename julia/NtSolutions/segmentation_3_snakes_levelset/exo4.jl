n = 200
Y, X = meshgrid(collect(1 : n), collect(1 : n))
r = n/3.
c = [n, n]./2
phi0 = max(abs(X - c[1]), abs(Y - c[2])) - r
