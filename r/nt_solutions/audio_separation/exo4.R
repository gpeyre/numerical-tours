nbins = 100
t = seq(from=pi/200, to = pi, by=(pi-pi/200)/(200-1))
s1 = c(seq(from=2, to=nbins),nbins-1)-1
s2 = c(2, seq(from=1, to=nbins-1))-1
I = which((h[s1+1]<h) & (h[s2+1]<h))
v = h[I][order(h[I])]
u = order(h[I], decreasing=TRUE)
theta1 = t[I[u[1:3]]]
M1 = rbind(cos(theta1), sin(theta1))
print(M, main="M")
print(M1, main="M1")