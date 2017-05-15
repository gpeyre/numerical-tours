d = sum(P.^2, 2)[:]
rho = .1
v = sort(d)
I = reverse(sortperm(d))

#transformed points
I = I[collect(1:Int(round(rho*length(I))+1))]
P1 = P[I,:]

#compute Theta
nrow = size(P1)[1]
Theta = zeros(nrow)
for i in 1:nrow
    Theta[i] = mod(atan2(P1[i,2],P1[i,1]),pi)
end

nbins = 200
h,t = plt[:hist](Theta,nbins)
h=h/sum(h)
clf()
bar(t[1:end-1], h, width = pi/nbins)
xlim(0,pi);
