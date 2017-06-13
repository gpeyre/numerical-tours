figure(figsize = (10,10))
glist = [.1,.01,.005,.001]
niter = 300

for k in 1:length(glist)
    gamma = glist[k]
    xi = exp(-C/gamma)
    b = ones(N[2])

    for i in 1:niter
        a = p./(xi*b)
        b = q./(xi'*a)
    end
    Pi = (diagm(a)*xi)*diagm(b)
    imageplot(clamp(Pi,0,minimum(1./Array(N))*.3),"\gamma=$gamma", [2,2,k])
end
