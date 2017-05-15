figure(figsize = (10, 10))
slist = [1, 2, 4, 6]

for i in 1:length(slist)
    sigma = slist[i]
    d = sqrt(sum(nabla(blur(f0, sigma)).^2, 3))
    t = maximum(d)*1./5
    imageplot((d .> t)[:, :, 1], @sprintf("sigma = %.1f", sigma) , [2, 2, i])
end
