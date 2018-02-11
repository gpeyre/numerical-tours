figure(figsize = (10, 10))
t_list = maximum(d).*[1./4, 1./5, 1./10, 1./20]

for i in 1:length(t_list)
    t = t_list[i]
    imageplot((d .> t)[:, :, 1], @sprintf("t = %.1f", t) , [2, 2, i])
end
