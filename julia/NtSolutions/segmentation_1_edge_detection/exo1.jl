figure(figsize = (10, 10))
slist = [1, 2, 5, 10]

for i in 1:length(slist)
    sigma = slist[i]
    imageplot(blur(f0, sigma), @sprintf("sigma = %i", sigma), [2, 2, i])
end
