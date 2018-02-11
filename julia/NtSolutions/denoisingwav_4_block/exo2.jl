figure(figsize = (10, 10))
tlist = linspace(.3, .9, 4)

for i in 1:length(tlist)
    T = tlist[i]
    imageplot(clamP(ThreshBlock(f, T)), @sprintf("T = %.1f", T), [2, 2, i])
end
