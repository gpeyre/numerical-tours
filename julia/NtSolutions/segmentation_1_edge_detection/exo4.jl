figure(figsize = (10, 10))
slist = [4, 6, 10, 15]

for i in 1:length(slist)
    sigma = slist[i]
    subplot(2, 2, i)
    NtToolBox.plot_levelset(delta(blur(f0, sigma)) , 0, f0)
    title(@sprintf("sigma = %i", sigma))
end
