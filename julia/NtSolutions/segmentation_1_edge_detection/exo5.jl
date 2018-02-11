figure(figsize = (10, 10))
slist = [4, 6, 10, 15]

for i in 1:length(slist)
    sigma = slist[i]

    g = grad(blur(f0, sigma))
    H = hessian(blur(f0, sigma))
    a = H[:, :, 1:2] .* repeat(g[:, :, 1], inner = [1, 1, 2]) + H[:, :, 2:3] .* repeat(g[:, :, 2], inner = [1, 1, 2])

    subplot(2, 2, i)
    plot_levelset(sum(a.*g, 3)[:, :, 1], 0, f0)
    title(@sprintf("sigma = %i", sigma))
end
