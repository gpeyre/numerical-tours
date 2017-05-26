fig=figure(figsize=(20,15))
X1 = Array(X)
err = [pnoisy]
tW = Array(tW)
for i in 1:12
    X1 = X1*tW'
    err = [err; [snr(X0, X1)]]

    if i%2 == 0
        plot_mesh(X1, F, sub=[2,3,Base.div(i,2)])
    end
end
