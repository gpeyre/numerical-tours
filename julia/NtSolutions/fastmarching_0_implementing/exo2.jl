method = "fm"
svg_rate = n*6
(D,Dsvg,Ssvg) = perform_dijkstra_fm(W, x0, method=method, svg_rate=svg_rate)
for i in 1:4
    subplot(2,2,i)
    d = Dsvg[:,:,i]
    d[d.==Inf] = 0
    imageplot(d)
    set_cmap("jet")
end
