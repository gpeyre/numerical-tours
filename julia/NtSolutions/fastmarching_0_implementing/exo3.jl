method = "fm"
(D,Dsvg,Ssvg) = perform_dijkstra_fm(W, x0, method=method)
k = 8
displ = D -> cos(2*pi*k*D/maximum(D))
imageplot(displ(D))
set_cmap("jet")
