#
a = 5; b = 3.5
tx = linspace(-a,a,q); ty = linspace(-b,b,q)
B,A = meshgrid( ty,tx )
G = [A[:] B[:]]
#
offs = [0.3 1 3 5]
niter = 10000
for io=1:length(offs)
    # generate data
    omega = offs[io]*[1 .5]
    X = [randn(div(n,2),2)-ones(div(n,2),1)*omega/2; randn(div(n,2),2)+ones(div(n,2),1)*omega/2]
    X = [randn(div(n,2),2)-ones(div(n,2),1)*omega/2; randn(div(n,2),2)+ones(div(n,2),1)*omega/2]
    # run gradient descent
    w = zeros(p+1,1)
    for i=1:niter
        w = w - tau * nablaE(w,AddBias(X),y)
    end
    # display
    Theta = reshape(theta(AddBias(G)*w), q, q)
    subplot(2,2,io)
    imshow(Theta'[:,end:-1:1], extent=[-a, a, -b, b]);
    plot_multiclasses(X,y,disp_legend=0, ms=5)
    axis([-a, a, -b, b])
    axis("off")
end
