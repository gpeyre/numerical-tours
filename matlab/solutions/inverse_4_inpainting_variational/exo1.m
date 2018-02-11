niter = 50;
E = [];
k = 1; ndisp = [1 5 10 niter];
norm1 = @(f)norm(f(:));
f = y;
for i=1:niter
    E(i) = norm1(grad(f));
    f = Pi( f + tau*Delta(f) );
    if i==ndisp(k)
        imageplot(f, ['iter=' num2str(i)], 2, 2, k);
        k = k+1;
    end
end
