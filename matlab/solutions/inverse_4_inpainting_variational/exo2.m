niter = 700;
J = [];
k = 1; ndisp = ceil(niter*[.001 1/10 1/5 1]);
f = y;
for i=1:niter
    J(i) = sum(sum(Amplitude(grad(f))));
    f = Pi( f - tau*G(f) );
    if i==ndisp(k)
        imageplot(f, ['iter=' num2str(i)], 2, 2, k);
        k = k+1;
    end
end
