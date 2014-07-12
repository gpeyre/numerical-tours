GradT = @(f)real(ifft2(W.^2.*fft2(f)));
T = @(v)1/2*norm( W.*fft2(v)/n, 'fro' ).^2;
clf;
imageplot(GradT(f), 'Grad(T)', 1,2,1);
imageplot( f-GradT(f), 'f-Grad(T)', 1,2,2);
