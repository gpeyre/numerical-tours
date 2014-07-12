Xproj = encode(zeros(n));
niter = 50000;
err = [];
for i=1:niter
    err(end+1) = norm(A*Xproj(:)-1,'fro');
	Xproj = max(projector(Xproj), 0);
end
clf;
plot(log10(err/err(1))); 
axis('tight');
