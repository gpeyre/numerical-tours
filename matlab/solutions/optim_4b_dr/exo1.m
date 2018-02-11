lun = []; err = [];
tx = zeros(n,1);
for i=1:niter
    tx = (1-mu/2)*tx + mu/2*rproxG( rproxF(tx,y),gamma );
    x = proxF(tx,y);
    lun(i) = norm(x,1);
    err(i) = norm(y-A*x);
end
clf;
plot(lun);
axis tight;
