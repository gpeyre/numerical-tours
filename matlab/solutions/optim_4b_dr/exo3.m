lun = [];
tx = x0;
niter = 2000;
for i=1:niter
    % progressbar(i,niter);
    tx = (1-mu/2)*tx + mu/2*rproxG( rproxF(tx,y),gamma );
    x = proxF(tx,y);
    lun(i) = norm(x,1);
end
proba = [];
E = mean(abs(x-x0))<.05;
for j=1:length(slist)
    s = slist(j);
    proba(j) = mean(E(Slist==s));
end
clf;
h = plot(slist, proba, 'k');
if using_matlab()
    set(h, 'LineWidth', 2);
end
axis([min(slist) max(slist) -.03 1.03]);
