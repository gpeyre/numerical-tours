Tmax = 100;
niter = ceil(Tmax/tau);
X = X0;
k = 1; 
displist = round([.05 .1 .5 1]*niter); 
clf;
for i=1:niter
    % step
    X = X - tau*X*tL';    
    if displist(k)==i
        subplot(2,2,k);
        plot_mesh(X,F);
        k = k+1;
    end
end
