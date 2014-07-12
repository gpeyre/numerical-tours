X1 = X;
clf;
err = [pnoisy];
for i=1:12
    X1 = X1*tW';  
    err(i+1) = snr(X0,X1);
    if mod(i,2)==0
        subplot(2,3,i/2);
        plot_mesh(X1,F, options); 
        axis('tight'); shading('interp'); 
    end    
    if err(length(err))>max(err(1:length(err)-1))
        Xbest = X1;
    end
end
