Xt = X;
k = 0; sob = []; err = [];
clf;
for i=1:niter
    % step
    Xt = Xt - tau*Xt*tL';   
    % error
    err(i) = snr(X0,Xt);
    if mod(i,floor(niter/4))==0
        k = k+1;
        subplot(2,2,k);
        plot_mesh(Xt,F, options);
        shading('interp'); axis('tight');
        % title(strcat(['T=' num2str(Tmax*k/4,3)]));
    end
end
