gamma = gamma1;
displist = round(linspace(1,niter,10));
k = 1;
clf; hold on;
for i=1:niter
    gamma = resample( gamma + dt * normalC(gamma) );
    if i==displist(k)
        % display
        h = plot(gamma([1:end 1]), 'r');
        if i==1 || i==niter
            set(h, 'LineWidth', 2);
        end
        k = k+1;
        drawnow;
        axis('tight');  axis('off');
    end
end
