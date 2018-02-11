niter = 50;
stress = [];
Y = V;
ndisp = [1 5 10 niter Inf];
s = [];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k); 
        plot_mesh(Y,F, options);
        axis('equal'); axis('off');
        k = k+1;
    end
    Y = Y * B(D(Y))' / N;  
    % update
    % Y = Y-repmat(mean(Y,2), [1 N]);
    % record stress
    s(end+1) = Stress(D(Y));
end
axis('equal'); axis('off'); axis('ij');
