niter = 50;
s = [];
Y = X/q;
ndisp = [1 5 niter Inf];
clf; k =1; 
hold on;
Y = Y-repmat(mean(Y,2), [1 N]);
h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1]));
for i=1:niter
    Y = Y * B(D(Y))' / N;  
    % update
    Y = Y-repmat(mean(Y,2), [1 N]);
    % record stress
    s(end+1) = Stress(D(Y));
    if ndisp(k)==i
        plot(Y(2,[1:N0 1]), Y(1,[1:N0 1])); 
        axis('equal'); axis('off');
        k = k+1;
    end
end
axis('equal'); axis('off'); axis('ij');
