cols = {'r' 'g' 'b' 'c' 'm' 'y'};
Xinit = {[-1.5;2.5] [1.7;2.7] [-.3;.85]};
niter = 10;
clf; hold on;
imagesc(x2, x1, perform_hist_eq(F, 'linear') );
colormap jet(256);
for k=1:length(Xinit)    
    x = Xinit{k}; % [ 4*rand-2; 3.5*rand-.5];
    X = [];
    for i=1:niter
        X(:,end+1) = x;
        x = x - pinv(Hessf(x))*Gradf(x);
    end
    plot(X(2,:), X(1,:), [cols{k} '.-'], 'MarkerSize', 20, 'LineWidth', 2);
end
axis([-.5 3 -2 2]);
