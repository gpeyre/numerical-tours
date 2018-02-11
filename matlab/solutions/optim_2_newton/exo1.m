niter = 10;
x = [-1.5;2.5];
x = [1.7;2.7];
E = []; X = [];
for i=1:niter
    X(:,end+1) = x;
    x = x - pinv(Hessf(x))*Gradf(x);
    E(i) = f(x(1),x(2));
end
% display
myplot = @(E)plot(log10(E(E>eps)));
Xs = [1;1];
e = sqrt( sum( (X-repmat(Xs, [1 niter])).^2) );
clf;
subplot(2,1,1);
myplot(E);  axis tight;
title('log_{10}|E(x_k)-E^*|');
subplot(2,1,2);
myplot(e);  axis tight;
title('log_{10}|x_k-x^*|');
