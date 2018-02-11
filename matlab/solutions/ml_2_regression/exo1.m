q = 50;
lmax = norm(X0)^2;
lambda_list = lmax*linspace(.3,1e-3,q);
W = []; E = [];
for i=1:q
    lambda = lambda_list(i);
    w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
    W(:,i) = w; % bookkeeping
    E(i) = norm(X1*w-y1) / norm(y1);
end
% find optimal lambda
[~,i] = min(E);
lambda0 = lambda_list(i);
wRidge = W(:,i);
% Display error evolution.
clf; hold on;
plot(lambda_list/lmax, E, 'LineWidth', 2);
plot( lambda0/lmax*[1 1], [min(E) max(E)], 'r--', 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda/|X|^2');
ylabel('E'); 
box on;
fprintf('Ridge: %.2f%%\n', 100*min(E));
