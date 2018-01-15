q = 200;
lmax = max(abs(X0'*y0));
lambda_list = lmax*linspace(.6,1e-3,q);
W = []; E = [];
w = zeros(p,1);
niter = 500;
for iq=1:q
    lambda = lambda_list(iq);
    % ISTA %
    for i=1:niter
        w = ISTA(w,lambda,tau);
    end
    W(:,iq) = w; % bookkeeping
    E(iq) = norm(X1*w-y1) / norm(y1);
end
% find optimal lambda
[~,i] = min(E);
lambda0 = lambda_list(i);
wSparse = W(:,i);
fprintf('Lasso: %.2f%%\n', 100*min(E));
% Display error evolution.
Il = find(lambda_list);
clf; hold on;
plot(lambda_list/lmax, E, 'LineWidth', 2);
plot( lambda0/lmax*[1 1], [min(E) max(E)], 'r--', 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda/|X^* y|_\infty');
ylabel('E'); box on;
