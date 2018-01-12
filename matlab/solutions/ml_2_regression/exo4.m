q = 200;
lmax = max(X0'*y0);
lmax = 3.5;
lambda_list = linspace(lmax,1e-3,q);
W = []; E = [];
w = zeros(p,1);
niter = 1000;
for iq=1:q
    lambda = lambda_list(iq);
    % ISTA %
    for i=1:niter
        w = Soft( w-tau*X0'*(X0*w-y0), lambda*tau );
    end
    W(:,iq) = w; % bookkeeping
    E(iq) = sqrt( sum( (X1*w-y1).^2 ) / n1 ) / mean(y1.^2);
end
% find optimal lambda
[~,i] = min(E);
lambda0 = lambda_list(i);
% Display error evolution.
Il = find(lambda_list<2);
clf; hold on;
plot(lambda_list(Il), E(Il), 'LineWidth', 2);
plot( lambda0*[1 1], [min(E(Il)) max(E(Il))], 'r--', 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda');
ylabel('E'); box on;
