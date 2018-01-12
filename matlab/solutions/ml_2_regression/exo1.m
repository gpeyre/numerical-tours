q = 50;
lambda_list = linspace(0,.8,q);
W = []; E = [];
for i=1:q
    lambda = lambda_list(i);
    w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
    W(:,i) = w; % bookkeeping
    E(i) = sqrt( sum( (X1*w-y1).^2 ) / n1 ) / mean(y1.^2);
end
% find optimal lambda
[~,i] = min(E);
lambda0 = lambda_list(i);
% Display error evolution.
clf; hold on;
plot(lambda_list, E, 'LineWidth', 2);
plot( lambda0*[1 1], [min(E) max(E)], 'r--', 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda');
ylabel('E'); 
box on;
