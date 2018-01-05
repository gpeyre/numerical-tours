q = 50;
lambda_list = linspace(0,20,q);
W = []; E = [];
for i=1:q
    lambda = lambda_list(i);
    w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
    W(:,i) = w; % bookkeeping
    E(i) = sqrt( sum( (X1*w-y1).^2 ) / n1 );
end
% Display error evolution.
clf;
plot(lambda_list, E, 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda');
ylabel('E');
