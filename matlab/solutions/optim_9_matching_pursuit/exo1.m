M = 2*N;
x = zeros(N,1);
X = []; lambda = []; E = [];
for k=1:M
    E(k) = norm(y-Phi*x);
    c = Phi'*(y-Phi*x);
    [lambda(k),i] = max(abs(c));
    x(i) = x(i) + c(i);
    % record
    X(:,k) = x;
end
clf;
h = plot(log10(E)); axis tight;
set(h, 'LineWidth', 2);
title('log(E)');
