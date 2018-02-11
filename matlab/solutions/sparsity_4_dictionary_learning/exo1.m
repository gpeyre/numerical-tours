niter = 100;
gamma = 1.6/norm(D)^2;
E = [];
X = zeros(p,m);
for i=1:niter
    R = D*X-Y;
    E(end+1,:) = sum(R.^2);
    X = ProjX(X - gamma * D'*R, k);
end
sel = 1:5;
clf;
plot(log10(E(1:end/2,sel) - repmat(min(E(:,sel),[],1),[niter/2 1])  ));
axis tight;
title('log_{10}(J(x_j) - J(x_j^*))');
