niter_learning = 10;
niter_dico = 50;
niter_coef = 100;
E0 = [];
X = zeros(p,m);
D = D0;
for iter = 1:niter_learning
    % --- coefficient update ----
    E = [];
    gamma = 1.6/norm(D)^2;
    for i=1:niter
        R = D*X - Y;
        E(end+1,:) = sum(R.^2);
        X = ProjX(X - gamma * D'*R, k);
    end
    E0(end+1) = norm(Y-D*X, 'fro')^2;
    % --- dictionary update ----
    E = [];
    tau = 1/norm(X*X');
    for i=1:niter_dico
        R = D*X - Y;
        E(end+1) = sum(R(:).^2);
        D = ProjC( D - tau * (D*X - Y)*X' );
    end
    E0(end+1) = norm(Y-D*X, 'fro')^2;
end
clf; hold on;
plot(1:2*niter_learning, E0);
plot(1:2:2*niter_learning, E0(1:2:2*niter_learning), '*');
plot(2:2:2*niter_learning, E0(2:2:2*niter_learning), 'o');
axis tight;
legend('|Y-DX|^2', 'After coefficient update', 'After dictionary update');
