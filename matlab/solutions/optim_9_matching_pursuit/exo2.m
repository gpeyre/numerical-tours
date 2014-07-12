% --- MP --- %
M = 2*N;
x = zeros(N,1);
X_mp = []; lambda_mp = []; E_mp = [];
for k=1:M
    E_mp(k) = norm(y-Phi*x);
    c = Phi'*(y-Phi*x);
    [lambda_mp(k),i] = max(abs(c));
    x(i) = x(i) + c(i);
    % record
    X_mp(:,k) = x;
end
% --- OMP --- %
M = P;
x = zeros(N,1);
X_omp = []; lambda_omp = []; E_omp = [];
for k=1:M
    E_omp(k) = norm(y-Phi*x);
    c = Phi'*(y-Phi*x);
    [lambda_omp(k),i] = max(abs(c));
    x(i) = x(i) + c(i);
    % projection    
    I = find(x~=0);
    x(I) = pinv(Phi(:,I))*y;
    % record
    X_omp(:,k) = x;
end
% display
I = find(x0~=0); J = setdiff(1:N,I);
clf;
subplot(2,1,1); hold on;
h = plot(E_mp, X_mp(I,:)', '-'); set(h, 'LineWidth', 2);
h = plot(E_mp, X_mp(J,:)', 'k-'); set(h, 'LineWidth', 2);
title('MP'); axis tight;  box on;
subplot(2,1,2); hold on;
h = plot(E_omp, X_omp(I,:)', '-'); set(h, 'LineWidth', 2);
h = plot(E_omp, X_omp(J,:)', 'k-'); set(h, 'LineWidth', 2);
title('OMP'); axis tight; box on;
