niter = 400;
lambda_list = linspace(1e-6,.01,20);
tau = 1.9 / ( 1 + max(lambda_list) * 8 / epsilon);
fBest = y; fTV = y;
err = [];
for it=1:length(lambda_list)
    lambda = lambda_list(it);
    for i=1:niter
        % Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV);
        d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
        G = -div( Gr./repmat(d, [1 1 2])  );
        % step
        e = Phi(fTV,h)-y;
        fTV = fTV - tau*( Phi(e,h) + lambda*G);
    end
    err(it) = snr(f0,fTV);
    if err(it)>snr(f0,fBest)
        fBest = fTV;
    end
end
clf;
plot(lambda_list,err);
axis('tight');
xlabel('\lambda'); ylabel('SNR');
fTV = fBest;
