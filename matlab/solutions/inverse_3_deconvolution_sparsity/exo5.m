switch setting
    case 1
        lambda_max = 0.002;
    case 2
        lambda_max = .02;
    otherwise
        error('Unknown setting 5');
end
epsilon = 1e-2;
niter = 1000;
lambda_list_TV = linspace(lambda_max/4,lambda_max,niter);
tau = 1.9 / ( 1 + max(lambda_list_TV) * 8 / epsilon);
fBestTV = y; fTV = y;
errTV = []; E = [];
% warmup stage
for i=1:1000
    lambda = lambda_list_TV(1);    
    % Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );
    % step
    e = Phi(fTV)-y;
    fTV = fTV - tau*( Phi(e) + lambda*G);
end
% descent
for i=1:niter
    lambda = lambda_list_TV(i);    
    % Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );
    % step
    e = Phi(fTV)-y;
    fTV = fTV - tau*( Phi(e) + lambda*G);
    % record error
    errTV(i) = snr(f0,fTV);
    E(i) = 1/2*norm(e, 'fro')^2 + lambda * sum(d(:));
    if errTV(i)>snr(f0,fBestTV)
        fBestTV = fTV;
    end
end
% clf; plot(log10(E(1:end/2)/E(end)-1));
clf;
plot(lambda_list_TV,errTV);
axis('tight');
xlabel('\lambda'); ylabel('SNR');
