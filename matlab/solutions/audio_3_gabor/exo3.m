niter = 100;
lambda_max = 1.4*sigma;
lambda_min = .8*sigma;
lambda_list = linspace(lambda_max, lambda_min, niter);
err = [];
for i=1:niter
    % progressbar(i,niter);
    lambda = lambda_list(i);
    % gradient
    r = x - x1;
    Sr = perform_stft(r, wlist, qlist, options);
    S1 = cell_add(S1, Sr);
    % threshold
    S1 = perform_thresholding(S1, lambda, 'soft');
    % update the value of lambda to match noise
    x1 = perform_stft(S1,wlist,qlist, options);
    % lambda = lambda * sqrt(n)*sigma / norm( x-x1, 'fro' ); 
    err(i) = snr(x0,x1);
    if i>1 & err(i)>max(err(1:i-1))
        xbp = x1; Sbp = S1;
    end
end
clf;
plot(lambda_list/sigma, err); axis('tight');
set_label('\lambda/\sigma', 'SNR');
