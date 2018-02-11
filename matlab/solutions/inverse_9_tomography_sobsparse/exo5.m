lambda_max = .3;
niter = 200;
fSpars = PhiS(y);
lambda_list_TI = linspace(lambda_max,0,niter);
errTI = [];
% Warmup
for i=1:niter
    fSpars = fSpars + tau * PhiS( y-Phi(fSpars) );
    fSpars = SoftThreshPsi(fSpars, lambda_list_TI(1)*tau );
end
% Iterations
for i=1:niter
    fSpars = fSpars + tau * PhiS( y-Phi(fSpars) );
    fSpars = SoftThreshPsi(fSpars, lambda_list_TI(i)*tau );
    % record the error
    errTI(i) = snr(f0,fSpars);
    if i>1 && errTI(i)>max(errTI(1:i-1))
        fBestTI = fSpars;
    end
end
clf;
h = plot(lambda_list_TI,errTI); 
axis tight;
set_label('lambda', 'SNR');
set(h, 'LineWidth', 2);
