lambda_max = .3;
niter = 200;
fSpars = PhiS(y);
lambda_list_Ortho = linspace(lambda_max,0,niter);
errOrtho = [];
% Warmup
for i=1:niter
    fSpars = fSpars + tau * PhiS( y-Phi(fSpars) );
    fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(1)*tau );
end
% Iterations
for i=1:niter
    fSpars = fSpars + tau * PhiS( y-Phi(fSpars) );
    fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(i)*tau );
    % record the error
    errOrtho(i) = snr(f0,fSpars);
    if i>1 && errOrtho(i)>max(errOrtho(1:i-1))
        fBestOrtho = fSpars;
    end
end
clf;
h = plot(lambda_list_Ortho,errOrtho); 
axis tight;
set_label('lambda', 'SNR');
set(h, 'LineWidth', 2);
