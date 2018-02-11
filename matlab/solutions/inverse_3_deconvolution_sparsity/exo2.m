switch setting
    case 1
        lambda_max = .01/4;
    case 2
        lambda_max = .02;
    otherwise
        error('Unknown setting 2');
end
niter = 1000;
fSpars = y;
lambda_list_Ortho = linspace(lambda_max,0,niter);
errOrtho = [];
% warmup
for i=1:niter
    fSpars = fSpars + tau * Phi( y-Phi(fSpars) );
    fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(1)*tau );
end
% descent
for i=1:niter
    % descent
    fSpars = fSpars + tau * Phi( y-Phi(fSpars) );
    % thresholding
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
