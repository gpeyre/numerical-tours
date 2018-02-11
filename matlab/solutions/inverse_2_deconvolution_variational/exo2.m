lambda_list = linspace(.03, .2, 40);
err = [];
for i=1:length(lambda_list)
    lambda = lambda_list(i);
    fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
    err(i) = snr(f0,fSob);
end
clf;
hh = plot(lambda_list, err); axis tight;
set_label('lambda', 'SNR');
if using_matlab()
    set(hh, 'LineWidth', 2);
end
[tmp,i] = max(err);
lambda = lambda_list(i);
fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
