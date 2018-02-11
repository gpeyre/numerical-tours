lambda_list = linspace(.01, 12, 40);
err = [];
for i=1:length(lambda_list)
    lambda = lambda_list(i);
    fSob = real( ifft2( y .* xi ./ ( xi + lambda*S) ) )*n;
    err(i) = snr(f0,fSob);
end
clf;
h = plot(lambda_list, err); axis tight;
set_label('lambda', 'SNR');
set(h, 'LineWidth', 2);
[tmp,i] = max(err);
lambda = lambda_list(i);
fSob = real( ifft2( y .* xi ./ ( xi + lambda*S) ) )*n;
