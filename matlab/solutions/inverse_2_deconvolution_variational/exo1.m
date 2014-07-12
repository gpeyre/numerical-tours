eps_list = linspace(2*1e-3, .03, 40);
err = [];
for i=1:length(eps_list)
    lambda = eps_list(i);
    fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
    err(i) = snr(f0,fL2);
end
clf;
hh = plot(eps_list, err); axis tight;
set_label('lambda', 'SNR');
if using_matlab()
    set(hh, 'LineWidth', 2);
end
[tmp,i] = max(err);
lambda = eps_list(i);
fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
