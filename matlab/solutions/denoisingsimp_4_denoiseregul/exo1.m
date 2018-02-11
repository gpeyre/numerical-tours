lambda_list = linspace(1,30,60);
err = [];
for i=1:length(lambda_list)
    lambda = lambda_list(i);
    fSob = real( ifft2( yF ./ ( 1 + lambda*S) ) );
    err(i) = snr(f0,fSob);
end
[tmp,i] = max(err);
lambda = lambda_list(i);
fSob0 = real( ifft2( yF ./ ( 1 + lambda*S) ) );
% plot
clf;
plot(lambda_list, err); axis('tight');
set_label('lambda', 'SNR');
