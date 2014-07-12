lambda_list = linspace(0,.3,100);
for i=1:length(lambda_list)
    lambda = lambda_list(i);
    err(i) = snr(M0, (1-lambda)*Mnl0 + lambda*M);
end
clf;
plot(lambda_list, err); axis('tight');
set_label('lambda', 'SNR');
[tmp,I] = max(err);
lambda = lambda_list(I);
Mnl1 = (1-lambda)*Mnl0 + lambda*M;
