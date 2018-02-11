lambda_list = linspace(2*1e-3, .7, 40)
err = zeros( (length(lambda_list), 1) )
for i in 1 : length(lambda_list)
    Lambda = lambda_list[i]
    fSob = real( plan_ifft( yF .* hF ./ ( abs(hF).^2 + Lambda.*S) )*(yF .* hF ./ ( abs(hF).^2 + Lambda.*S)) )
    err[i] = snr(f0, fSob)
end
clf;
plot(lambda_list, err)
axis("tight")
xlabel("lambda")
ylabel("SNR")

val, i = findmax(err);
Lambda = lambda_list[i];

fSob = real( plan_ifft( yF .* hF ./ ( abs(hF).^2 + Lambda.*S) )*(yF .* hF ./ ( abs(hF).^2 + Lambda.*S)) )
