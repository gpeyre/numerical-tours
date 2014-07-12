f = f0 + sigma*randn(n);
ntrials = 30;
lambda_list = linspace(.3, 2, ntrials);
E0 = []; E = [];
for i=1:ntrials
    lambda = lambda_list(i);
    hf = h(f,lambda);
    E0(i) = norm(hf-f0, 'fro')^2;
    E(i) = SURE(f,hf,lambda);
end
clf;
plot(lambda_list, [E0; E]);
xlabel('\lambda'); legend('E_0', 'E');
axis tight;
