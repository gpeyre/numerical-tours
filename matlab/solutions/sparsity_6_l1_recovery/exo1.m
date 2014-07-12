eta_list = linspace(.1,1,10);
ntrials = 20;
mu_mean = [];
mu_std = [];
N = 500;
for i=1:length(eta_list)
    eta = eta_list(i);
    P = round(eta*N);
    c = []; 
    for k=1:ntrials
        c(end+1) = mu(PhiRand(P,N));
    end
    mu_mean(end+1) = mean(c);
    mu_std(end+1) = std(c);
end
k_mean = 1/2*(1+1./mu_mean);
clf; 
subplot(2,1,1);
h = plot(eta_list, mu_mean); set(h, 'LineWidth', 2);
axis tight;
xlabel('\eta'); ylabel('\mu');
subplot(2,1,2);
%h = plot(eta_list, k_mean ); set(h, 'LineWidth', 2);
%xlabel('\eta'); ylabel('1/2(1+1/E(\mu))');
h = plot(log10(eta_list), log10(k_mean) ); set(h, 'LineWidth', 2);
xlabel('log_{10}(\eta)'); ylabel('log_{10}(1/2(1+1/\mu))');
axis tight;
