eta_list = linspace(.1, 1, 10)
ntrials = 20
mu_mean = []
mu_std = []
N = 500

for i in 1:length(eta_list)
    eta = eta_list[i]
    P = Int(round(eta*N))
    c = []

    for k in 1:ntrials
        c = [c; mu(PhiRand(P, N))]
    end

    mu_mean = [mu_mean; mean(c)]
#     mu_std = [mu_std; std(c)] # Tere is a problem with std in Julia

end

k_mean = 1/2.*(1 + 1./mu_mean)

figure(figsize = (10, 7))

subplot(2, 1, 1)
plot(eta_list, mu_mean, linewidth = 2)
xlabel(L"$\eta$")
ylabel(L"$\mu$")

subplot(2, 1, 2)
plot(log10(eta_list), log10(k_mean), linewidth = 2)
xlabel(L"$\log_{10}(\eta)$")
ylabel(L"$\log_{10}(1/ 2(1 + 1/ \mu))$")

show()
