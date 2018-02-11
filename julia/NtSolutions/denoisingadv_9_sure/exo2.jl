f = f0 + sigma.*rand(Normal(), n, n)
ntrials = 30
lambda_list = linspace(.3, 2, ntrials)
E0 = []
E = []

for i in 1:ntrials
    lambd = lambda_list[i]
    hf = h(f, lambd)
    E0 = [E0; vecnorm(hf - f0)^2]
    E  = [E; SURE(f, hf, lambd)]
end

figure(figsize = (7, 5))

plot(lambda_list, E0)
plot(lambda_list, E)
legend([L"$E_0$", L"$E$"])
xlabel(L"$\lambda$")
xlim(minimum(lambda_list), maximum(lambda_list))
ylim([minimum(E), 400])

show()
