niter = 300
lambda_list = linspace(1e-6,.01,20)
tau = 1.9 / ( 1 + maximum(lambda_list) * 8 / epsilon)
fBest = y
fTV = y
err = zeros(length(lambda_list), 1)
for it in  1 : length(lambda_list)
    # print('.',end="",flush=True)
    Lambda = lambda_list[it]
    for i in 1 : niter
        # Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon^2 + sum(Gr.^2, 3))
        G = -NtToolBox.div(Gr ./ repeat3(d) )
        # step
        e = Phi(fTV, h) - y
        fTV = fTV - tau.*( Phi(e, h) + Lambda.*G)
    end
    err[it] = snr(f0, fTV)
    if err[it] > snr(f0, fBest)
        fBest = fTV
    end
end
clf
plot(lambda_list, err)
axis("tight")
xlabel(L"\lambda")
ylabel("SNR")
fTV = fBest
