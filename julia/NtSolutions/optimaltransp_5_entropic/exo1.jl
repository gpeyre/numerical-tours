b = ones(N[2])
niter = 300
Err_p = []
Err_q = []

for i in 2:niter
    a = p./(xi*b)
    append!(Err_q, norm(b.*(xi'*a) - q)/norm(q))
    b = q ./(xi'*a)
    append!(Err_p, norm(a.*(xi*b) - p)/norm(p))
end

figure(figsize = (10,7))

subplot(2,1,1)
title("||\pi -p||")
plot(log(Array(Err_p) + 1e-5), linewidth = 2)
xlim(0,niter)

subplot(2,1,2)
title("||\pi^T -q||")
plot(log(Array(Err_q) + 1e-5), linewidth = 2)
xlim(0,niter)
