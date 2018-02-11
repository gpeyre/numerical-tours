b = ones(N)
niter = 2000
Err_p = []
Err_q = []

for i in 1:niter
    a = p./(xi*b)
    append!(Err_q, norm(b.*(xi*a) - q)/norm(q))
    b = q./(xi'*a)
    append!(Err_p, norm(a.*(xi'*b) - p)/norm(p))
end

figure(figsize = (10,7))

subplot(2,1,1)
title("||\pi -p||")
plot(log.(Array(Err_p)), linewidth = 2)
xlim(0,niter)

subplot(2,1,2)
title("||\pi^T -q||")
plot(log.(Array(Err_q)), linewidth = 2)
xlim(0,niter)
