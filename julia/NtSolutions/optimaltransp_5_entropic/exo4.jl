niter = 800
b = ones(N,N,K)
a = copy(b)
Err_q = zeros(niter)


for i in 1:niter

    for k in 1:K
        Err_q[i] = Err_q[i] + norm(a[:,:,k].*xi(b[:,:,k]) - P[:,:,k])/norm(P[:,:,k])
        a[:,:,k] = P[:,:,k]./xi(b[:,:,k])
    end

    q = zeros(N,N)

    for k in 1:K
        q = q + lambd[k] * log(max(1e-19, b[:,:,k].*xi(a[:,:,k])))

        # if sum(1e-19 .> b[:,:,k].*xi(a[:,:,k]))>0
        #     x = b[:,:,k].*xi(a[:,:,k])
        #     #print(x[1e-19.> b[:,:,k].*xi(a[:,:,k])])
        # end
        # if (sum(b[:,:,k].*xi(a[:,:,k]).==0)>0)
        #     println((b[:,:,k].*xi(a[:,:,k])))
        # end
        # println((sum(b[:,:,k].*xi(a[:,:,k]).<0)))
        # x = b[:,:,k].*xi(a[:,:,k])
        # x[!(x.>0)]=0
        # x = log(x)
        # x[x.==-Inf]=0
        #println(sum(x))
        # q = q + lambd[k] * x
    end
    #println(sum(q))
    q = exp(q)

    for k in 1:K
        b[:,:,k] = q./xi(a[:,:,k])
    end

end
figure(figsize=(7,5))
plot(log(Err_q),linewidth = 2)
xlim(0,niter)
