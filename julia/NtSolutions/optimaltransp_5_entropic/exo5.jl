# Insert your code here.

m = 5
T,S = meshgrid(collect(linspace(0,1,m)), collect(linspace(0,1,m)))
T = vec(T')
S = vec(S')
niter = 800

for j in 1:m^2
    # weights
    lambd = [S[j]*T[j] (1-S[j])*T[j] S[j]*(1-T[j]) (1-S[j])*(1-T[j])]
    # computation
    b = ones(N,N,K)
    a = copy(b)

    for i in 1:niter

        for k in 1:K
            a[:,:,k] = P[:,:,k]./xi(b[:,:,k])
        end

        q = zeros(N,N)

        for k in 1:K
            q = q + lambd[k] * log(max(1e-15, b[:,:,k].*xi(a[:,:,k])))
        end

        q = exp(q)

        for k in 1:K
            b[:,:,k] = q./xi(a[:,:,k])
        end
    end
    # display
    subplot(m,m,j)
    imageplot(q)
end
