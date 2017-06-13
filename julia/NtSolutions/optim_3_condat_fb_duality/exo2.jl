gamma = 1/8 # put the right value here
a = 3   # try different values
Lambda = 25
nbiter = 400
(N1,N2) = size(y)
u=zeros(N1,N2,2)
Ep_array_fista = zeros(nbiter) # array for the primal energy E_p
Ed_array_fista = zeros(nbiter) # array for the dual energy E_d
sqnormy = norm(vec(y))^2/2

z = u

for iter in 1:400

    uprevious = u
    u = prox_g_conj(z + gamma*D(grad_f_conj(-Dadj(z),y)), Lambda)
    alpha = (iter-1)/(iter+a)
    z = u + alpha * (u - uprevious)
    x = grad_f_conj(-Dadj(u),y)
    
    Ep_array_fista[iter] = norm(vec(x-y))^2/2 + Lambda*sum(sqrt(sum(D(x).^2,3))) 
    Ed_array_fista[iter] = norm(vec(y-Dadj(u)))^2/2 - sqnormy
end
xdenoised = x;