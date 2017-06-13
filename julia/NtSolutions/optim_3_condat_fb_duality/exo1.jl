Lambda = 25
gamma = 1.9/8 # we must have 0 < gamma < 2/8
nbiter = 400
(N1,N2) = size(y)
u = zeros(N1,N2,2)
Ep_array = zeros(nbiter) # array for the primal energy E_p
Ed_array = zeros(nbiter) # array for the dual energy E_d
sqnormy = norm(vec(y))^2/2
x=zeros(size(y))
for iter in 1:400
    x = grad_f_conj(-Dadj(u),y)
    u = prox_g_conj(u + gamma*D(x), Lambda)
    Ep_array[iter] = norm(vec(x-y)).^2/2 + Lambda*sum(sqrt(sum(D(x).^2,3))) 
    Ed_array[iter] = norm(vec(y-Dadj(u)))^2/2 - sqnormy
end
xdenoised = x;