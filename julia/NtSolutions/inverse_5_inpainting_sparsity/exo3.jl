niter = 1000
a = U.*PsiS(fSpars)
E = []
for i in 1 : niter
    fTI = Psi(a)
    d = y - Phi(fTI, Omega)
    E = [E; 1/2.*vecnorm(d )^2 + lambd*sum(abs(a))]
    # step 
    a = SoftThresh(a + tau.*PsiS(Phi(d, Omega)), lambd*tau)
end

figure(figsize = (7, 5))    
plot(E)
show()