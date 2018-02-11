niter = 1000;
a = U.*PsiS(fSpars);
E = [];
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI,Omega);
    E(i) = 1/2*norm(d , 'fro')^2 + lambda * sum( abs(a(:)) );   
    % step 
    a = SoftThresh( a + tau*PsiS(Phi(d,Omega)), lambda*tau );
end
clf;
plot(E); axis('tight');
