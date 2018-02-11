niter = 500;
a = PsiS(y);
E = [];
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI);
    % measure energy
    E(end+1) = 1/2*norm(d, 'fro')^2 + lambda*sum(sum(sum( abs(a.*U) )));
    % step
    a = SoftThresh( a + tau*PsiS(Phi(d)), lambda*tau );
end
clf;
plot(E); axis('tight');
