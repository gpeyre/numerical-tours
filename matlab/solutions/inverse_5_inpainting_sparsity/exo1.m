fSpars = y;
energy = [];
niter = 1000;
for i=1:niter
    fSpars = SoftThreshPsi( ProjC(fSpars,Omega), lambda );
    % record the energy
    fW = PsiS(fSpars);
    energy(i) = 1/2 * norm(y-Phi(fSpars,Omega), 'fro')^2 + lambda * sum(abs(fW(:)));
end
clf;
h = plot(energy);
axis('tight');
set_label('Iteration', 'E');
if using_matlab()
    set(h, 'LineWidth', 2);
end
