fSpars = y;
energy = [];
niter = 100;
for i=1:niter
    fSpars = fSpars + tau * Phi( y-Phi(fSpars) );
    % thresholding
    fSpars = SoftThreshPsi( fSpars, lambda*tau );
    % record the energy
    fW = PsiS(fSpars);
    energy(end+1) = 1/2 * norm(y-Phi(fSpars), 'fro')^2 + lambda * sum(abs(fW(:)));
end
clf;
h = plot(energy); axis([1,niter,min(energy)*1.05+max(energy)*(-.05) max(energy)]);
set_label('Iteration', 'E');
set(h, 'LineWidth', 2);
