niter = 200;
fTV = y;
energy = [];
for i=1:niter
    Gr = grad(fTV);
    d = sqrt(sum3(Gr.^2,3));
    deps = sqrt( epsilon^2 + d.^2 );
    G0 = -div( Gr ./ repmat(deps, [1 1 2]) );
    G = fTV-y+lambda*G0;
    energy(i) = 1/2*norm( y-fTV,'fro' )^2 + lambda*sum(deps(:));
    fTV = fTV - tau*G;
end
clf;
plot(1:niter, energy); axis('tight');
set_label('iteration', 'Energy');
