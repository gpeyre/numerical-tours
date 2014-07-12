niter = 200;
u = f;
energy = [];
for i=1:niter    
    energy(i) = 1/2*norm( f-u,'fro' )^2 + lambda*J(u);
    u = u - tau*( u - f + lambda*GradJ(u) );
end
clf;
plot(1:niter, energy); axis('tight');
set_label('iteration', 'Energy');
