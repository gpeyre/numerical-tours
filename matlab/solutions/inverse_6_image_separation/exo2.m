niter = 500;
u = f;
E = [];
for i=1:niter
    u = u - tau * ( GradT(u-f) + lambda*GradJ(u) );
    E(end+1) = T(f-u) + lambda*J(u);
end
clf;
plot(E); 
axis('tight');
