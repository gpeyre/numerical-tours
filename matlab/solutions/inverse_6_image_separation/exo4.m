lambda = .05;
tau = 1.9 /( max(W(:))^2 + 8*lambda/epsilon );
%
niter = 500;
E = [];
u = f;
for i=1:niter
    u = u - tau * ( GradT(u-f) + lambda*GradJ(u) );
    E(end+1) = T(f-u) + lambda*J(u);
end
clf;
plot(E); 
axis('tight');
