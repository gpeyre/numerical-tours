v = zeros(n,n,2);
niter = 100;
energy = [];
for i=1:niter
    U = Dt + sum(v.*D,3);
    L = cat(4,grad(v(:,:,1)),grad(v(:,:,2)));
    G = D.*repmat(U, [1 1 2])  - lambda * cat(3, div(L(:,:,:,1)), div(L(:,:,:,2)) );
    v = v - tau*G;
    energy(end+1) = sum(U(:).^2) + lambda*sum(L(:).^2);
end
clf;
plot(energy, '.-'); axis('tight');
xlabel('Iteration');
ylabel('Energy');
