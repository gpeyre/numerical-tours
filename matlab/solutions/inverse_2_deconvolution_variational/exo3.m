tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
E = []; 
for i=1:niter
    % Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );
    % step
    e = Phi(fTV,h)-y;
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
    % energy
    E(i) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
end
% display energy
clf;
plot(E); axis('tight');
set_label('Iteration 3', 'Energy');
