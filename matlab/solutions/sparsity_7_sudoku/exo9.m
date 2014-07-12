% fill the constraint matrix
[I,J] = ind2sub( [n n], find(x1(:)~=0) ); v = x1(x1(:)~=0);
Ainp = [];
for i=1:length(I)
    Z = zeros(n,n,n);
    Z(I(i), J(i), v(i)) = 1; % double( ==k );
    Ainp(end+1,:) = Z(:)';
end
A = [Aenc; Arow; Acol; Ablock; Ainp];
pA = pinv(A);
projector = @(u)reshape( u(:) - pA*(A*u(:)-1), [n n n]);
% POCS
Xproj = encode(zeros(n));
err = [];
for i=1:niter
    err(end+1) = norm(A*Xproj(:)-1,'fro');
	Xproj = clamp(projector(Xproj),0,1);
end
clf;
plot(log10(err/err(1))); 
axis('tight');
% Check wether this is a valid solution.
[tmp,xproj] = min( abs(Xproj-1), [], 3 );
Xproj1 = encode(xproj); 
disp(['Number of violated constraints: ' num2str(sum(A*Xproj1(:)~=1)) '.']);
