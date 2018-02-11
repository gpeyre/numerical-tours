% fill matrix.
[I,J] = ind2sub( [n n], find(x1(:)~=0) ); v = x1(x1(:)~=0);
Ainp = [];
for i=1:length(I)
    Z = zeros(n,n,n);
    Z(I(i), J(i), v(i)) = 1; % double( ==k );
    Ainp(end+1,:) = Z(:)';
end
A = [Aenc; Arow; Acol; Ablock; Ainp];
% solve
niter = 6;
u = ones(n,n,n);
err = [];
for i=1:niter
    Xrw = solvel1( A*diag(u(:)) ) .* u;
    u = (abs(Xrw).^(1-alpha)+epsilon);
    % 
    [tmp,xrw] = min( abs(Xrw-1), [], 3 );
    Xrw1 = encode(xrw); 
    err(end+1) = sum(A*Xrw1(:)~=1);
end
% 
clf;
h = plot(err, '.-');
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', 20);
title('Number of invalidated constraints');
