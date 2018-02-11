function L = compute_combinatorial_laplacian(A)

% compute_laplacian - return the combinatorial laplacian of a given adjacency matrix
%
%   L = compute_combinatorial_laplacian(A);
%
%   L(i,j) = -1 if i \neq j and i is connected to j
%   L(i,i) = - sum_j L(i,j)
%   L(i,j) = 0  otherwise
%
%   This laplacian is symmetric, L=G'*G where G is the graph gradient,
%   computed with adjacency2incidence.
%
%   Copyright (c) 2004 Gabriel Peyre

error('Not used anymore');

Ic = adjacency2incidence(A);
L = Ic*Ic';

return;

w0 = sum( abs(Ic') );
w0 = diag( w0.^(-1/2) );
L = w0 * A * w0;
L = diag( sum(L) ) - L;
L = (L+L')/2; 


return;

%% old code

n = length(A);
L = zeros(n,n);

for i=1:n
    ni = sum( A(i,:) );
    for j=1:n
        nj = sum( A(:,j) );
        % number of neighbor
        if A(i,j) ==1
            L(i,j) = -1/sqrt(ni*nj);
        end
    end
end

for i=1:n
    L(i,i) = -sum( L(i,:) );
end