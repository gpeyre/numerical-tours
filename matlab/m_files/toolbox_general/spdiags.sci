function A = spdiags(x, u, n,n)

// Sparse diagonal matrix

A = make_sparse(1:n,1:n,x);


endfunction