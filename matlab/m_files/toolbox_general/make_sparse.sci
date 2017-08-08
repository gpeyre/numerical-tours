function A = make_sparse(i,j,x,n,m)

// make_sparse - wrapper for sparse matrix creation.

if argn(2)<3
    A = spzeros(n,n);
    return;
end

if argn(2)<4
    n = max(i);
end
if argn(2)<5
    m = max(j);
end

A = sparse( [i(:) j(:)], x(:), [n;m] );

endfunction