function M = sum3(M,k)

// sum3 - sum along dimension 3
//
//  Apparently there is a strange bug with the windows version of scilab.
//
//  Gabriel Peyre (c) 2008

n = size(M,1); p = size(M,2); q = size(M,3);
M = reshape(M, [n*p, q]);
M = reshape( sum(M',1)', [n, p] );

endfunction