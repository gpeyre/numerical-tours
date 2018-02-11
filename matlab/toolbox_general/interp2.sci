function M1 = interp2(Y,X,M,Yi,Xi)

// interp2 - bilinear interpolation
//
//  M1 = interp2(X,Y,M,Xi,Yi);
//
//  Copyright (c) 2008 Gabriel Peyre

X = X(:,1); Y = Y(1,:)';
M1 = linear_interpn( Yi,Xi, Y,X, M)';

endfunction
