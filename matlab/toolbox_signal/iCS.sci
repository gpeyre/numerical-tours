function M = iCS(y)

// iCS - dual of Fourier measurements
//
//  M = iCS(y);
//
//  Copyright (c) 2009 Gabriel Peyre

global Ics;
M = zeros(n,n);
M(Ics) = y;
M = iFourier(M);

endfunction