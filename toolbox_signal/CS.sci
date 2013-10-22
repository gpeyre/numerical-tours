function y = CS(M)

// iCS - Fourier measurements
//
//  y = iCS(M);
//
//  Copyright (c) 2009 Gabriel Peyre


global Ics;
y = Fourier(M);
y = y(Ics);

endfunction