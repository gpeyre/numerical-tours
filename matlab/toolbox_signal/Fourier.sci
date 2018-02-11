function F = Fourier(M)

// Fourier - centered Fourier transform
//
//  M = Fourier(F);
//
//  Copyright (c) 2009 Gabriel Peyre

n = size(M,1);
F = fftshift(fft2(M)/n);

endfunction