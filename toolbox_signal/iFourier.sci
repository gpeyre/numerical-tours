function M = iFourier(F)

// iFourier - inverse centered Fourier transform
//
//  M = iFourier(F);
//
//  Copyright (c) 2009 Gabriel Peyre


n = size(F,1);
M = real( ifft2( fftshift(F)*n ) );

endfunction