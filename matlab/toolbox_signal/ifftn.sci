function y = ifftn(x)

// ifftn - Inverse Fourier transform (for Matlab compatibility)
//
//   y = ifftn(x);
//
//  Copyright (c) 2008 Gabriel Peyre

y = ifft(x);

endfunction