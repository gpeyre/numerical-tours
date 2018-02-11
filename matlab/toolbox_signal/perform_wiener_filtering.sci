function [y,Hwien] = peform_wiener_filtering(x0,x,sigma)

// peform_wiener_filtering - perform Wiener filtering
//
//   [y,Hwien] = peform_wiener_filtering(x0,x,sigma);
//
//   x0 is the original (clean) image
//   x is the noisy image
//   sigma is the noise level.
//
//   Copyright (c) 2008 Gabriel Peyre

n = prod(size(x));

    // FFT-based wiener filtering (using the oracle fourier coefficients)
    x0f = fft(x0);
    Pxf = abs(x0f).^2; // power spectra
    Hf = Pxf./(Pxf + n*sigma^2); // filter fourier transform
    // compute convolution
    xf = fft(x);
    y = real( ifft(xf.*Hf) );
    if argn(2)>1
        Hwien = real( fftshift( ifft(Hf) ) );
    end

endfunction