function [y,Hwien] = peform_wiener_filtering(x0,x,sigma)

% peform_wiener_filtering - perform Wiener filtering
%
%   [y,Hwien] = peform_wiener_filtering(x0,x,sigma);
%
%   x0 is the original (clean) image
%   x is the noisy image
%   sigma is the noise level.
%
%   Copyright (c) 2008 Gabriel Peyre

n = prod(size(x));



if nb_dims(x)==1
    % FFT-based wiener filtering (using the oracle fourier coefficients)
    x0f = fft(x0);
    Pxf = abs(x0f).^2; % power spectra
    Hf = Pxf./(Pxf + n*sigma^2); % filter fourier transform
    % compute convolution
    xf = fft(x);
    y = real( ifft(xf.*Hf) );
    if nargout>1
        Hwien = real( fftshift( ifft(Hf) ) );
    end
else
    % FFT-based wiener filtering (using the oracle fourier coefficients)
    x0f = fft2(x0);
    Pxf = abs(x0f).^2; % power spectra
    Hf = Pxf./(Pxf + n*sigma^2); % filter fourier transform
    % compute convolution
    xf = fft2(x);
    y = real( ifft2(xf.*Hf) );
    if nargout>1
        Hwien = real( fftshift( ifft2(Hf) ) );
    end
end