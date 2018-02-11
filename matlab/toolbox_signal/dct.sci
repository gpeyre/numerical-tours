function b=dct(a,n)

//DCT  Discrete cosine transform.
//
//   Y = DCT(X) returns the discrete cosine transform of X.
//   The vector Y is the same size as X and contains the
//   discrete cosine transform coefficients.
//
//   Y = DCT(X,N) pads or truncates the vector X to length N 
//   before transforming.
//
//   If X is a matrix, the DCT operation is applied to each
//   column.  This transform can be inverted using IDCT.
//
//   See also FFT, IFFT, IDCT.

//   Author(s): C. Thompson, 2-12-93
//              S. Eddins, 10-26-94, revised
//   Copyright 1988-2004 The MathWorks, Inc.
//   $Revision: 1.7.4.2 $  $Date: 2004/12/26 22:15:37 $

//   References: 
//   1) A. K. Jain, "Fundamentals of Digital Image
//      Processing", pp. 150-153.
//   2) Wallace, "The JPEG Still Picture Compression Standard",
//      Communications of the ACM, April 1991.


if argn(2) == 0,
	error('Not enough input arguments.');
end

if isempty(a)
   b = [];
   return
end


if argn(2)==1,
  n = size(a,1);
end

// If input is a vector, make it a column:
do_trans = (size(a,1) == 1);
if do_trans, a = a(:); end

m = size(a,2);

// Pad or truncate input if necessary
if size(a,1)<n,
  aa = zeros(n,m);
  aa(1:size(a,1),:) = a;
else
  aa = a(1:n,:);
end

// Compute weights to multiply DFT coefficients
ww = (exp(-sqrt(-1)*(0:n-1)*pi()/(2*n))/sqrt(2*n)).';
ww(1) = ww(1) / sqrt(2);

if pmodulo(n,2)==1 | ~isreal(a), // odd case
  // Form intermediate even-symmetric matrix
  y = zeros(2*n,m);
  y(1:n,:) = aa;
  y(n+1:2*n,:) = flipud(aa);
  
  // Compute the FFT and keep the appropriate portion:
  yy = myfft(y);  
  yy = yy(1:n,:);

else // even case
  // Re-order the elements of the columns of x
  y = [ aa(1:2:n,:); aa(n:-2:2,:) ];
  yy = myfft(y);  
  ww = 2*ww;  // Double the weights for even-length case  
end

// Multiply FFT by weights:
b = ww(:,ones(1,m)) .* yy;

if isreal(a), b = real(b); end
if do_trans, b = b.'; end

endfunction

function y = flipud(x)
y = x($:-1:1,:);
endfunction

function y = myfft(x)
y = x;
for i=1:size(x,2)
    y(:,i) = fft(x(:,i));
end
endfunction