function b=dct(a,n)
%DCT  Discrete cosine transform.
%
%   Y = DCT(X) returns the discrete cosine transform of X.
%   The vector Y is the same size as X and contains the
%   discrete cosine transform coefficients.
%
%   Y = DCT(X,N) pads or truncates the vector X to length N 
%   before transforming.
%
%   If X is a matrix, the DCT operation is applied to each
%   column.  This transform can be inverted using IDCT.
%
%   See also FFT, IFFT, IDCT.

%   Author(s): C. Thompson, 2-12-93
%              S. Eddins, 10-26-94, revised
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/04/15 01:10:40 $

%   References: 
%   1) A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
%   2) Wallace, "The JPEG Still Picture Compression Standard",
%      Communications of the ACM, April 1991.


if nargin == 0,
	error('Not enough input arguments.');
end

if isempty(a)
   b = [];
   return
end

% If input is a vector, make it a column:
do_trans = (size(a,1) == 1);
if do_trans, a = a(:); end

if nargin==1,
  n = size(a,1);
end
m = size(a,2);

% Pad or truncate input if necessary
if size(a,1)<n,
  aa = zeros(n,m);
  aa(1:size(a,1),:) = a;
else
  aa = a(1:n,:);
end

% Compute weights to multiply DFT coefficients
ww = (exp(-i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
ww(1) = ww(1) / sqrt(2);

if rem(n,2)==1 | ~isreal(a), % odd case
  % Form intermediate even-symmetric matrix
  y = zeros(2*n,m);
  y(1:n,:) = aa;
  y(n+1:2*n,:) = flipud(aa);
  
  % Compute the FFT and keep the appropriate portion:
  yy = fft(y);  
  yy = yy(1:n,:);

else % even case
  % Re-order the elements of the columns of x
  y = [ aa(1:2:n,:); aa(n:-2:2,:) ];
  yy = fft(y);  
  ww = 2*ww;  % Double the weights for even-length case  
end

% Multiply FFT by weights:
b = ww(:,ones(1,m)) .* yy;

if isreal(a), b = real(b); end
if do_trans, b = b.'; end
