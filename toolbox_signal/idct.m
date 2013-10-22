function a = idct(b,n)
%IDCT Inverse discrete cosine transform.
%
%   X = IDCT(Y) inverts the DCT transform, returning the
%   original vector if Y was obtained using Y = DCT(X).
%
%   X = IDCT(Y,N) pads or truncates the vector Y to length N 
%   before transforming.
%
%   If Y is a matrix, the IDCT operation is applied to
%   each column.
%
%   See also FFT, IFFT, DCT.

%   Author(s): C. Thompson, 2-12-93
%              S. Eddins, 10-26-94, revised
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/04/15 01:11:58 $

%   References: 
%   1) A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
%   2) Wallace, "The JPEG Still Picture Compression Standard",
%      Communications of the ACM, April 1991.

if nargin == 0,
	error('Not enough input arguments.');
end

if isempty(b),
   a = [];
   return
end

% If input is a vector, make it a column:
do_trans = (size(b,1) == 1);
if do_trans, b = b(:); end
   
if nargin==1,
  n = size(b,1);
end
m = size(b,2);

% Pad or truncate b if necessary
if size(b,1)<n,
  bb = zeros(n,m);
  bb(1:size(b,1),:) = b;
else
  bb = b(1:n,:);
end

% Compute wieghts
ww = sqrt(2*n) * exp(j*(0:n-1)*pi/(2*n)).';

if rem(n,2)==1 | ~isreal(b), % odd case
  % Form intermediate even-symmetric matrix.
  ww(1) = ww(1) * sqrt(2);
  W = ww(:,ones(1,m));
  yy = zeros(2*n,m);
  yy(1:n,:) = W.*bb;
  yy(n+2:2*n,:) = -j*W(2:n,:).*flipud(bb(2:n,:));
  
  y = ifft(yy);

  % Extract inverse DCT
  a = y(1:n,:);

else % even case
  % Compute precorrection factor
  ww(1) = ww(1)/sqrt(2);
  W = ww(:,ones(1,m));
  yy = W.*bb;

  % Compute x tilde using equation (5.93) in Jain
  y = ifft(yy);
  
  % Re-order elements of each column according to equations (5.93) and
  % (5.94) in Jain
  a = zeros(n,m);
  a(1:2:n,:) = y(1:n/2,:);
  a(2:2:n,:) = y(n:-1:n/2+1,:);
end

if isreal(b), a = real(a); end
if do_trans, a = a.'; end
