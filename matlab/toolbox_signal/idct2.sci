function a = idct2(arg1,mrows,ncols)
//IDCT2 Compute 2-D inverse discrete cosine transform.
//   B = IDCT2(A) returns the two-dimensional inverse discrete
//   cosine transform of A.
//
//   B = IDCT2(A,[M N]) or B = IDCT2(A,M,N) pads A with zeros (or
//   truncates A) to create a matrix of size M-by-N before
//   transforming. 
//
//   For any A, IDCT2(DCT2(A)) equals A to within roundoff error.
//
//   The discrete cosine transform is often used for image
//   compression applications.
//
//   Class Support
//   -------------
//   The input matrix A can be of class double or of any
//   numeric class. The output matrix B is of class double.
//
//   See also DCT2, DCTMTX, FFT2, IFFT2.

//   Copyright 1993-2003 The MathWorks, Inc.  
//   $Revision: 5.17.4.1 $  $Date: 2003/01/26 05:55:39 $

//   References: 
//   1) A. K. Jain, "Fundamentals of Digital Image
//      Processing", pp. 150-153.
//   2) Wallace, "The JPEG Still Picture Compression Standard",
//      Communications of the ACM, April 1991.

// checknargin(1,3,nargin,mfilename);

[m, n] = size(arg1);
// Basic algorithm.
if (argn(2) == 1),
  if (m > 1) & (n > 1),
    a = idct(idct(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

// Padding for vector input.

b = arg1;
if argn(2)==2, 
    ncols = mrows(2); 
    mrows = mrows(1); 
end

mpad = mrows; npad = ncols;
if m == 1 & mpad > m, b(2, 1) = 0; m = 2; end
if n == 1 & npad > n, b(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   // For row vector.

// Transform.

a = idct(b, mpad);
if m > 1 & n > 1, a = idct(a.', npad).'; end


endfunction
