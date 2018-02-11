function b=dct2(arg1,mrows,ncols)
//DCT2 Compute 2-D discrete cosine transform.
//   B = DCT2(A) returns the discrete cosine transform of A.
//   The matrix B is the same size as A and contains the
//   discrete cosine transform coefficients.
//
//   B = DCT2(A,[M N]) or B = DCT2(A,M,N) pads the matrix A with
//   zeros to size M-by-N before transforming. If M or N is
//   smaller than the corresponding dimension of A, DCT2 truncates
//   A. 
//
//   This transform can be inverted using IDCT2.
//
//   Class Support
//   -------------
//   A can be numeric or logical. The returned matrix B is of 
//   class double.
//
//   Example
//   -------
//       RGB = imread('autumn.tif');
//       I = rgb2gray(RGB);
//       J = dct2(I);
//       imshow(log(abs(J)),[]), colormap(jet), colorbar
//
//   The commands below set values less than magnitude 10 in the
//   DCT matrix to zero, then reconstruct the image using the
//   inverse DCT function IDCT2.
//
//       J(abs(J)<10) = 0;
//       K = idct2(J);
//       imview(I)
//       imview(K,[0 255])
//
//   See also FFT2, IDCT2, IFFT2.

//   Copyright 1993-2003 The MathWorks, Inc.  
//   $Revision: 5.22.4.2 $  $Date: 2003/05/03 17:50:23 $

//   References: 
//        1) A. K. Jain, "Fundamentals of Digital Image
//           Processing", pp. 150-153.
//        2) Wallace, "The JPEG Still Picture Compression Standard",
//           Communications of the ACM, April 1991.

[m, n] = size(arg1);
// Basic algorithm.
if (argn(2) == 1),
  if (m > 1) & (n > 1),
    b = dct(dct(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

// Padding for vector input.
a = arg1;
if argn(2)==2, ncols = mrows(2); mrows = mrows(1); end
mpad = mrows; npad = ncols;
if m == 1 & mpad > m, a(2, 1) = 0; m = 2; end
if n == 1 & npad > n, a(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   // For row vector.

// Transform.

b = dct(a, mpad);
if m > 1 & n > 1, b = dct(b.', npad).'; end


endfunction
