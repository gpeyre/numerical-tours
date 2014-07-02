function f = compute_gaussian_filter(n,s,N);

// compute_gaussian_filter - compute a 1D or 2D Gaussian filter.
//
//   f = compute_gaussian_filter(n,s,N);
//
//   'n' is the size of the filter, odd for no phase in the filter.
//       (if too small it will alterate the filter).
//       use n=[n1,n2] for a 2D filter
//   's' is the standard deviation of the filter.
//   'N' is the size of the big signal/image (supposed to lie in [0,1] or [0,1]x[0,1]).
//       use N=[N1,N2] for a 2D filter.
//
//   The equation (in 1D) is
//       f[k] = exp( -(x(k)^2/(2*s^2)) );
//   where x span [-1/2,1/2].
//
//   The filter is normalised so that it sums to 1.
//
//   Copyright (c) 2004 Gabriel Peyre

nd = 1;
if length(n)>1 & n(2)>1
    nd = 2;
end

if nd==2 & length(s)==1
    s = [s s];
end
if nd==2 & length(N)==1
    N = [N N];
end

if nd==1
    f = build_gaussian_filter_1d(n,s,N);
else
    f = build_gaussian_filter_2d(n,s,N);
end

endfunction


// build_gaussian_filter_2d - compute a 2D Gaussian filter.
//
//   f = build_gaussian_filter_2d(n,s,N);
//
//   'n' is the size of the filter, odd for no phase in the filter.
//       (if too small it will alterate the filter).
//   's' is the standard deviation of the filter.
//   'N' is the size of the big image (supposed to lie in [0,1]x[0,1]).
//
//   The filter is normalised so that it sums to 1.
//
//   Copyright (c) 2004 Gabriel Peyr?

function f = build_gaussian_filter_2d(n,s,N)

if argn(2)<2
    error('Not enough arguments.');
end
if argn(2)<3
    N = n;
end

if length(N)==1 | N(1)==1
    N = N(:); N = [N N];
end
if length(s)==1 | s(1)==1
    s = s(:); s = [s s];
end

if s<=0
    f = zeros(n);
    f(round((n-1)/2),round((n-1)/2)) = 1;
    return;
end

x = ( (0:n(1)-1)-(n(1)-1)/2 )/(N(1)-1);
y = ( (0:n(2)-1)-(n(2)-1)/2 )/(N(2)-1);
[Y,X] = meshgrid(y,x);
f = exp( -(X.^2/ (2*s(1)^2)) - (Y.^2/ (2*s(2)^2)) );
f = f / sum(f(:));

endfunction

// build_gaussian_filter_1d - compute a Gaussian filter.
//
//   f = build_gaussian_filter_1d(n,s,N);
//
//   Copyright (c) 2004 Gabriel Peyr?

function f = build_gaussian_filter_1d(n,s,N)

if argn(2)<2
    error('Not enough arguments.');
end
if argn(2)<3
    N = n;
end

if s<=0
    f = zeros(n,1);
    f(round((n-1)/2)) = 1;
    return;
end

x = ( (0:n-1)-(n-1)/2 )/(N-1);
f = exp( -x.^2/(2*s^2) );
f = f / sum(f(:));

endfunction