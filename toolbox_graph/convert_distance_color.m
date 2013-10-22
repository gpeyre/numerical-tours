function v = convert_distance_color(D,M, equalize)

% convert_distance_color - convert a distance function to a color image
%
%   A = convert_distance_color(D,M);
%
%   M is optional: background image.
%
%   Very useful to save a result of distance computation to an image file
%   with nice colors.
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    equalize = 0;
end

n = size(D,1);
if nargin<2
    M = ones(n);
end

c = jet(256);

U = D; 

if equalize
    U(U~=Inf) = perform_hist_eq(U(U~=Inf), 'linear');
end

U(U==Inf) = min(U(U~=Inf));
I = floor(255*rescale(U))+1;
v = c(I(:),:); 
v = reshape(v, [n n 3]);
A = repmat(rescale(M), [1 1 3]);
B = repmat(D, [1 1 3]);
v(B==Inf) = A(B==Inf);