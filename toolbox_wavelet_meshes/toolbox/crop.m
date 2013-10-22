function M = crop(M,n,c)

% crop - crop an image to reduce its size
%
%   M = crop(M,n,c);
%
%   n is the new size of the image
%   c is the center of the grop
%
%   Copyright (c) 2007 Gabriel Peyre

n0 = size(M);
n0(3) = size(M,3);

d = 2;
if size(M,3)>3
    % multi-dimensional array
    d = 3;
end

if nargin<2
    n = round( n0/2 );
end
if nargin<3 || isempty(c)
    c = round( n0/2 );
end

if isempty(n)
    return;
end

if length(n)==1
    n = repmat(n,d,3);
end
if length(c)==1
    c = repmat(c,d,3);
end
c = round(c);
if length(n)<3
    n(3) = 1;
end

selx = c(1)-ceil(n(1)/2)+1:c(1)+floor(n(1)/2);
sely = c(2)-ceil(n(2)/2)+1:c(2)+floor(n(2)/2);
if d==2
    selz = 1:size(M,3);
else
    selz = c(3)-ceil(n(3)/2)+1:c(3)+floor(n(3)/2);
end

selx(selx<1) = []; 
selx(selx>n0(1)) = [];
sely(sely<1) = []; 
sely(sely>n0(2)) = [];
selz(selz<1) = []; 
selz(selz>n0(3)) = [];

M = M(selx,sely,selz);