function M = crop(M,n,c)

// crop - crop an image to reduce its size
//
// M = crop(M,n,c)
//
//  Copyright (c) 2008 Gabriel Peyre

n0 = size(M);
// number of dimensions
d = nb_dims(M);

if argn(2)<2
    n = round( n0/2 );
end
if argn(2)<3
    c = round( n0/2 );
end

if isempty(n)
    return;
end

if length(n)==1
    n = n * ones(d,1);
end
if length(c)==1
    c = c * ones(d,1);
end


c = round(c);
if length(n)<3
    n(3) = 1;
end

selx = c(1)-ceil(n(1)/2)+1:c(1)+floor(n(1)/2);
sely = c(2)-ceil(n(2)/2)+1:c(2)+floor(n(2)/2);
selz = 1;
if nb_dims(M)==3
    // 3D cropping
    selz = c(3)-ceil(n(3)/2)+1:c(3)+floor(n(3)/2);
    sely(selz<1) = []; 
    sely(selz>n0(3)) = [];
end

selx(selx<1) = []; 
selx(selx>n0(1)) = [];
sely(sely<1) = []; 
sely(sely>n0(2)) = [];

M = M(selx,sely,selz);

endfunction