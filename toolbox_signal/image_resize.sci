function M = image_resize(M,n,p)

// image_resize - resize an image using bicubic interpolation
//
//   M = image_resize(M,n,p)
//
//  n,p is the new size (works for n=p for now).
//
//   Copyright (c) 2008 Gabriel Peyre

if nb_dims(M)==3
    // color images
    M0 = M;
    d = size(M,3);
    M = zeros(n,n,d);
    for i=1:d
        M(:,:,i) = image_resize(M0(:,:,i),n,n);
    end
    return;
end

n0 = size(M,1);
p0 = size(M,2);
X = linspace(0,1,n0);
Y = linspace(0,1,p0);
Xi = linspace(0,1,n);
Yi = Xi;

[Xi,Yi] = meshgrid( Xi,Yi );
M = linear_interpn( Xi,Yi, X,Y, M)';


endfunction