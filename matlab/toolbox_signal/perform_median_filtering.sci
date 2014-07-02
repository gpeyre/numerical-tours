function M = perform_median_filtering(M,k)

// perform_median_filtering - perform moving average median
//
//   M = perform_median_filtering(M,k);
//
//   k is the half width of the window (detult k=1).
//
//   This filtering method is very efficient to remove impulsive or salt and
//   peper noise.
//
//   Copyright (c) 2007 Gabriel Peyre

if argn(2)<2
    k = 1;
end
w = 2*k+1;
n = size(M,1);

options.sampling = 'uniform';
H = compute_patch_library(M,k,options);
H = reshape(H, [w*w n*n]);
H = median(H, 'r'); 
M = reshape(H, [n,n]);


endfunction

function H = compute_patch_library(M,w,options)

// [H,X,Y] = compute_patch_library(M,w,options);
//
//   M is the texture
//   w is the half-width of the patch, so that each patch
//       has size (2*w+1,2*w+1,s) where s is the number of colors.
//
//   H(:,:,:,i) is the ith patch (can be a color patch).
//   X(i) is the x location of H(:,:,:,i) in the image M.
//   Y(i) is the y location of H(:,:,:,i) in the image M.
//
//   options.sampling can be set to 'random' or 'uniform'
//   If options.sampling=='random' then you can define
//       options.nbr as the number of random patches and/or
//       option.locations_x and option.locations_y.
//
//   If w is a vector of integer size, then the algorithm 
//   compute a set of library, one for each size, 
//   and use the same location X,Y for each patch.
//   H{k} is the k-th library.
//
//   You can define options.wmax to avoid too big patch.
//   In this case, if w>wmax, the texture M will be filtered, and the patches
//   will be subsampled to match the size.
//
//   Copyright (c) 2006 Gabriel Peyre

options.null = 0;


n = size(M,1);
s = 1;
ww = 2*w+1;

// do some padding to avoid boundary problems
M = symmetric_extension(M,w);


[Y,X] = meshgrid(w+1:w+n, w+1:w+n);
X = X(:);
Y = Y(:);
p = length(X);

H = zeros(ww,ww,s,p);

    ////// in this case, a fast sampling can be used //////
    [dY,dX] = meshgrid(-w:w,-w:w);
    Xp = repmat( reshape(X,[1,1,1,p]) ,[ww ww s 1]) + repmat(dX,[1 1 s p]);
    Yp = repmat( reshape(Y,[1,1,1,p]) ,[ww ww s 1]) + repmat(dY,[1 1 s p]);
//    I = sub2ind([size(M) 1], Xp,Yp,Cp);
    I = Xp + (Yp-1) * size(M,1);
    H = reshape(M(I(:)),size(I));

endfunction

function M_padded = symmetric_extension(M,k)

// symmetric_extension - perform a symmetric extension of the signal.
//
//   M_padded = symmetric_extension(M,k);
//
// M can be 1D or 2D array
// If size(M)=[n,p], the result is of size [n+2*k,p+2*k]
//
//   Copyright (c) 2004 Gabriel Peyre

n1 = size(M,1);
n2 = size(M,2);
m1 = n1+2*k;
m2 = n2+2*k;

if nb_dims(M)==1
    M = M(:);
    M_padded = [ M(k:-1:1); M; M(n1:-1:n2-k+1) ];
elseif nb_dims(M)==2
    M_padded = zeros(n1+2*k,n2+2*k);
    M_padded(k+1:m1-k,k+1:m2-k) = M;
    // extension
    M_padded(1:k,:) = M_padded(2*k:-1:k+1,:);
    M_padded(m1-k+1:m1,:) = M_padded(m1-k:-1:m1-2*k+1,:);
    M_padded(:,1:k) = M_padded(:,2*k:-1:k+1);
    M_padded(:,m2-k+1:m2) = M_padded(:,m2-k:-1:m2-2*k+1);
else
    error('Only supported for array of dimension less than 2.')
end


endfunction
