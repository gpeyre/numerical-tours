function [fx,fy,fz] = grad(M, options)

// grad - gradient, forward differences
//
//   [gx,gy] = grad(M, options);
// or
//   g = grad(M, options);
//
//   options.bound = 'per' or 'sym'
//   options.order = 1 (backward differences)
//                 = 2 (centered differences)
//
//   Works also for 3D array.
//   Assme that the function is evenly sampled with sampling step 1.
//
//   See also: div.
//
//   Copyright (c) Gabriel Peyre


options.null = 0;
bound = getoptions(options, 'bound', 'sym');
order = getoptions(options, 'order', 1);


// retrieve number of dimensions
nbdims = nb_dims(M);

n = size(M,1);


if strcmp(bound, 'sym')    
    if order==1
        fx = M([2:n n],:,:)-M;
    else
        fx = ( M([2:n n],:,:)-M([1 1:n-1],:,:) )/2;
        // boundary
        fx(1,:,:) = M(2,:,:)-M(1,:,:);
        fx(n,:,:) = M(n,:,:)-M(n-1,:,:);
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:n n],:)-M;
        else
            fy = ( M(:,[2:n n],:)-M(:,[1 1:n-1],:) )/2;
            // boundary
            fy(:,1,:) = M(:,2,:)-M(:,1,:);
            fy(:,n,:) = M(:,n,:)-M(:,n-1,:);
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:n n])-M;
        else
            fz = ( M(:,:,[2:n n])-M(:,:,[1 1:n-1]) )/2;
            // boundary
            fz(:,:,1) = M(:,:,2)-M(:,:,1);
            fz(:,:,n) = M(:,:,n)-M(:,:,n-1);
        end
    end
else
    if order==1
        fx = M([2:n 1],:,:)-M;
    else
        fx = ( M([2:n 1],:,:)-M([n 1:n-1],:,:) )/2;
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:n 1],:)-M;
        else
            fy = ( M(:,[2:n 1],:)-M(:,[n 1:n-1],:) )/2;
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:n 1])-M;
        else
            fz = ( M(:,:,[2:n 1])-M(:,:,[n 1:n-1]) )/2;
        end
    end
end

if argn(1)==1
    if nbdims==2
        fx = cat_3d(fx,fy);
    elseif nbdims==3
        error();
        fx = cat(4,fx,fy,fz);
    end
end

endfunction

// --- //
function c = cat_3d(a,b)
// I do not know why, but cat is not working properly on scilab ...
c = zeros(size(a,1), size(a,2), 2);
c(:,:,1) = a;
c(:,:,2) = b;
endfunction