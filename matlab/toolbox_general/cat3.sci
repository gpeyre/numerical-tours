function c = cat3(d, a,b)

// Debugged version of cat function
//
//  c = cat3(dimension, a,b);
//
// I do not know why, but cat is not working properly on scilab for large matrices ...
//
//  Copyright (c) Gabriel Peyre

if nb_dims(a)<d
    n = 1;
else
    n = size(a,d); 
end

if nb_dims(b)<d
    p = 1;
else
    p = size(b,d); 
end


s = size(a);
if nb_dims(a)<=2
    s(3) = 1;
end

if d==1
    c = zeros(n+p,s(2),s(3));
    c(1:n,:,:) = a;
    c(n+1:n+p,:,:) = b;
elseif d==2
    c = zeros(s(1), n+p, s(3));
    c(:,1:n,:) = a;
    c(:,n+1:n+p,:) = b;
elseif d==3
    c = zeros(s(1), s(2), n+p);
    c(:,:,1:n) = a;
    c(:,:,n+1:n+p) = b;
else
    error('Works only for 3D arrays');
end

c = squeeze(c);

endfunction