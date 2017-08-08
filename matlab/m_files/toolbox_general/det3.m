function d = det3(A)

% det3 - 3x3 determinant
%
%   d = det3(A);
%
%   A is a [3 3 n] matrix, d is a vector of size d where
%   d(i)=det(A(:,:,i)).
%
%   It computes explicitely the det by expanding along the 1st columns 2x2
%   determinant. Faster than recursive loops.
%
%   Note: works also for 2x2 det.
%   Note: for higher order det, use loop
%
%   Copyright (c) 2008 Gabriel Peyre

if size(A,1)==2 && size(A,2)==2
    %% 2x2 Det %%
    d = A(1,1,:).*A(2,2,:) - A(1,2,:).*A(2,1,:);
elseif size(A,1)==3 && size(A,2)==3
    %% 3x3 Det %%
    d = A(1,1,:).*( A(2,2,:).*A(3,3,:) - A(2,3,:).*A(3,2,:) ) - ...
        A(1,2,:).*( A(2,1,:).*A(3,3,:) - A(2,3,:).*A(3,1,:) ) + ...
        A(1,3,:).*( A(2,1,:).*A(3,2,:) - A(2,2,:).*A(3,1,:) );
else
    n = size(A,3);
    d = zeros(n,1);
    for i=1:n
        d(i) = det(A(:,:,i));
    end
end
d = d(:);
