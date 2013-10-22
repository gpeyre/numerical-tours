function M1 = perform_color_matching(M1,M,niter,finishmatch)

% perform_color_matching - perform histogram matching of color images
%
%   M1 = perform_color_matching(M1,M,niter);
%
%   Perform iterative matching over randomized colorspaces.
%
%   Copyright (c) 2009 Gabriel Peyre

n = size(M1,1);
p = size(M1,3);

if nargin<3
    niter = 10;
end
if nargin<4
    finishmatch = 0;
end
% iterates matching
for k=1:niter
    [U,R] = qr(randn(p));
    d = reshape(M,[n^2 p])*U;
    d1 = reshape(M1,[n^2 p])*U;
    for c=1:p
        d1(:,c) = perform_hist_eq(d1(:,c),d(:,c));
    end
    M1 = reshape(d1*U',[n n p]);
end

if finishmatch
    % finish the match
    u = randn(1,1,p);
    d  = sum( M.*repmat(u,[n n]), 3);
    d1 = sum( M1.*repmat(u,[n n]), 3);
    [tmp,I] = sort(d(:)); [tmp,I1] = sort(d1(:));
    if p==2
        M1([I1;I1+n^2]) = M([I;I+n^2]);
    elseif p==3
        M1([I1;I1+n^2;I1+2*n^2]) = M([I;I+n^2;I+2*n^2]);
    elseif p==4
        M1([I1;I1+n^2;I1+2*n^2;I1+3*n^2]) = M([I;I+n^2;I+2*n^2;I+3*n^2]);
    else
        error('Not implemented');
    end
end
