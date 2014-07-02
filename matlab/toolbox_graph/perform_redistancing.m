function D1 = perform_redistancing(D, options)

% perform_redistancing - redistance a function
%
%   D1 = perform_redistancing(D, options);
%
%   Compute a signed distance function D1 that have the same 0 level set as D.
%   You can turn off interpolation by using options.use_interpolation=0.
%
%   Note that the distance function is computed with 1 pixel = distance of
%   1, so the overall image range over a 1...n=size(D,1)
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
use_interpolation = getoptions(options, 'use_interpolation', 1);

n = size(D,1);

% horizontal
P1 = D(1:end-1,:); P2 = D(2:end,:);
P = (P1.*P2)<=0;
d = abs(P1-P2); d(d<eps) = 1;
v1 = abs(P1)./d;
v2 = abs(P2)./d;
Ah = ([P; zeros(1,n)] + [zeros(1,n); P])>0;
Vh = max([v1; zeros(1,n)], [zeros(1,n); v2]);
% vertical
P1 = D(:,1:end-1); P2 = D(:,2:end);
P = (P1.*P2)<=0;
d = abs(P1-P2); d(d<eps) = 1;
v1 = abs(P1)./d;
v2 = abs(P2)./d;
Av = ([P, zeros(n,1)] + [zeros(n,1), P])>0;
Vv = max([v1, zeros(n,1)],[zeros(n,1), v2]);

V = zeros(n);
I = find(Ah>0);
V(I) = Vh(I);
I = find(Av>0);
V(I) = max(V(I),Vv(I));

I = find(V~=0);
[x,y] = ind2sub(size(D),I); 
start_points = [x(:)'; y(:)'];
vals = V(I);

% with interpolation
options.nb_iter_max = Inf;
if use_interpolation
    options.values = vals/n;
else
    options.values = [];    
end
D1 = perform_fast_marching(ones(n), start_points, options);
D1 = D1*n;
D1(D<0) = -D1(D<0);