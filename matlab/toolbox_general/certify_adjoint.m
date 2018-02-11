function e = certify_adjoint(A,AS,d)

% certify_adjoint - test of adjointess
%
%   e = certify_adjoint(A,AS,d);
%
%   d is the input dimension of A.
%   e = |<A(x),y> - <x,AS(y)>|/|<A(x),y>|
%
%   for x=randn(d) and y=randn(size(A(x)).
%
%   Copyright (c) 2012 Gabriel Peyre

dotp = @(a,b)sum(a(:).*b(:));

x = randn(d);
Ax = A(x);
y = randn(size(Ax));
ASy = AS(y);

e = abs( dotp(Ax,y)-dotp(x,ASy) ) ./ abs(dotp(Ax,y));