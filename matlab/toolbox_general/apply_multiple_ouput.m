function y = apply_multiple_ouput(f, g, x)

% apply_multiple_ouput - compose multiple arguments functions
%
%   y = apply_multiple_ouput(f, g, x);
%
% compute y=f(g(x)) where f takes several input and g computes several outputs.
%
%   Copyright (c) 2010 Gabriel Peyre

if not(nargout(g)==nargin(f))
    error('one must have nargout(g)==nargin(f)');
end

switch nargout(g)    
    case 1
        [a{1}] = g(x);
        y = f(a{1});
    case 2
        [a{1}, a{2}] = g(x);
        y = f(a{1}, a{2});
    case 3
        [a{1}, a{2}, a{3}] = g(x);
        y = f(a{1}, a{2}, a{3});
    case 4
        [a{1}, a{2}, a{3}, a{4}] = g(x);
        y = f(a{1}, a{2}, a{3}, a{4});
    otherwise
        error('Not implemented');
end
        
