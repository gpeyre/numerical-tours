function y = randn(n,p,q)

// randn - Gaussian number generator
//
// y = randn(n,p);
// y = randn([n p]);
// y = randn(n,p,q);
// y = randn([n p q]);
//
//  Copyright (c) 2008 Gabriel Peyre


if isstr(n)
    y = 0;
    return;
end

if argn(2)==2
    n = [n p];
end
if argn(2)==3
    n = [n p q];
end

y = squeeze(mtlb_randn(n));

if n(1)==1 & n(2)~=1 & length(n)==2
    y = y(:)';
end

endfunction