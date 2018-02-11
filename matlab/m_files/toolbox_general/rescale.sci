function y = rescale(x,a,b)

// rescale - rescale data in [a,b]
//
// y = rescale(x,a,b)
//
//  Copyright (c) 2008 Gabriel Peyre



if argn(2)<2
    a = 0;
end

if argn(2)<3
    b = 1;
end

m = min(x(:));
M = max(x(:));

if M-m<1e-9
    y = x;
else
    y = (b-a) * (x-m)/(M-m) + a;
end




endfunction