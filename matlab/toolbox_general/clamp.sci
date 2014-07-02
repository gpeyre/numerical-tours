function y = clamp(x,a,b)

// clamp - clamp a value or an array
//
//  y = clamp(x,a,b);
//
//  Copyright (c) 2008 Gabriel Peyre

if argn(2)<2
    a = 0;
end

if argn(2)<3
    b = 1;
end

if iscell(x)
    for i=1:length(x)
        y(i) = clamp(x(i),a,b);
    end
    return;
end

y = max(x(:),a);
y = min(y,b);
y = reshape(y,size(x));

endfunction