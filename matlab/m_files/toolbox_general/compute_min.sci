function [Y,I] = compute_min(X,d)

// compute_min - compute min along dimension d
//
//   [Y,I] = compute_min(X,d);
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)<2
    d = 1;
end

if d==1
    [Y,I] = min(X,'r');
elseif d==2
    [Y,I] = min(X,'c');
else
    error('Not implemented');
end

endfunction