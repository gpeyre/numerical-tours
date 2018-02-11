function [Y,I] = compute_max(X,d)

// compute_max - compute maximum along dimension d
//
//   [Y,I] = compute_max(X,d);
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)<2
    d = 1;
end

if d==1
    [Y,I] = max(X,'r');
elseif d==2
    [Y,I] = max(X,'c');
else
    error('Not implemented');
end

endfunction