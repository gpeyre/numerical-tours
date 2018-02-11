function A = squeeze(A)

// squeeze - remove singleton dimensions
//
//  A = squeeze(A);
//
//  Copyright (c) 2008 Gabriel Peyre

d = size(A);

d1 = [];
for i=1:length(d)
    if d(i)~=1
        d1($+1) = d(i);
    end
end
A = matrix(A,d1);

endfunction