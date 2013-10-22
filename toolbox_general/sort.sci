function [y,I] = sort(x)


// sort - emulation of matlab function
//
//  sort in increasing order
//
// Copyright (c) 2011 Gabriel Peyre

[y,I] = gsort(x, 'r', 'i');

endfunction