function s = randperm(n)

// randperm - random permuation of n points
// 
// s = randperm(n);
//
// Copyright (c) 2008 Gabriel Peyre

[tmp,s] = sort(rand(n,1));

endfunction