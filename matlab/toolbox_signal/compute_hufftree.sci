function T = compute_hufftree(p)

// compute_hufftree - build a huffman tree
//
//   T = compute_hufftree(p);
//
//   p(i) is the probability of having token i.
//
//   Copyright (c) 2008 Gabriel Peyre

m = length(p);

// build initial trees
T = cell(0);
for i=1:m
    T = cell_set(T,i,i);
end

// iterative merging of the leading probabilities
while length(p)>1   
    // sort in decaying order the probabilities
    [v,I] = sort(p);
    if v(1)>v(length(v))
        v = reverse(v); I = reverse(I);
    end 
    q = sum(v(1:2));
    t = cell_sub(T, I(1:2));
    // trimed tree
    T = cell_sub(T, I(3:length(I)) );
    p = v(3:length(v));
    // add a new node with the corresponding probability
    p(length(p)+1) = q;
    T = cell_set(T, length(p), t);
end

endfunction