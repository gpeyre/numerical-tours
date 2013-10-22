function T = compute_hufftree(p)

% compute_hufftree - build a huffman tree
%
%   T = compute_hufftree(p);
%
%   p(i) is the probability of having token i.
%
%   Copyright (c) 2008 Gabriel Peyre

m = length(p);

% build initial trees
T = {};
for i=1:m
    T{i} = i;
end

% build Huffman tree
while length(T)>1
    [v,I] = sort(p);
    q = sum(v(1:2));
    t = {T{I(1:2)}};
    % trimed tree
    T = {T{I(3:end)}};
    p = v(3:end);
    % add new nodes
    p(end+1) = q;
    T{end+1} = t;
end