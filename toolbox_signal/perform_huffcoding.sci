function y = perform_huffcoding(x,T,direction)

// perform_huffcoding - perform huffman coding
//
//   y = perform_huffcoding(x,T,direction);
//
//   direction=+1 for coding
//   direction=-1 for decoding
//
//   T is a Huffman tree, computed with compute_hufftree
//
//   Copyright (c) 2008 Gabriel Peyre

if direction==1
    ////// CODING //////
    C = huffman_gencode(T);
    m = size(C,1);
    x = round(x); 
    if min(x)<1 | max(x)>m
        error('Too small or too large token');
    end
    y = [];
    for i=1:length(x)
        y = [y cell_get(C, x(i))];
    end
else
    ////// DE-CODING //////  
 
// initial pointer on the tree: on the root
t = cell_get(T,1);
// initial empty decoded stream
x1 = [];
// initial stream buffer
y1 = x;
while not(isempty(y1))
    // go down in the tree
    if y1(1)==0
        t = cell_get(t,1);
    else
        t = cell_get(t,2);
    end
    // remove the symbol from the stream buffer
    y1(1) = [];
    if not(iscell(t))
        // we are on a leaf of the tree: output symbol
        x1 = [x1 t];
        t = cell_get(T,1);
    end
end
y = x1(:);


end



endfunction

