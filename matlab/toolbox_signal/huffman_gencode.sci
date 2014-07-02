function C = huffman_gencode(T)

// huffman_gencode - generate a code associated to a Hufmann tree T
//
//   C = huffman_gencode(T);
//
//   Copyright (c) 2008 Gabriel Peyre

if not(iscell(T)) 
    C = cell(0);
    C = cell_set(C, T, -1);
elseif size(T,1)==1
    C = huffman_gencode( cell_get(T,1) );
    // remove traling -1
    for i=1:size(C,1)
        // C{i} = C{i}(1:end-1);
        v = cell_get(C,i);
        v = v(1:length(v)-1);
        C = cell_set(C,i, v );
    end
elseif size(T,1)==2
    C1 = huffman_gencode( cell_get(T,1) );
    C2 = huffman_gencode( cell_get(T,2) );
    C = cell(0);
    for i=1:size(C1,1)
        // C1{i}
        C1i = cell_get(C1,i);
        if not(isempty(C1i))
            C = cell_set(C, i, [0 C1i]);
        end    
    end
    for i=1:size(C2,1)
        // C2{i}
        C2i = cell_get(C2,i);
        if not(isempty(C2i))
            C = cell_set(C, i, [1 C2i]);
        end
    end
else
    error('Problem');
end


endfunction