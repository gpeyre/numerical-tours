function y = perform_huffcoding(x,T,dir)

% perform_huffcoding - perform huffman coding
%
%   y = perform_huffcoding(x,T,dir);
%
%   dir=+1 for coding
%   dir=-1 for decoding
%
%   T is a Huffman tree, computed with compute_hufftree
%
%   Copyright (c) 2008 Gabriel Peyre

if dir==1
    %%% CODING %%%
    C = huffman_gencode(T);
    m = length(C);
    x = round(x); 
    if min(x)<1 || max(x)>m
        error('Too small or too large token');
    end
    y = [];
    for i=1:length(x)
        y = [y C{x(i)}];        
    end
else
    %%% DE-CODING %%%    
    t = T{1};
    y = [];
    while not(isempty(x))
        if x(1)==0
            t = t{1};
        else
            t = t{2};
        end
        x(1) = [];
        if not(iscell(t))
            y = [y t];
            t = T{1};
        end        
    end
    y = y(:);
end



%%

function C = huffman_gencode(T)

if not(iscell(T)) % || not(length(T)==2)
    C = {};
    C{T} = -1;
elseif length(T)==1
    C = huffman_gencode(T{1});
    % remove traling -1
    for i=1:length(C)
        C{i} = C{i}(1:end-1);
    end
elseif length(T)==2
    C1 = huffman_gencode(T{1});
    C2 = huffman_gencode(T{2});
    C = {};
    for i=1:length(C1)
        if not(isempty(C1{i}))
            C{i} = [0 C1{i}];
        end    
    end
    for i=1:length(C2)
        if not(isempty(C2{i}))
            C{i} = [1 C2{i}];
        end
    end
else
    error('Problem');
end
