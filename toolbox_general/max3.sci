function x = max3(x,v)

// max3 - debuggued version of max
//
//  x = max3(x,v);
//
//  Copyright (c) 2009 Gabriel Peyre

if nb_dims(v)==1
    x = reshape( max(x(:),v), size(x) );
else
    x = max(x,y);
end

endfunction