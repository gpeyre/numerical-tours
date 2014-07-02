function x = repmat(x, s)

// repmat - duplicate an array
//
// x = repmat(x, s)
//
//  Copyright (c) 2008 Gabriel Peyre



x = mtlb_repmat(x,s);

return;

// warning : very innefficient
for i=1:length(s)
    y = x;
    for j=1:s(i)-1
       y = cat(i,y,x);
    end
    x = y;
end


endfunction