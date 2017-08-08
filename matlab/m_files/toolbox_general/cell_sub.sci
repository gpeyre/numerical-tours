function y = cell_sub(x,sel)

// cell_sub - extract a sub-cell array
//
//   y = cell_sub(x,sel);
// is equivalent with
//   y = {x{sel}};
//
// Copyright (c) Gabriel Peyre 2008

y = cell(0);

for i=1:length(sel)
    y(i).entries = x(sel(i)).entries;
end

endfunction