function v = cell_get(x,i)

// cell_get - get value
//
//   v = cell_get(x,i,v);
// is equivalent to 
//   v = x{i};
//
//   Copyright (c) 2008 Gabriel Peyre

v = x(i).entries;

endfunction
