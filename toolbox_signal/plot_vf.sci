function plot_vf(v,M)

// plot_vf - plot 2D vector field
//
//  plot_vf(v,M);
//
//  Copyright (c) 2008 Gabriel Peyre


sc = 1;
n = size(v,1);
champ(1:n,1:n, v(:,:,1), v(:,:,2), rect=[1,1,n,n], arfact=sc, strf="021");

endfunction