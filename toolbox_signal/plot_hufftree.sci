function plot_hufftree(T)

// plot_hufftree - plot a huffman tree 
//
//   plot_hufftree(T);
//
//   Copyright (c) 2008 Gabriel Peyre

plot_tree(T(1).entries,[0,0],1);
axis('tight');
axis('off');

endfunction

////
function plot_tree(T,x,j)

tw = 20;
lw = 1.5;
ms = 20;

if not(iscell(T))
    de = [-.02 -.2];
    de = [0 -0.1];
    // display a string
	xstring(x(1)+de(1),x(2)+de(2),num2str(T));
    // u = text(x(1)+de(1),x(2)+de(2),num2str(T));
    // set(u,'FontSize',tw);
else
    // display a split
    h0 = 1/2^j;
    h = 1 * h0;
    u = plot( [x(1) x(1)-h], [x(2) x(2)-1], '.-' );
    // set(u,'MarkerSize',ms);
    // set(u,'LineWidth',lw);
    u = plot( [x(1) x(1)+h], [x(2) x(2)-1], '.-' );
    // set(u,'MarkerSize',ms);
    // set(u,'LineWidth',lw);
    plot_tree(T(1).entries,[x(1)-h x(2)-1],j+1);
    plot_tree(T(2).entries,[x(1)+h x(2)-1],j+1);
end

endfunction
