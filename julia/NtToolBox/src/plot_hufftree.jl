function plot_hufftree(T)
    """
        plot_hufftree - plot a huffman tree

        plot_hufftree(T);

        Copyright (c) 2008 Gabriel Peyre
    """
    figure(figsize=(10,10))
    plot_tree(T,[0,0],1)
    xlim(-1,1)
    axis("off")
end

function plot_tree(T,x,j)
    tw = 15;
    lw = 1.5;
    ms = 10;

    if typeof(T)==String
        de = [-.02 -.2]
        # display a string
        u = text(x[1]+de[1],x[2]+de[2],T);
        #set(u,"FontSize",tw);
    else
        # display a split
        h0 = 1/2^j;
        h = 1 * h0;
        plot( [x[1]; x[1]-h], [x[2]; x[2]-1], ".-" ,ms=ms,lw=lw,c="b")
        #set(u,"MarkerSize",ms);
        #set(u,"LineWidth",lw);
        plot( [x[1]; x[1]+h], [x[2]; x[2]-1], ".-",ms=ms,lw=lw,c="b")
        #set(u,"MarkerSize",ms);
        #set(u,"LineWidth",lw);
        #xlim(-1,1)
        plot_tree(T[1],[x[1]-h x[2]-1],j+1)
        plot_tree(T[2],[x[1]+h x[2]-1],j+1)
    end
end
