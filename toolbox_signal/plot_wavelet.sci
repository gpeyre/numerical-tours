function Mout = plot_wavelet(M, Jmin, options)

// plot_wavelet - plot 2D wavelet transform stored using Mallat's ordering.
//
//   plot_wavelet(MW, Jmin, options);
//
//   'MW' is the wavelet transform (in Mallat's ordering, not lifting inplace
//       ordering).
//   'Jmin' is the minimum scale of the transform.
//
//   You can set options.style, options.edge_color, options.renormalize
//       options.line_width.
//
//   Works with 1D and 2D wavelets
//   
//   Copyright (c) 2006 Gabriel Peyr?e

if argn(2)<2
    Jmin = 1;
end

options.null = 0;
style = getoptions(options, 'style', 'real');
edge_color = getoptions(options, 'edge_color', 'r');
renormalize = getoptions(options, 'renormalize', 0);
lw = getoptions(options, 'line_width', 2);
gam = getoptions(options, 'gamma', .8);

n = size(M,1);
Jmax = log2(n)-1;

if size(M,1)==1 | size(M,2)==1
    // 1D Signal
//    hold on;
    plot(M); axis tight;
    a = min(M); b = max(M);
    for j=Jmin:Jmax
        plot(2^j*[1 1], [a b], 'r--');
    end
//    hold off;
	Mout = [];
    return;
end


separable = getoptions(options, 'separable',0);
if separable
    a = [1 Jmin:Jmax];
    for j1=a
        sel1 = 2^j1+1:2^(j1+1);
        if j1==1
            sel1 = 1:2^Jmin;
        end
        for j2=a
            sel2 = 2^j2+1:2^(j2+1);
            if j2==1
                sel2 = 1:2^Jmin;
            end
            M1 = M(sel1,sel2,:);
            M(sel1,sel2,:) = M(sel1,sel2,:) / max(max( abs(M(sel1,sel2,:)) ));
            st = stdev(M1(:));
            if st<1e-9
                st=1;
            end
            M1 = 0.3 * M1/st;
           	M1 = clamp(M1, -1,1);
            M(sel1,sel2,:) = (M1+1)/2; // rescale( M1 );
        end
    end
    // hold on;
    s = [1/2 n-1/2];
//    imagesc(s,s,M);
//    axis image; axis off;
//    colormap gray(256);
    xset('thickness',3);
    imageplot(M);
    // display axis separation
    for j=Jmax:-1:Jmin;        
        x = [0 n];
        y = [2^j 2^j];
        h = myplot(x,y, edge_color);
        // set(h, 'LineWidth', lw);
        y = [0 n];
        x = [2^j 2^j];
        h = myplot(x,y, edge_color);
        // set(h, 'LineWidth', lw);
    end
    // plot overall square
    x = [0 0 n n 0];
    y = [0 n n 0 0];
    h = myplot(x,y,edge_color);
//    set(h, 'LineWidth', lw*1.2);
//    hold off;
//    axis ij;
    set_axis(0);
    if argn(1)>=1
        Mout = M;
    end
    return;
end

for j=Jmin:Jmax
    qmin = double(~(j==Jmin));
    for q=qmin:3
        [selx,sely] = compute_quadsel(j,q);
        M1 = M(selx,sely,:);
        if strcmp(style, 'abs')
            M1 = rescale( abs(M1) );
            if renormalize
                M1 = 1.5*M1;
                I = find(abs(M1)>1);
                M1(I) = M1(I)./abs(M1(I));
            end
        elseif strcmp(style, 'absinv')
            M1 = rescale(abs(M1));
            if renormalize
                M1 = 2*M1;
                I = find(abs(M1)>1);
                M1(I) = M1(I)./abs(M1(I));
            end
            M1 = rescale( -abs(M1) );
        elseif strcmp(style, 'real')
            if q>0
                st = stdev(M1(:));
                if st<1e-10
                    st=1;
                end
                M1 = 0.3 * M1/st;
                M1 = clamp(M1, -1,1);
                M1 = (M1+1)/2; // rescale( M1 );
            else
                M1 = rescale(M1);
            end
        elseif strcmp(style, 'bw')
            M1 = rescale(abs(M1));
            M1 = double( M1<20/255 );
        elseif strcmp(style, 'centered')
            if q==0
                M1 = rescale(M1);
            else
                M1 = M1 ./ max(abs(M1(:)));
                M1 = sign(M1).*abs(M1).^gam;
                M1 = (M1+1)/2;
            end
        else
            error('unknown style');    
        end
        M(selx,sely,:) = M1;
    end
    // white borders
    if 0 & edge_color>=0
        M(2^j,1:2^(j+1)) = edge_color; 
        M(1:2^(j+1),2^j) = edge_color;
    end
end

// hold on;
s = [1/2 n-1/2];
imageplot(M);

xset('thickness',3);

// imageplot(s,s,M);
// axis image; axis off;
// colormap gray(256);
// display axis separation
for j=Jmax:-1:Jmin;
    x = [0 2^(j+1)];
    y = [2^j 2^j];
    h = myplot(x,y, edge_color);
//    set(h, 'LineWidth', lw);
    y = [0 2^(j+1)];
    x = [2^j 2^j];
    h = myplot(x,y, edge_color);
//    set(h, 'LineWidth', lw);
end
// plot overall square
x = [0 0 n n 0];
y = [0 n n 0 0];
h = myplot(x,y,edge_color);
// set(h, 'LineWidth', lw*1.2);
// hold off;
// axis ij;

set_axis(0);

if argn(1)>=1
    Mout = M;
end

endfunction


function myplot(x,y,str,n)

plot(x,n-y,str);

endfunction

