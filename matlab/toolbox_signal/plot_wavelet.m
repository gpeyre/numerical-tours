function Mout = plot_wavelet(M, Jmin, options)

% plot_wavelet - plot 2D wavelet transform stored using Mallat's ordering.
%
%   plot_wavelet(MW, Jmin, options);
%
%   'MW' is the wavelet transform (in Mallat's ordering, not lifting inplace
%       ordering).
%   'Jmin' is the minimum scale of the transform.
%
%   You can set options.style, options.edge_color, options.renormalize
%       options.line_width.
%
%   Works with 1D and 2D wavelets
%   
%   Copyright (c) 2006 Gabriel Peyrée

if nargin<2
    Jmin = 1;
end

M = double(M);

options.null = 0;
style = getoptions(options, 'style', 'real');
edge_color = getoptions(options, 'edge_color', 'r');
renormalize = getoptions(options, 'renormalize', 0);
lw = getoptions(options, 'line_width', 2);
gamma = getoptions(options, 'gamma', .8);

n = size(M,1);
Jmax = log2(n)-1;

if size(M,1)==1 || size(M,2)==1
    % 1D Signal
    hold on;
    h = plot(M); axis tight;
    a = min(M); b = max(M);
    for j=Jmin:Jmax
        plot(2^j*[1 1], [a b], 'r--');
    end
    hold off;
	Mout = [];
    box on;
    return;
end

separable = getoptions(options, 'separable',0);
if separable
    A = M*0;
    for j1=[-1 Jmin:Jmax]
        sel1 = 2^j1+1:2^(j1+1);
        if j1==-1
            sel1 = 1:2^Jmin;
        end
        for j2=[-1 Jmin:Jmax]
            sel2 = 2^j2+1:2^(j2+1);
            if j2==-1
                sel2 = 1:2^Jmin;
            end
            M1 = M(sel1,sel2,:);
            if j1<0 && j2<0
                M1 = rescale(M1);
            else
            vmax = max(abs(M1(:)));
            if vmax<eps
                vmax = 1;
            end
            M1 = M1 / vmax;
            if 1
                st = std(M1(:));
                if st<eps
                    st=1;
                end
                M1 = 0.3 * M1/st;
            end
            M1 = clamp(M1, -1,1);
            M1 = (M1+1)/2;
            end
            M(sel1,sel2,:) = M1; % rescale( M1 );
        end
    end
    hold on;
    s = [1/2 n-1/2];
    imagesc(s,s,M);
    axis image; axis off;
    colormap gray(256);
    % display axis separation
    for j=Jmax:-1:Jmin;        
        x = [0 n];
        y = [2^j 2^j];
        h = plot(x,y, edge_color);
        set(h, 'LineWidth', lw);
        y = [0 n];
        x = [2^j 2^j];
        h = plot(x,y, edge_color);
        set(h, 'LineWidth', lw);
    end
    % plot overall square
    x = [0 0 n n 0];
    y = [0 n n 0 0];
    h = plot(x,y,edge_color);
    set(h, 'LineWidth', lw*1.2);
    hold off;
    axis ij;
    if nargout>=1
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
                st = std(M1(:));
                if st<eps
                    st=1;
                end
                M1 = 0.3 * M1/st;
                M1 = clamp(M1, -1,1);
                M1 = (M1+1)/2; % rescale( M1 );
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
                M1 = M1 ./ mmax(abs(M1));
                M1 = sign(M1).*abs(M1).^gamma;
                M1 = (M1+1)/2;
            end
        else
            error('unknown style');    
        end
        M(selx,sely,:) = M1;
    end
    % white borders
    if 0 && edge_color>=0
        M(2^j,1:2^(j+1)) = edge_color; 
        M(1:2^(j+1),2^j) = edge_color;
    end
end

hold on;
s = [1/2 n-1/2];
imagesc(s,s,M);
axis image; axis off;
colormap gray(256);
% display axis separation
for j=Jmax:-1:Jmin;
    x = [0 2^(j+1)];
    y = [2^j 2^j];
    h = plot(x,y, edge_color);
    set(h, 'LineWidth', lw);
    y = [0 2^(j+1)];
    x = [2^j 2^j];
    h = plot(x,y, edge_color);
    set(h, 'LineWidth', lw);
end
% plot overall square
x = [0 0 n n 0];
y = [0 n n 0 0];
h = plot(x,y,edge_color);
set(h, 'LineWidth', lw*1.2);
hold off;
axis ij;

if nargout>=1
    Mout = M;
end

