function plot_vf(vf, M, options)

% plot_vf - plot a vector field with 
%   an optional image in the background.
%
% plot_vf(vf, M, options);
%
%   WORKS ONLY FOR 2D VECTOR FIELDS
%
%   set options.display_streamlines=1 to display streamlines.
%   
%   See also: plot_tensor_field.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    M = [];
end

if iscell(vf)
    nrows = 1;
    if length(vf)>3
        nrows = 2;
    end
    if length(vf)>6
        nrows = 3;
    end
    ncols = ceil(length(vf)/nrows);
    lgd = getoptions(options, 'lgd', []);
    clf;
    for i=1:length(vf)
        if iscell(M)
            Mi = M{i};
        else
            Mi = M;
        end
        subplot(nrows,ncols,i);
        plot_vf(vf{i}, Mi, options);
        if iscell(lgd)
            title(lgd{i});
        else
            title(lgd);
        end        
        axis tight;
    end
    return;
end

options.null = 1;
is_oriented = getoptions(options, 'is_oriented', 1);
strech_factor = getoptions(options, 'strech_factor', .6);
reorient = getoptions(options, 'reorient', 0);
linestyle = getoptions(options, 'linestyle', 'b');
display_streamlines = getoptions(options, 'display_streamlines', 0);
streamline_density = getoptions(options, 'streamline_density', 8);
streamline_width = getoptions(options, 'streamline_width', 1);
line_width = getoptions(options, 'line_width', 1);
display_arrows = getoptions(options, 'display_arrows', 1);
subsampling = getoptions(options, 'subsampling', []);
normalize_flow = getoptions(options, 'normalize_flow', 0);

if display_streamlines && ~isfield(options, 'display_arrows')
    display_arrows = 0;
end

if display_arrows==1 && not(isempty(subsampling))
    vf = vf(1:subsampling:end,1:subsampling:end,:);
end

if size(vf,3)~=2
    warning('Dimension >2, cropping ...');
    vf = vf(:,:,1:2);
end

if reorient
    % reorient the vf to x>0
    epsi = sign(vf(:,:,1));
    I = find( epsi==0 );
    epsi(I) = 1;
    vf(:,:,1) = vf(:,:,1).*epsi;
    vf(:,:,2) = vf(:,:,2).*epsi;
end

if normalize_flow
    vf = perform_vf_normalization(vf);
end

n = size(vf,1);
p = size(vf,2);

x = 0:1/(n-1):1;
y = 0:1/(p-1):1;
[Y,X] = meshgrid(y,x);

hold on;

if display_arrows
    imagesc(x,y,M');
    if is_oriented
        h = quiver(X,Y,vf(:,:,1),vf(:,:,2), strech_factor, linestyle);
    else
        h = quiver(X,Y,vf(:,:,1),vf(:,:,2), strech_factor*0.7, linestyle);
        h = quiver(X,Y,-vf(:,:,1),-vf(:,:,2), strech_factor*0.7, linestyle);
    end
    axis xy;
    axis equal;
    set(h, 'LineWidth', line_width);
end


if display_streamlines
    if not(isempty(M))
        imagesc(M');
    end
    [X,Y] = meshgrid(1:n,1:p);
    [XY,tmp] = streamslice(X,Y, vf(:,:,2), vf(:,:,1) ,streamline_density);
    % reverse stream
    for i=1:length(XY)
        XY{i} = XY{i}(:,2:-1:1);
    end
    h = streamline(XY);
    set(h, 'LineWidth', streamline_width);
    axis ij;
end

axis off;
hold off;