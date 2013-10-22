function plot_spherical_function(vertex,face,f, options)

% plot_spherical_function - display a function on the sphere
%
%   plot_spherical_function(vertex,face,f, options);
%
%   Set:
%       options.use_color=1: add color
%       options.use_elevation=1: extrude the mesh
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
rho = getoptions(options, 'rho', .5);
color = getoptions(options, 'color', 'rescale');

use_color       = getoptions(options, 'use_color', 1);
use_elevation   = getoptions(options, 'use_elevation', 1);

if isempty(f)
    f = vertex(1,:);
end

if size(f,1)<size(f,2)
    f = f';
end
f = f(:,1);

if iscell(vertex)
    vertex = vertex{end};
end
if iscell(face)
    face = face{end};
end

if use_color
    % scale the colors
    switch color
        case 'rescale'
            fv = rescale(f);
        case 'wavelets'
            fv = f;
            if sum(fv<0)>0
                fv(fv<0) =  -fv(fv<0) / min(fv(fv<0));
            end
            fv(fv>0) =  fv(fv>0) / max(fv(fv>0));
        otherwise
            error('Unknown color type.');
    end
else
    fv = f*0+1;
end

if use_elevation
    normals = vertex;
    d = sum( vertex.^2, 1 );
    if std(d)>.01 % we are not on a sphere
         normals = compute_normal(vertex,face);
    end
%    vertex = vertex + normals .* repmat( rescale(f(:)',0,rho), [3 1] );
    vertex = vertex + normals .* repmat( f(:)'/max(f)*rho, [3 1] );
end

% display
options.face_vertex_color = fv;
plot_mesh(vertex, face, options);
colormap jet(256);

func = getoptions(options, 'func', '');
switch func
    case {'image', 'earth'}
        image_name = getoptions(options, 'image_name', 'lena');
        if strcmp(image_name, 'lena')
            view(125,-15);
        end
        if strcmp(func, 'earth')
            view(125,35);
        end
        colormap gray(256);
        lighting none;
    case {'mesh'}
        if use_color
            colormap jet(256);
        else
            colormap gray(256);
        end
end
axis tight;
if strcmp(color, 'wavelets')
    colormap jet(256);
end