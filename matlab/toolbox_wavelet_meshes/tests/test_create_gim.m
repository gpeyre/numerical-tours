% Creat a Geometry Image (GIM) from either
%       A geometry image stored in an .m file (H. Hoppe data)
%       A geometry image from a spherical parameterization of a mesh 
%           (for instance gathered from H. Hoppe data).
%   
%   Copyright (c) 2007 Gabriel Peyre

type = 'gim';
type = 'sph';

if not(exist('name'))
    name = 'bunny';
    name = 'venus';
    name = 'armadillo';
    name = 'gargoyle';
    name = 'horse';
end

rep = '../toolbox_graph_data/gim/';
filename = [rep name '-' type '.m'];

disp('--> Loading file.');
[vertex, face, normal, uv, sphparam] = read_mfile(filename, type);

if strcmp(type, 'mesh')
    M = reshape(vertex, [3, sqrt(size(vertex,2)), sqrt(size(vertex,2))]);
    M = shiftdim(M, 1);
    Mnormal = reshape(normal, [3, sqrt(size(normal,2)), sqrt(size(normal,2))]);
    Mnormal = shiftdim(Mnormal, 1);
else
    % build a spherical gim
    n = 256;
    M = perform_sgim_sampling(vertex, sphparam, face, n);
    Mnormal = perform_sgim_sampling(normal, sphparam, face, n);
end

% save geometry image
write_gim(M, [rep name '-' type], 0);
write_gim(Mnormal, [rep name '-' type '.normal'], 0);

plot_geometry_image(M);
colormap gray(256);
