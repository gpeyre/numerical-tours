function [vertex,face] = compute_semiregular_sphere(J,options)

% compute_semiregular_sphere - compute a semi-regular sphere
%
%   [vertex,face] = compute_semiregular_sphere(J,options);
%
%   J is the level of subdivision
%
%   set options.keep_subdivision if you want to use this mesh for
%   semi-regular wavelet transform.
%
%   set options.relaxation>0 to enhance the positions
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
base_mesh = getoptions(options, 'base_mesh', 'ico');
keep_subdivision = getoptions(options, 'keep_subdivision', 0);
relaxation = getoptions(options, 'relaxation', 0);

if strcmp(base_mesh, 'triangle') || strcmp(base_mesh, 'rand')
    options.spherical = 0;
    options.sanity_check = 0;
else
    options.spherical = 1;
end
[vertex{1},face{1}] = compute_base_mesh(base_mesh, 0, options);
for j=2:J
    [vertex{j},face{j}] = perform_mesh_subdivision(vertex{j-1},face{j-1}, 1, options);
    if relaxation>0
        % do some smoothing
        options.averaging_type = 'combinatorial';
        for i=1:relaxation
            vertex{end} = perform_mesh_smoothing(face{end},vertex{end},vertex{end},options)';
            d = sqrt( sum(vertex{end}.^2,1) );
            vertex{end} = vertex{end} ./ repmat( d, [3 1]);
        end
    end
end

if not(keep_subdivision)
    vertex = vertex{end};
    face = face{end};
end