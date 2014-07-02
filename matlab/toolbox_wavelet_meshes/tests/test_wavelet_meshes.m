% test for wavelet on semi regular mesh

path(path, 'toolbox/');
path(path, 'gim/');
path(path, 'data/');

name = 'bunny';
name = 'gargoyle';
name = 'sphere';
name = 'triangle'; 
name = 'rand';

% number of scales
J = 4;
J = 8;

if strcmp(name, 'sphere') || strcmp(name, 'triangle') || strcmp(name, 'rand')
    %% spherical data compression
    options.keep_subdivision = 1;
    if strcmp(name, 'triangle')
        options.base_mesh = 'triangle';
    elseif strcmp(name, 'rand')
        options.base_mesh = 'rand';
        options.nverts_base = 8;
    else
        options.base_mesh = 'tetra';
        options.base_mesh = 'ico';
        options.relaxation = 1;
    end
    disp('--> Loading multiresolution sphere.');
    [vertex,face] = compute_semiregular_sphere(J,options);
    if strcmp(name, 'triangle') || strcmp(name, 'rand')
        f = sum(vertex{end}(1:2,:).^2)';
    else
        %% load a function
        func = 'linear';
        func = 'cos';
        func = 'singular';
        func = 'image';
        func  = 'earth';
        options.func = func;
        disp('--> Loading spherical function');
        f = load_spherical_function(func, vertex{end}, options);
    end
    options.use_elevation = 1;
    options.use_color = 1;
else
    %% semi regular GIM
    M = read_gim([name '-sph.gim']);
    options.func = 'mesh';
    [vertex,face,vertex0] = compute_semiregular_gim(M,J,options);
    %% take the 3 coordinates as functions
    f = vertex{end}';
    options.use_elevation = 0;
    options.use_color = 0;
end

rep = 'results/wavelets-meshes/';
if not(exist(rep))
    mkdir(rep);
end

options.rho = .3;
options.color = 'rescale';
options.use_elevation = 0;
clf;
plot_spherical_function(vertex,face,f, options);
Jmin = 3;
if strcmp(name, 'triangle') || strcmp(name, 'rand')
    Jmin = 1;
end

if not(strcmp(name, 'sphere'))
    for j=Jmin:min(J,6)
        clf;
        if strcmp(name, 'triangle') || strcmp(name, 'rand')
            plot_graph(triangulation2adjacency(face{j}), vertex{j}(1:2,:));
        else
            plot_mesh(vertex{j},face{j},options);
            shading faceted; camlight;
            lighting flat; axis tight;
        end                
        saveas(gcf, [rep name '-semiregular-' num2str(j) '.png'], 'png')
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compress a function or a mesh
if not( strcmp(name, 'triangle') )
disp('--> Computing forward wavelet transform.');
fw = perform_wavelet_mesh_transform(vertex,face, f, +1, options);
% approximation
rlist = [.02 .05 .1 .2 1];
for i = 1:length(rlist)
    if i==length(rlist)
        f1 = f;
    else
        r = rlist(i);
        fwT = keep_biggest( fw, round(r*length(fw)) );
        disp('--> Computing backward wavelet transform.');
        f1 = perform_wavelet_mesh_transform(vertex,face, fwT, -1, options);
    end
  	clf;
    if not(strcmp(name, 'sphere'))
        plot_spherical_function(vertex,face,f1, options);
    else
        plot_mesh(f1,face{j},options);
    end
    if strcmp(name, 'triangle')
        lighting none; axis tight;
    end
    saveas(gcf, [rep name '-compression-' num2str(i) '.png'], 'png');
end
end


if not(strcmp(name, 'triangle'))
    return;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display a wavelet
nverts = size(vertex{end}, 2);
for j = [3 4 5]

    nj = size(vertex{j},2); nj1 = size(vertex{j+1},2);
    sel = nj+1:nj1-1;
    d = sum( abs(vertex{end}(:,sel)) );
    [tmp,k] = min(d); k = sel(k);
    fw2 = zeros(nverts,1); fw2(k) = 1;
    disp('--> Computing backward wavelet transform.');
    f2 = perform_wavelet_mesh_transform(vertex,face, fw2, -1, options);


    options.color = 'wavelets';
    options.use_color = 1;
    options.rho = .4;

    clf;
    options.use_elevation = 1;
    options.view_param = [-5,50];
    plot_spherical_function(-vertex{end},face{end},f2, options);
    saveas(gcf, [rep name '-wavelets-' num2str(j) '.png'], 'png');
end
