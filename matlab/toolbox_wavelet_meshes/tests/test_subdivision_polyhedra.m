% Test of subdivision of a base polyhedron.

mesh_types = {'tetra','oct','ico'};
nmesh = length(mesh_types);
nsub = 3;

clf;
for i=1:nmesh
    mesh_type = mesh_types{i};
    [vertex,face] = compute_base_mesh(mesh_type);
    for s=0:nsub
        subplot(nmesh,nsub+1,s+1+(nsub+1)*(i-1));
        plot_mesh(vertex,face);
        lighting flat;
        if s~=nsub     
			options.spherical = 1;
            options.relaxation = 3;
            [vertex,face] = perform_mesh_subdivision(vertex,face,1,options);
        end
    end
end
    