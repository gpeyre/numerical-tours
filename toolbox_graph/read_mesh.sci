function [vertex,face] = read_mesh(filename)

// read_mesh - read data from OFF file.
//
//   [vertex,face] = read_off(filename);
//
//   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
//   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
//
//   Copyright (c) 2003 Gabriel Peyré

if isempty(strfind(filename, '.'))
    filename = strcat(['toolbox_graph/' filename '.off']);
end

[vertex,face] = read_off(filename);


endfunction