function h = plot_surf_texture(M, T)

// plot_surf_texture - plot a surface with a texture on it.
//
//   h = plot_surf_texture(M, T);
//
//   M(:,:,i) for i=1,2,3 are the 3D coordonate of the surface.
//   T is a 2D image mapped on the surface.
//
//   Copyright (c) 2010 Gabriel Peyre


n = size(M,1);
p = size(T,1);



// interpolate
if n~=p
    M1 = M;
    M = zeros(p,p,3);
    for i=1:3
        M(:,:,i) = image_resize(M1(:,:,i),p);
    end
end

colormap(gray(256));

surf(M(:,:,1), M(:,:,2), M(:,:,3), T' );

h=gce();
h.color_flag=3; // interpolated
h.thickness=0;

endfunction