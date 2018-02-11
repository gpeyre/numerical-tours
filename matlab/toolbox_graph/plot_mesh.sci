function plot_mesh(vertex,face,options)

// plot_mesh - plot a 3D mesh.
//
//   plot_mesh(vertex,face, options);
//
//   'options' is a structure that may contains:
//       - 'face_vertex_color' : a color per vertex or face.
//
//  Set options.lighting=1 to display with lighting/
//
//   Copyright (c) 2008 Gabriel Peyre

options.null = 0;


[vertex,face] = check_face_vertex(vertex,face);

n = size(vertex,2);
m = size(face,2);


colv = getoptions(options, 'face_vertex_color', []);
colv = colv(:);

if ~isempty(colv)
	// rescale to 0,1
	colv = (colv-min(colv))/(max(colv)-min(colv));
end

if myisfield(options, 'lighting')
	lighting = options.lighting;
	if lighting==1
		// use lighting
		normal = compute_normal(vertex,face);
		n = size(vertex,2);
		if myisfield(options, 'light')
			light = options.light;
		else
			light = [1 1 1];
		end
		light = light/norm(light,'fro');
		c = sum( normal.*mtlb_repmat(light(:),[1 n]), 1);
		c = max(c(:),0);
		if ~isempty(colv)
			colv = colv.*c;
		else
			colv = c;
		end
	end
else
	lighting = 0;
end

// compute coordinates
face = face(3:-1:1,:);
x = matrix(vertex(1,face), 3, m);
y = matrix(vertex(2,face), 3, m);
z = matrix(vertex(3,face), 3, m);

// assign colors
if ~isempty(colv)
	colv = colv*255+1; // rescale to [1,256]
	colors = matrix(colv(face), 3, m);
	Z = list(z,colors);
	f = [-1 2 0];
else
	Z = z;
	f = [0 2 0];
end


if myisfield(options, 'faceted')
if options.faceted==1
	f(1) = abs(f(1));
end
end

t = 0; a = 0; l = '';

// display
// clf();
h=gcf();
if lighting==1
	h.color_map = hotcolormap(256);
else
	h.color_map = graycolormap(256);
end
plot3d(x,y,Z,t,a,l,f);

endfunction