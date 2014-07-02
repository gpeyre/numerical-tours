function h = plot_mesh(vertex,face,options)

% plot_mesh - plot a 3D mesh.
%
%   plot_mesh(vertex,face, options);
%
%   'options' is a structure that may contains:
%       - 'normal' : a (nvertx x 3) array specifying the normals at each vertex.
%       - 'edge_color' : a float specifying the color of the edges.
%       - 'face_color' : a float specifying the color of the faces.
%       - 'face_vertex_color' : a color per vertex or face.
%       - 'vertex'
%       - 'texture' : a 2-D image to be mapped on the surface
%       - 'texture_coords' : a (nvertx x 2) array specifying the texture
%           coordinates in [0,1] of the vertices in the texture.
%       - 'tmesh' : set it to 1 if this corresponds to a volumetric tet mesh.
%
%   See also: mesh_previewer.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<2
    error('Not enough arguments.');
end

options.null = 0;

name            = getoptions(options, 'name', '');
normal          = getoptions(options, 'normal', []);
face_color      = getoptions(options, 'face_color', .7);
edge_color      = getoptions(options, 'edge_color', 0);
normal_scaling  = getoptions(options, 'normal_scaling', .8);
sanity_check    = getoptions(options, 'sanity_check', 1);
view_param      = getoptions(options, 'view_param', []);
texture         = getoptions(options, 'texture', []);
texture_coords  = getoptions(options, 'texture_coords', []);
tmesh           = getoptions(options, 'tmesh', 0);

if size(vertex,1)==2
    % 2D triangulation
    % vertex = cat(1,vertex, zeros(1,size(vertex,2)));
    plot_graph(triangulation2adjacency(face),vertex);
    return;
end

% can flip to accept data in correct ordering
[vertex,face] = check_face_vertex(vertex,face);

if size(face,1)==4 && tmesh==1
    %%%% tet mesh %%%%
    % normal to the plane <x,w><=a
    w = getoptions(options, 'cutting_plane', [0.2 0 1]');
    w = w(:)/sqrt(sum(w.^2));
    t = sum(vertex.*repmat(w,[1 size(vertex,2)]));
    a = getoptions(options, 'cutting_offs', median(t(:)) );
    b = getoptions(options, 'cutting_interactive', 0);
    plot_points = getoptions(options, 'plot_points', 0);
    
    while true;

        % in/out
        I = ( t<=a );
        % trim
        e = sum(I(face));
        J = find(e==4);
        facetrim = face(:,J);

        % convert to triangular mesh
        hold on;
        if not(isempty(facetrim))
            face1 = tet2tri(facetrim, vertex, 1);
            % options.method = 'slow';
            face1 = perform_faces_reorientation(vertex,face1, options);
            h{1} = plot_mesh(vertex,face1, options);
        end
        view(3); % camlight;
        shading faceted;
        if plot_points
            K = find(e==0);
            K = face(:,K); K = unique(K(:));
            h{2} = plot3(vertex(1,K), vertex(2,K), vertex(3,K), 'k.');
        end
        hold off;
        
        if b==0
            break;
        end

        [x,y,b] = ginput(1);
        
        if b==1
            a = a+.03;
        elseif b==3
            a = a-.03;
        else
            break;
        end
    end
    return;    
end

vertex = vertex';
face = face';

if strcmp(name, 'bunny') || strcmp(name, 'pieta')
%    vertex = -vertex;
end
if strcmp(name, 'armadillo')
    vertex(:,3) = -vertex(:,3);
end
if sanity_check && ( (size(face,2)~=3 && size(face,2)~=4) || (size(vertex,2)~=3 && size(vertex,2)~=2))
    error('face or vertex does not have correct format.');
end
if ~isfield(options, 'face_vertex_color') || isempty(options.face_vertex_color)
    options.face_vertex_color = zeros(size(vertex,1),1);
end
face_vertex_color = options.face_vertex_color;


if not(isempty(texture))
    %%% textured mesh %%%
    if isempty(texture_coords)
        error('You need to provide texture_coord.');
    end
    if size(texture_coords,2)~=2
        texture_coords = texture_coords';
    end
    opts.EdgeColor = 'none';
    patcht(face,vertex,face,texture_coords,texture',opts);
    if size(texture,3)==1
        colormap gray(256);
    else
        colormap jet(256);
    end
    set_view(name, view_param);
    axis off; axis equal;
    camlight;
    shading faceted;
    return;
end


shading_type = 'interp';
if isempty(face_vertex_color)
    h = patch('vertices',vertex,'faces',face,'facecolor',[face_color face_color face_color],'edgecolor',[edge_color edge_color edge_color]);
else
    nverts = size(vertex,1);
    % vertex_color = rand(nverts,1);
    if size(face_vertex_color,1)==size(vertex,1)
        shading_type = 'interp';
    else
        shading_type = 'flat';
    end
    h = patch('vertices',vertex,'faces',face,'FaceVertexCData',face_vertex_color, 'FaceColor',shading_type);
end
colormap gray(256);
lighting phong;
% camlight infinite; 
camproj('perspective');
axis square; 
axis off;

if ~isempty(normal)
    %%% plot the normals %%%
    n = size(vertex,1);
    subsample_normal = getoptions(options, 'subsample_normal', min(4000/n,1) );
    sel = randperm(n); sel = sel(1:floor(end*subsample_normal));    
    hold on;
    quiver3(vertex(sel,1),vertex(sel,2),vertex(sel,3),normal(1,sel)',normal(2,sel)',normal(3,sel)',normal_scaling);
    hold off;
end

cameramenu;
set_view(name, view_param);
shading(shading_type);
camlight;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_view(name, view_param)


switch lower(name)
    case 'hammerheadtriang'
        view(150,-45);
    case 'horse'
        view(134,-61);
    case 'skull'
        view(21.5,-12);
    case 'mushroom'
        view(160,-75);
    case 'bunny'
%        view(0,-55);
        view(0,90);
    case 'david_head'
        view(-100,10);
    case 'screwdriver'
        view(-10,25);
    case 'pieta'
        view(15,31);
    case 'mannequin'
        view(25,15);
        view(27,6);
    case 'david-low'
        view(40,3);
    case 'david-head'
        view(-150,5);
    case 'brain'
        view(30,40);
    case 'pelvis'
        view(5,-15);
    case 'fandisk'
        view(36,-34);
    case 'earth'
        view(125,35);
    case 'camel'
        view(-123,-5);
        camroll(-90);
    case 'beetle'
        view(-117,-5);
        camroll(-90);
        zoom(.85);
    case 'cat'
        view(-60,15);
    case 'nefertiti'
        view(-20,65);
end

if not(isempty(view_param))
    view(view_param(1),view_param(2));
end

axis tight;
axis equal;

if strcmp(name, 'david50kf') || strcmp(name, 'hand')
    zoom(.85);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function patcht(FF,VV,TF,VT,I,Options)
% This function PATCHT, will show a triangulated mesh like Matlab function
% Patch but then with a texture.
%
% patcht(FF,VV,TF,VT,I,Options);
%
% inputs,
%   FF : Face list 3 x N with vertex indices
%   VV : Vertices 3 x M
%   TF : Texture list 3 x N with texture vertex indices
%   VT : Texture Coordinates s 2 x K, range must be [0..1] or real pixel postions
%   I : The texture-image RGB [O x P x 3] or Grayscale [O x P] 
%   Options : Structure with options for the textured patch such as
%           EdgeColor, EdgeAlpha see help "Surface Properties :: Functions"
%
%   Options.PSize : Special option, defines the image texturesize for each 
%           individual  polygon, a low number gives a more block 
%           like texture, defaults to 64;
%
% note: 
%   On a normal PC displaying 10,000 faces will take about 6 sec.
%
% Example,
%
%  % Load Data;
%   load testdata;
%  % Show the textured patch
%   figure, patcht(FF,VV,TF,VT,I);
%  % Allow Camera Control (with left, right and center mouse button)
%   mouse3d
%
% Function is written by D.Kroon University of Twente (July 2010)

% FaceColor is a texture
Options.FaceColor='texturemap';

% Size of texture image used for every triangle
if(isfield(Options,'PSize'))
    sizep=round(Options.PSize(1));
    Options=rmfield(Options,'PSize');
else
    sizep=64;
end

% Check input sizes
if(size(FF,2)~=size(TF,2))
    error('patcht:inputs','Face list must be equal in size to texture-index list');
end

if((ndims(I)~=2)&&(ndims(I)~=3))
    error('patcht:inputs','No valid Input texture image');
end

% Detect if grayscale or color image
switch(size(I,3))
    case 1
        iscolor=false;
    case 3
        iscolor=true;
    otherwise
        error('patcht:inputs','No valid Input texture image');
end

   
if(max(VT(:))<2)
    % Remap texture coordinates to image coordinates
    VT2(:,1)=(size(I,1)-1)*(VT(:,1))+1;
    VT2(:,2)=(size(I,2)-1)*(VT(:,2))+1;
else
    VT2=VT;
end

% Calculate the texture interpolation values
[lambda1 lambda2 lambda3 jind]=calculateBarycentricInterpolationValues(sizep);
 
% Split texture-image in r,g,b to allow fast 1D index 
Ir=I(:,:,1); if(iscolor), Ig=I(:,:,2); Ib=I(:,:,3); end

% The Patch used for every triangle (rgb)
Jr=zeros([(sizep+1) (sizep+1) 1],class(I));
if(iscolor)
    Jg=zeros([(sizep+1) (sizep+1) 1],class(I));
    Jb=zeros([(sizep+1) (sizep+1) 1],class(I));
end

hold on;
% Loop through all triangles of the mesh
for i=1:size(FF,1)
    % Get current triangle vertices and current texture-vertices
    V=VV(FF(i,:),:); 
    Vt=VT2(TF(i,:),:); 
    
    % Define the triangle as a surface
    x=[V(1,1) V(2,1); V(3,1) V(3,1)];
    y=[V(1,2) V(2,2); V(3,2) V(3,2)];
    z=[V(1,3) V(2,3); V(3,3) V(3,3)];
    
    % Define the texture coordinates of the surface
    tx=[Vt(1,1) Vt(2,1) Vt(3,1) Vt(3,1)];
    ty=[Vt(1,2) Vt(2,2) Vt(3,2) Vt(3,2)] ;
    xy=[tx(1) ty(1); tx(2) ty(2); tx(3) ty(3); tx(3) ty(3)];

    % Calculate texture interpolation coordinates
    pos(:,1)=xy(1,1)*lambda1+xy(2,1)*lambda2+xy(3,1)*lambda3;
    pos(:,2)=xy(1,2)*lambda1+xy(2,2)*lambda2+xy(3,2)*lambda3;
    pos=round(pos); pos=max(pos,1); pos(:,1)=min(pos(:,1),size(I,1)); pos(:,2)=min(pos(:,2),size(I,2));
    posind=(pos(:,1)-1)+(pos(:,2)-1)*size(I,1)+1;
    
    % Map texture to surface image
    Jr(jind)=Ir(posind);
    J(:,:,1)=Jr; 
    if(iscolor)
        Jg(jind)=Ig(posind); 
        Jb(jind)=Ib(posind);
        J(:,:,2)=Jg; 
        J(:,:,3)=Jb;
    end
    
    % Show the surface
    surface(x,y,z,J,Options);
end
hold off; 

function [lambda1 lambda2 lambda3 jind]=calculateBarycentricInterpolationValues(sizep)
% Define a triangle in the upperpart of an square, because only that
% part is used by the surface function
x1=sizep; y1=sizep; x2=sizep; y2=0; x3=0 ;y3=0;
% Calculate the bary centric coordinates (instead of creating a 2D image
% with the interpolation values, we map them directly to an 1D vector)
detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
[x,y]=ndgrid(0:sizep,0:sizep); x=x(:); y=y(:);
lambda1=((y2-y3).*(x-x3)+(x3-x2).*(y-y3))/detT;
lambda2=((y3-y1).*(x-x3)+(x1-x3).*(y-y3))/detT;
lambda3=1-lambda1-lambda2;
% Make from 2D (surface)image indices 1D image indices
[jx jy]=ndgrid(sizep-(0:sizep)+1,sizep-(0:sizep)+1);
jind=(jx(:)-1)+(jy(:)-1)*(sizep+1)+1;



    
 