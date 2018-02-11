function f = perform_wavelet_mesh_transform(vertex,face, f, dir, options)

% perform_wavelet_mesh_transform - compute a wavelet tranform on a mesh
%
%   f = perform_wavelet_mesh_transform(vertex,face, f, dir, options);
%
%   vertex,face must be a semi-regular cell array of meshes.
%
%   Compute the wavelet transform of a function defined on a semi-regular
%   mesh. The full sized mesh is stored as vertex{end},face{end}, and f is
%   a vector with f(i) being the value of the function at vertex i.
%
%   This transform is implemented using the lifting scheme and a butterfly
%   predictor, as explained in
%
%       Peter Schroder and Wim Sweldens
%       Spherical Wavelets: Texture Processing 
%       Rendering Techniques 95, Springer Verlag
%
%       Peter Schrodder and Wim Sweldens
%    	Spherical Wavelets: Efficiently Representing Functions on the Sphere
%       Siggraph 95
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
J = length(vertex);

if size(f,1)<size(f,2)
    f = f';
end

if size(f,2)>1
    for i=1:size(f,2)
        f(:,i) = perform_wavelet_mesh_transform(vertex,face, f(:,i), dir, options);
    end
    return;
end

global vring;
global e2f;
global fring;
global facej;
    
jlist = J-1:-1:1;
if dir==-1
    jlist = jlist(end:-1:1);
end

do_update = getoptions(options,'do_update', 1);
do_scaling = getoptions(options,'do_scaling', 1);

if do_update
    % compute the integral of the scaling function
    n = length(f);
    I = zeros(n,1);
    for j=J-1:-1:1
        % compute navigator helpers
        vring = compute_vertex_ring(face{j+1});
        e2f = compute_edge_face_ring(face{j});
        fring = compute_face_ring(face{j});
        facej = face{j};

        % number of coarse points
        nj = size(vertex{j},2);
        % number of fine points
        nj1 = size(vertex{j+1},2);
        
        if j==J-1
            I(nj+1:end) = 1;
        end

        %%%% PREDICT STEP %%%%
        for k=nj+1:nj1  % for all the fine point
            % retrieve coarse neighbors
            [e,v,g] = compute_butterfly_neighbors(k, nj);
            % do interpolation with butterfly
            I(e) = I(e) + 1/2 * I(k);
            I(v) = I(v) + 1/8 * I(k);
            I(g) = I(g) - 1/16 * I(k);
        end
        
    end    
end


for j=jlist
    % compute navigator helpers
    vring = compute_vertex_ring(face{j+1});
    e2f = compute_edge_face_ring(face{j});
    fring = compute_face_ring(face{j});
    facej = face{j};

    % number of coarse points
    nj = size(vertex{j},2);
    % number of fine points
    nj1 = size(vertex{j+1},2);
    
    %%%% SCALING %%%%%
    if dir==-1 && do_scaling
        f(nj+1:nj1) = f(nj+1:nj1) / sqrt(2^(J-1-j));
    end

    %%%% Reverse UPDATE STEP %%%%
    if dir==-1 && do_update        
        for k=nj+1:nj1  % for all the fine point
            % retrieve coarse neighbors
            e = compute_butterfly_neighbors(k, nj);
            f(e) = f(e) + I(k)./(2*I(e))*f(k);
        end
    end
    %%%% PREDICT STEP %%%%
    for k=nj+1:nj1  % for all the fine point
        % retrieve coarse neighbors
        [e,v,g] = compute_butterfly_neighbors(k, nj);
        % do interpolation with butterfly
        fi = 1/2*sum(f(e)) + 1/8*sum(f(v)) - 1/16*sum(f(g));
        if dir==1
            f(k) = f(k) - fi;
        else
            f(k) = f(k) + fi;
        end
    end
    %%%% Direct UPDATE STEP %%%%
    if dir==1 && do_update
        for k=nj+1:nj1  % for all the fine point
            % retrieve coarse neighbors
            e = compute_butterfly_neighbors(k, nj);
            f(e) = f(e) - I(k)./(2*I(e))*f(k);
        end
    end    
    
    %%%% SCALING %%%%%
    if dir==1 && do_scaling
        f(nj+1:nj1) = f(nj+1:nj1) * sqrt(2^(J-1-j));
    end
end
