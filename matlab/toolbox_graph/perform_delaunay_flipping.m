function [face, flips, flipsinv] = perform_delaunay_flipping(vertex,face,options)

% perform_delaunay_flipping - compute Dalaunay triangulation via flipping
%
%   [face1, flips, flipsinv] = perform_delaunay_flipping(vertex,face,options);
%
%   Set options.display_flips = 1 for graphical display.
%
%   face is turned into a Delaunay triangulation face1 of the points by flipping
%   the edges. 
%
%   'flips' list the flips needed to go from face to face1.
%   'flipsinv' list the flips needed to go from face1 to face.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
display_flips = getoptions(options, 'display_flips', 0);
verb = getoptions(options, 'verb', 1);

% edge topology
edge = compute_edges(face);
e2f = compute_edge_face_ring(face);

% initialize the stack
estack = find( check_incircle_edge(vertex,face,edge)==0 );

A = sort( compute_triangulation_angles(vertex,face) );

count = 0;
flips = []; flipsinv = [];
ntot = length(estack);
vbar = 0;
while not(isempty(estack))
    if verb
        vbar = max(vbar, ntot - length(estack) );
        progressbar( max(vbar,1), ntot);
    end
    count = count+1;
    if count>Inf % ntot
        warning('Problem with Delaunay flipping');
        break;
    end
    % pop non valid edge
    k = estack(1); estack(1) = [];
    i = edge(1,k); j = edge(2,k);
    f1 = e2f(i,j); f2 = e2f(j,i);
    %%% try to flip it %%%
    % find the two other points
    k1 = find_other( face(:,f1), i,j );
    k2 = find_other( face(:,f2), i,j );
    
    % new faces
    t1 = [k1;k2;i];
    t2 = [k1;k2;j];
    flips(:,end+1) = [i;j];
    flipsinv(:,end+1) = [k1;k2];
    face(:,f1) = t1;
    face(:,f2) = t2;
    % update e2f
    edge(:,k) = [k1;k2];
    e2f = compute_edge_face_ring(face);
    % display
    if display_flips
        clf; plot_graph(triangulation2adjacency(face),vertex);
        drawnow; % pause(.01);
    end    
    
    Anew = sort( compute_triangulation_angles(vertex,face) );
    if lexicmp(Anew,A)<=0
        warning('Problem with flipping.');
    end
    A = Anew;
    
    % update the stack
    estack = find( check_incircle_edge(vertex,face,edge)==0 );
    estack = [ estack(estack>k); estack(estack<=k) ];
end
flipsinv = flipsinv(:,end:-1:1);

progressbar(ntot, ntot);

%%
function f = find_other( f, i,j )

I = find(f==i);
if length(I)~=1
    error('Problem');
end
f(I) = [];
I = find(f==j);
if length(I)~=1
    error('Problem');
end
f(I) = [];