function faces = perform_faces_reorientation(vertex,faces, options)

% perform_faces_reorientation - reorient the faces with respect to the center of the mesh
%
%   faces = perform_faces_reorientation(vertex,faces, options);
%
% try to find a consistant reorientation for faces of a mesh.
%
%   if options.method = 'fast', then use a fast non accurate method
%   if options.method = 'slow', then use a slow exact method
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
method = getoptions(options, 'method', 'fast');
verb = getoptions(options, 'verb', 1);

[vertex,faces] = check_face_vertex(vertex,faces);
n = size(vertex,2);
m = size(faces,2);


if strcmp(method, 'fast')
    % compute the center of mass of the mesh
    G = mean(vertex,2);
    % center of faces
    Cf = (vertex(:,faces(1,:)) + vertex(:,faces(2,:)) + vertex(:,faces(3,:)))/3;
    Cf = Cf - repmat(G,[1 m]);
    % normal to the faces
    V1 = vertex(:,faces(2,:))-vertex(:,faces(1,:));
    V2 = vertex(:,faces(3,:))-vertex(:,faces(1,:));
    N = [V1(2,:).*V2(3,:) - V1(3,:).*V2(2,:) ; ...
        -V1(1,:).*V2(3,:) + V1(3,:).*V2(1,:) ; ...
        V1(1,:).*V2(2,:) - V1(2,:).*V2(1,:) ];
    % dot product
    s = sign(sum(N.*Cf));
    % reverse faces
    I = find(s>0);
    faces(:,I) = faces(3:-1:1,I);
    return
end

options.method = 'fast';
faces = perform_faces_reorientation(vertex,faces, options);

fring = compute_face_ring(faces);
tag = zeros(m,1)-1;
heap = 1;
for i=1:m
    if m>100 & verb==1
        progressbar(i,m);
    end
    if isempty(heap)
        I = find(tag==-1); 
        if isempty(I)
            error('Problem');
        end
        heap = I(1);
    end
    f = heap(end); heap(end) = [];
    if tag(f)==1
        warning('Problem');
    end
    tag(f) = 1;    % computed
    fr = fring{f};
    fr = fr(tag(fr)==-1);
    tag(fr)=0;
    heap = [heap fr];
    for k=fr(:)'
        % see if re-orientation is needed
        if check_orientation(faces(:,k), faces(:,f))==1
            faces(1:2,k) = faces(2:-1:1,k);
        end
    end
end

if not(isempty(heap)) || not(isempty(find(tag<1)))
    warning('Problem');
end

% try to see if face are facing in the correct direction
[normal,normalf] = compute_normal(vertex,faces);
a = mean(vertex,2);
a = repmat(a, [1 n]);
dp = sign(sum(normal.*a,1));
if sum(dp>0)<sum(dp<0)
    faces(1:2,:) = faces(2:-1:1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = check_orientation(f1,f2)

i1 = 0; j1 = 0;
for i=1:3
    for j=1:3
        if f1(i)==f2(j)
            i1 = i; j1 = j;
            break;
        end
    end
    if i1~=0
       break;
    end
end
if i1==0
    error('Problem.');
end
ia = mod(i1,3)+1;
ja = mod(j1,3)+1;
ib = mod(i1-2,3)+1;
jb = mod(j1-2,3)+1;
if f1(ia)==f2(ja) || f1(ib)==f2(jb)
    o=1;
else
    o = -1;
end
    
