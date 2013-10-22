function ic = check_incircle(vertex,face,points)

% check_incicle_edge - compute "empty circle" property for a set of faces
%
%   ic = check_incicle(vertex,points);
% -> ic(i)==1 if the circumcircle of the 3 points 
%   vertex(:,face(:,i)) does not contain the point vertex(:,points(i))
%
%   Copyright (c) 2008 Gabriel Peyre

n = size(face, 2);
if n~=length(points)
    error('Problem');
end

if 0
    %% Apparently there is a bug here in the code %%
    % set of points
    a1 = vertex(:,face(1,:));
    a2 = vertex(:,face(2,:));
    a3 = vertex(:,face(3,:));
    p = vertex(:,points);
    % build matrix for determinant
    if 1
        A = zeros(3,3,n);
        A(1:2,1,:) = reshape(a1-p, [2 1 n]);
        A(1:2,2,:) = reshape(a2-p, [2 1 n]);
        A(1:2,3,:) = reshape(a3-p, [2 1 n]);
        d = [sum((a1-p).^2); sum((a2-p).^2); sum((a3-p).^2)];
        A(3,:,:) = reshape( d, [1 3 n] );
        ic = det3(A)>0;
    else
        A = ones(4,4,n);
        A(1:2,1,:) = reshape(a1, [2 1 n]);
        A(1:2,2,:) = reshape(a2, [2 1 n]);
        A(1:2,3,:) = reshape(a3, [2 1 n]);
        A(1:2,4,:) = reshape(p, [2 1 n]);
        d = [sum(a1.^2); sum(a2.^2); sum(a3.^2); sum(p.^2)];
        A(4,:,:) = reshape( d, [1 4 n] );
        ic = det3(A)<0;
    end
end

% orthocenter
[q,r] = compute_orthocenter(vertex,face);
% distance to centerpoint
P = vertex(:,points);
dc = sum( (q-P).^2, 1);
% check in/out
ic = dc>r.^2;
